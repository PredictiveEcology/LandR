utils::globalVariables(c(
  "broadleaf", "conif", "conifer", "coverIndex", "decid", "keep", "LCC", "percTree",
  "sparseness", "totalCoverIndex", "treeRSF", "treeTypeCoverIndex", "updatedLCC"
))

#' Calculate landcover classes based on species cohorts
#'
#' Convert `cohortData` into landcover class based on leading species and density.
#' 1. Classify each pixel as "broadleaf", "mixedwood", or "conifer";
#' 2. Define the openness of each pixel as "dense", "open", or "sparse";
#' 3. Combine these into cover classes for forested cells.
#'
#' By default, landcover classes are coded to match the Land Cover Map of Canada 2005 at 250m
#' (Latifovic & Pouliot, 2005).
#'
#' Code modified from
#' <https://github.com/tati-micheletti/caribouRSF_NT/blob/main/R/makeLCCfromCohortData.R>.
#'
#' @references Latifovic, R. & Pouliot, D. (2005). Multitemporal land cover mapping for Canada:
#' methodology and products. Canadian Journal of Remote Sensing, 31, 347–363.
#'
#' @note Sparse classes:
#' Even though the default LCC2005 input does *not* have "mixed sparse" or "broadleaf sparse",
#' these are likely to be produced in simulations, starting from "conifer sparse" sites
#' (i.e. broadleaf trees might grow on these sites, and eventually convert these into mixed or
#' broadleaf sparse sites).
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @param lccTable `data.table` (with columns `leading` and `LCC`) defining the mappings
#' between leading species and landcover classes.
#'
#' @param deciduousCoverDiscount numeric between 0 and 1 that translates %% cover to %% biomass.
#' It assumes all hardwoods are equivalent; all softwoods are equivalent;
#' and that %% cover of hardwoods will be an overesimate of the %%% biomass of hardwoods.
#' Hardwoods in Canada have a much wider canopy than softwoods.
#' E.g., 30%% cover of hardwoods might translate to 20%% biomass of hardwoods.
#' **NB:** the default `deciduousCoverDiscount` value was estimated from NWT data in March 2020,
#' and may not be useful for other study areas.
#'
#' @template vegLeadingProportion
#'
#' @param decidousSpp character vector of deciduous species names. If not supplied, will
#' attempt to match deciduous species names in [sppEquivalencies_CA] to those in `cohortData`.
#'
#' @param rstLCC landcover raster
#'
#' @return raster of the same type as `pixelGroupMap`, containing LCC values
#'
#' @author Tati Micheletti and Alex Chubaty
#'
#' @export
lccMapGenerator <- function(cohortData, pixelGroupMap,
                            lccTable = NULL,
                            deciduousCoverDiscount = 0.8418911, ## from Biomass_borealDataPrep
                            vegLeadingProportion = 0.75,
                            decidousSpp = NULL,
                            rstLCC) {
  stopifnot(
    !missing(cohortData),
    !missing(pixelGroupMap),
    !missing(rstLCC),
  )

  if (is.null(lccTable)) {
    ## defaults based on LCC2005
    ##
    ## NOTE: there are LCCs that map to multiple cover/density types, and a cover/density type
    ## may map to multiple LCCs; will need to remove duplicates (simplify) downstream.
    ## `keep` defines the rows to use for simplifying. ensure one of each `leading` type is kept.
    lccTable <- rbind(
      c(1L, "conif_dense", TRUE),
      c(2L, "decid_dense", TRUE),
      c(3L, "mixed_dense", TRUE),
      c(4L, "mixed_dense", FALSE),
      c(5L, "mixed_dense", FALSE),
      c(6L, "conif_open", TRUE),
      c(7L, "conif_open", FALSE),
      c(8L, "conif_sparse", TRUE),
      c(9L, "conif_sparse", FALSE),
      c(10L, "conif_sparse", FALSE),

      ## 11 can be either decid open/sparse:
      c(11L, "decid_open", FALSE),
      c(11L, "decid_sparse", TRUE), ## this is the only decid_sparse in the table

      c(12L, "decid_open", TRUE), ## since 11L decid_open not kept, need to keep this one

      ## 13 can be either mixed open/sparse:
      c(13L, "mixed_open", TRUE),
      c(13L, "mixed_sparse", FALSE),

      ## 14 can be either mixed open/sparse:
      c(14L, "mixed_open", FALSE),
      c(14L, "mixed_sparse", TRUE),
      c(15L, "mixed_sparse", FALSE),
      c(32L, "conif_sparse", FALSE)
    ) |>
      as.data.table() |>
      setnames(c("V1", "V2", "V3"), c("LCC", "leading", "keep"))
    lccTable[, LCC := as.integer(LCC)]
    lccTable[, keep := as.logical(keep)]
  }

  if (is.null(decidousSpp)) {
    sppEquiv <- get(data("sppEquivalencies_CA", package = "LandR", envir = environment()),
                    inherits = FALSE)

    decidousSpp <- data.table(
      species = as.character(unique(cohortData[["speciesCode"]])),
      type = equivalentName(unique(cohortData[["speciesCode"]]), sppEquiv, "Type")
    )[type == "Deciduous"][["species"]] |> sort()
  }

  ## 1. Classify each pixel as "broadleaf", "mixedwood", or "conifer" ------------------------------

  ## a. create species table
  species <- as.character(unique(cohortData[["speciesCode"]]))
  treeSpecies <- data.table(
    speciesCode = c(decidousSpp, species[!species %in% decidousSpp]),
    treeRSF = c(
      rep("broadleaf", times = length(decidousSpp)),
      rep("conifer", times = length(species[!species %in% decidousSpp]))
    )
  )

  cohortData <- merge(cohortData, treeSpecies, by = "speciesCode", all.x = TRUE)
  if (NROW(cohortData) == 0) {
    stop(paste(
      "Empty 'cohortData' after merge with tree types.",
      "Ensure 'cohortData$speciesCode' values match those in 'decidousSpp'."
    ))
  }

  ## b. calculate species cover based on biomass
  ## NOTE: in itself, this cover makes no sense but as a percentage, we can use to apply
  ## the threshold to define the leading sp in the pixel correctly
  cohortData[, coverIndex := fifelse(treeRSF == "broadleaf", B / deciduousCoverDiscount, B)]

  ## c. define which ones are conifers and which are broadleaf
  cohortData[, totalCoverIndex := sum(coverIndex), by = c("pixelGroup")]

  ## d. calculate the sum of treeRSF per pixel group
  cohortData[, treeTypeCoverIndex := sum(coverIndex), by = c("pixelGroup", "treeRSF")]

  ## e. calculate the percentage of forest biomass that is
  cohortData[, percTree := treeTypeCoverIndex / totalCoverIndex, by = c("pixelGroup", "treeRSF")]

  ## f. simplify and dcast cohortData to be able to compare the percentages
  cohortDataSim <- unique(cohortData[, c("pixelGroup", "treeRSF", "percTree")])
  cohortDataD <- dcast(
    data = cohortDataSim, formula = pixelGroup ~ treeRSF,
    fill = 0, value.var = "percTree"
  )

  ## g. define first if a stand is pure or mixed:
  cohortDataD[, decid := fifelse(broadleaf >= vegLeadingProportion, 1, 0)]
  cohortDataD[, conif := fifelse(conifer >= vegLeadingProportion, 1, 0)]
  cohortDataD[, leading := colnames(.SD)[max.col(.SD, ties.method = "first")],
    .SDcols = c("decid", "conif")
  ]
  cohortDataD[, leading := fifelse(decid + conif == 0, "mixed", leading)]

  ## h. simplifying and recoding
  cohortDataSim <- unique(cohortDataD[, c("pixelGroup", "leading")])

  ## 2. Define the openness of each pixel as "dense", "open", or "sparse" --------------------------

  sparsenessMap <- rstLCC
  sparsenessMap[!sparsenessMap[] %in% lccTable[["LCC"]]] <- NA
  sparse <- lccTable[["LCC"]][grep(pattern = "sparse", x = lccTable[["leading"]])]
  open <- lccTable[["LCC"]][grep(pattern = "open", x = lccTable[["leading"]])]
  dense <- lccTable[["LCC"]][grep(pattern = "dense", x = lccTable[["leading"]])]
  ## dense = 1; open = 2; sparse = 3
  ## use neg values initially to be able to see changes from LCC values to sparseness
  sparsenessMap[sparsenessMap[] %in% dense] <- -1
  sparsenessMap[sparsenessMap[] %in% open] <- -2
  sparsenessMap[sparsenessMap[] %in% sparse] <- -3
  sparsenessMap <- -sparsenessMap

  sparsenessMapDT <- c(sparsenessMap, pixelGroupMap) |>
    setNames(c("sparseness", "pixelGroup")) |>
    values() |>
    data.table() |>
    na.omit() |>
    unique()
  sparsenessMapDT[, sparseness := as.factor(sparseness)]
  levels(sparsenessMapDT$sparseness) <- c("dense", "open", "sparse")
  sparsenessMapDT[, sparseness := as.character(sparseness)]

  finalDT <- merge(cohortDataSim, sparsenessMapDT, by = "pixelGroup", all.x = TRUE)
  finalDT[, leading := paste(leading, sparseness, sep = "_")]

  ## simplify mapping of LCC and cover/density type
  lccTable <- lccTable[keep == TRUE, ] |> set(NULL, "keep", NULL)

  finalDT <- merge(finalDT, lccTable, by = "leading", all.x = TRUE)

  ## get the new classes to the LCC where they are supposed to be
  newLCC <- rasterizeReduced(
    reduced = finalDT, fullRaster = pixelGroupMap,
    newRasterCols = "LCC", mapcode = "pixelGroup"
  )
  DT <- data.table(pixelID = 1:ncell(newLCC), values(c(rstLCC, newLCC)))
  names(DT) <- c("pixelID", "LCC", "newLCC")
  DT[, updatedLCC := fifelse(!is.na(newLCC), newLCC, LCC)]
  updatedLCCras <- rast(rstLCC)
  set.values(x = updatedLCCras, cells = DT[["pixelID"]], values = DT[["updatedLCC"]])

  return(updatedLCCras)
}
