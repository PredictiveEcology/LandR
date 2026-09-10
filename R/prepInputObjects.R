utils::globalVariables(c(
  "cover", "distYear", "ecoregionGroup", "establishprob", "fireYear",
  "harvestYear", "lcc", "logAge", "longevity", "maxB", "maxANPP", "newAge",
  "postfireregen", "resproutprob", "result", "SCANFIage", "speciesCode"
))

#' Check if all species in have trait values
#'
#' @template speciesLayers
#' @template species
#' @template sppColorVect
#'
#' @return
#' A `list` with the `speciesLayers` and `sppColorVect`
#'   containing only the species that have trait values in `species`
#'
#' @export
checkSpeciesTraits <- function(speciesLayers, species, sppColorVect) {
  missTraits <- setdiff(names(speciesLayers), species$species)
  missTraits <- c(missTraits, setdiff(species$species, species[complete.cases(species), species]))
  if (length(missTraits)) {
    message(cli::col_blue(
      "The following species in 'speciesLayers' have missing traits",
      "and will be excluded:\n",
      paste(missTraits, collapse = " "),
      "\n If this is wrong check if species synonyms are included in 'sppEquiv'"
    ))
    speciesLayers <- speciesLayers[[which(!names(speciesLayers) %in% missTraits)]]
    sppColorVect <- sppColorVect[c(names(speciesLayers), "Mixed")]
  }

  return(list(speciesLayers = speciesLayers, sppColorVect = sppColorVect))
}

#' Make `pixelTable` from biomass, age, land-cover and species cover data
#'
#' @template speciesLayers
#' @template standAgeMap
#' @param ecoregionFiles A list with two objects: the `ecoregionMap` and a table summarizing
#'   its information per `pixelID.` See `ecoregionProducer`.
#' @param biomassMap raster of total stand biomass in t/ha. Biomass units are
#'   converted to g/m^2.
#' @template rasterToMatch
#' @template rstLCC
#' @param printSummary Logical. If `TRUE`, the default, a print out of the
#'   `summary(pixelTable)` will occur.
#' @template doAssertion
#'
#' @return
#' A `data.table` as many rows as non-NA pixels in `rasterToMath` and
#'  the columns containing pixel data from the input raster layers, with
#'  biomass in g/m^2.
#'
#' @export
makePixelTable <- function(
    speciesLayers,
    standAgeMap,
    ecoregionFiles,
    biomassMap,
    rasterToMatch,
    rstLCC,
    # pixelGroupAgeClass = 1,
    printSummary = TRUE,
    doAssertion = getOption("LandR.assertions", TRUE)
) {
  if (missing(rasterToMatch)) {
    rasterToMatch <- rasterRead(speciesLayers[[1]])
    rasterToMatch[] <- 0
    rasterToMatch[is.na(speciesLayers[[1]])] <- NA
  }

  if (missing(ecoregionFiles)) {
    ecoregionFiles <- list()
    ecoregionFiles$ecoregionMap <- rasterRead(rasterToMatch)
    rtmNotNA <- which(!is.na(as.vector(rasterToMatch[])))
    ecoregionFiles$ecoregionMap[rtmNotNA] <- seq_along(rtmNotNA)
    initialEcoregionCodeVals <- as.vector(ecoregionFiles$ecoregionMap[])
  } else {
    initialEcoregionCodeVals <- factorValues2(
      ecoregionFiles$ecoregionMap,
      as.vector(values(ecoregionFiles$ecoregionMap)),
      att = "ecoregion_lcc"
    )
  }

  # message(cli::col_blue("Round age to nearest pixelGroupAgeClass, which is", pixelGroupAgeClass))
  coverMatrix <- matrix(asInteger(speciesLayers[]), ncol = length(names(speciesLayers)))
  colnames(coverMatrix) <- names(speciesLayers)

  ## faster to use as.factor, which is fine for a numeric.
  iec <- if (is.numeric(initialEcoregionCodeVals)) {
    as.factor(initialEcoregionCodeVals)
  } else {
    factor(initialEcoregionCodeVals)
  }
  pixelTable <- data.table(
    initialEcoregionCode = iec,
    cover = coverMatrix,
    pixelIndex = seq(ncell(rasterToMatch)),
    rasterToMatch = as.vector(values(rasterToMatch))
  )
  if (!missing(standAgeMap)) {
    set(pixelTable, NULL, "age", asInteger(as.vector(standAgeMap[])))
    set(pixelTable, NULL, "logAge", .logFloor(as.vector(standAgeMap[])))
  }

  if (!missing(biomassMap)) {
    set(pixelTable, NULL, "totalBiomass", asInteger(as.vector(biomassMap[]) * 100)) # change units)
  }

  if (!missing(rstLCC)) {
    set(pixelTable, NULL, "lcc", as.vector(rstLCC[]))
  }

  # pixelTable <- data.table(#age = asInteger(ceiling(asInteger(as.vector(standAgeMap[])) /
  #                           pixelGroupAgeClass) * pixelGroupAgeClass),
  # logAge = .logFloor(as.vector(standAgeMap[])),
  # initialEcoregionCode = factor(initialEcoregionCodeVals),
  # totalBiomass = asInteger(as.vector(biomassMap[]) * 100), # change units
  # cover = coverMatrix,
  # pixelIndex = seq(ncell(rasterToMatch)),
  # lcc = as.vector(rstLCC[]),
  # rasterToMatch = as.vector(rasterToMatch[]))

  # Remove NAs from pixelTable
  ## 1) If in rasterToMatch
  pixelTable1 <- na.omit(pixelTable, cols = c("rasterToMatch"))
  ## 2) If in rasterToMatch and initialEcoregionCode
  pixelTable2 <- na.omit(pixelTable, cols = c("rasterToMatch", "initialEcoregionCode"))
  ## 3) For species that we have traits for
  coverColNames <- paste0("cover.", names(speciesLayers))
  pixelTable <- na.omit(pixelTable2, cols = c(coverColNames))

  if (NROW(pixelTable1) != NROW(pixelTable)) {
    message(
      "Setting pixels to NA where there is NA in sim$speciesLayers. Vegetation succession",
      " parameters will only be calculated where there is data for species cover.",
      "\n  Check if rasterToMatch shouldn't also only have data where there is cover data,",
      " as this may affect other modules."
    )
  }
  if (NROW(pixelTable2) != NROW(pixelTable)) {
    message("Setting pixels to NA where there is NA in 'ecoregionMap'")
  }

  message(cli::col_blue("rm NAs, leaving", cli::col_magenta(NROW(pixelTable)), "pixels with data"))
  message(cli::col_blue(
    "This is the summary of the input data for age, ecoregionGroup, biomass, speciesLayers:"
  ))
  if (isTRUE(printSummary)) {
    print(summary(pixelTable))
  }

  return(pixelTable)
}

#' Create `speciesEcoregion`
#'
#' Use statistically estimated `maxB`, `maxANPP` and establishment probabilities
#' to generate `speciesEcoregion` table.
#'
#' See Details.
#'
#' @param cohortDataBiomass a subset of `cohortData` object
#' @param cohortDataShort a subset of `cohortData`
#' @param cohortDataShortNoCover a subset of `cohortData`
#' @template species
#' @param modelCover statistical model of species presence/absence
#' @param modelBiomass statistical model of species biomass
#' @template successionTimestep
#' @param currentYear `time(sim)`
#'
#' @section `establishprob`:
#' This section takes the cover as estimated from the mature tree cover and
#' partitions it between resprouting and seeds Unfortunately, establishment by
#' seed is not independent of resprouting, i.e., some pixels would have both
#' Since we don't know the level of independence, we can't correctly assess how
#' much to discount the two. If there is resprouting > 0, then this is the
#' partitioning:
#' `establishprob = f(establishprob + resproutprob + jointEstablishProbResproutProb)`
#' If `jointEstablishProbResproutProb` is 0, then these are independent events
#' and the total cover probability can be partitioned easily between seeds and
#' resprout. This is unlikely ever to be the case. We are picking 50% overlap as
#' a number that is better than 0 (totally independent probabilities, meaning no
#' pixel has both seeds and resprout potential) and  100% overlap (totally
#' dependent probabilities, i.e., every pixel where there is seeds will also be
#' a pixel with resprouting) This is expressed with the "* 0.5" in the code.
#'
#' #' @return
#' A `speciesEcoregion` `data.table` with added columns for parameters
#'   `maxB`, `maxANPP` and `establishprob`
#'
#' @export
makeSpeciesEcoregion <- function(
    cohortDataBiomass,
    cohortDataShort,
    cohortDataShortNoCover,
    species,
    modelCover,
    modelBiomass,
    successionTimestep,
    currentYear
) {
  if (!is.null(modelBiomass$scaledVarsModelB)) {
    if (!is(modelBiomass$scaledVarsModelB, "list")) {
      stop("modelBiomass$scaledVarsModelB must be a list")
    }

    if (!all(names(modelBiomass$scaledVarsModelB) %in% c("cover", "logAge"))) {
      stop("modelBiomass$scaledVarsModelB must be a list with 'cover' and 'logAge' entries")
    }
  }

  ## Create speciesEcoregion table
  joinOn <- c("ecoregionGroup", "speciesCode")
  speciesEcoregion <- unique(cohortDataBiomass, by = joinOn)
  speciesEcoregion[, c("B", "logAge", "cover") := NULL]
  species[, speciesCode := as.factor(species)]
  speciesEcoregion <- species[, .(speciesCode, longevity)][speciesEcoregion, on = "speciesCode"]
  speciesEcoregion[, ecoregionGroup := factor(as.character(ecoregionGroup))]

  ## establishProb ----------------------------------------------------------------------------
  predictedCoverVals <- if (is(modelCover, "numeric")) {
    modelCover
  } else {
    predict(modelCover$mod, newdata = cohortDataShort, type = "response")
  }
  establishprobBySuccessionTimestep <- 1 - (1 - predictedCoverVals)^successionTimestep
  cohortDataShort[, establishprob := establishprobBySuccessionTimestep]
  cohortDataShort <- species[, .(resproutprob, postfireregen, speciesCode)][
    cohortDataShort,
    on = "speciesCode"
  ]

  ## partitioning between seed and resprout. See documentation about the "* 0.5"
  cohortDataShort[, establishprob := pmax(0, pmin(1, (establishprob * (1 - resproutprob * 0.5))))]

  cohortDataShort <- rbindlist(
    list(cohortDataShort, cohortDataShortNoCover),
    use.names = TRUE,
    fill = TRUE
  )
  cohortDataShort[is.na(establishprob), establishprob := 0]

  ## join cohortDataShort with establishprob predictions to speciesEcoregion
  speciesEcoregion <- cohortDataShort[, .(ecoregionGroup, speciesCode, establishprob)][
    speciesEcoregion,
    on = joinOn
  ]

  ## maxB -------------------------------------------------------------------------------------
  ## set age to the age of longevity and cover to 100%
  speciesEcoregion[, `:=`(logAge = .logFloor(longevity), cover = 100)]

  ## rescale if need be (modelBiomass may have been fitted on scaled variables)
  if (!is.null(modelBiomass$scaledVarsModelB)) {
    speciesEcoregion2 <- copy(speciesEcoregion)
    speciesEcoregion2[, `:=`(
      logAge = scale(
        logAge,
        center = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:center"),
        scale = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:scale")
      ),
      cover = scale(
        cover,
        center = attr(modelBiomass$scaledVarsModelB$cover, "scaled:center"),
        scale = attr(modelBiomass$scaledVarsModelB$cover, "scaled:scale")
      )
    )]
    speciesEcoregion2[,
                      maxB := asInteger(predict(modelBiomass$mod, newdata = speciesEcoregion2, type = "response"))
    ]
    speciesEcoregion[, maxB := speciesEcoregion2$maxB]
  } else {
    speciesEcoregion[,
                     maxB := asInteger(predict(modelBiomass$mod, newdata = speciesEcoregion, type = "response"))
    ]
  }

  speciesEcoregion[maxB < 0L, maxB := 0L] # fix negative predictions

  ## maxANPP ----------------------------------------------------------------------------------
  message(cli::col_blue("Add maxANPP to speciesEcoregion -- currently --> maxB/30"))
  speciesEcoregion[, maxANPP := asInteger(maxB / 30)]

  ## clean up unneeded columns
  speciesEcoregion[, `:=`(logAge = NULL, cover = NULL, longevity = NULL, lcc = NULL)]

  speciesEcoregion[, year := currentYear]
  return(speciesEcoregion)
}

#' Create `biomassMap`
#'
#' This is a function that creates the `biomassMap` raster used  for simulations in
#' `Biomass_core` module, using estimated data based on `rawBiomassMap` contained in
#' `pixelCohortData`.
#'
#' @template pixelCohortData
#' @template rasterToMatch
#'
#' @return The `biomassMap`, a raster of total stand biomass per pixel.
#'
#' @export
makeBiomassMap <- function(pixelCohortData, rasterToMatch) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  biomassMap <- rasterRead(rasterToMatch)
  # suppress this message call no non-missing arguments to min;
  # returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(biomassMap[pixelData$pixelIndex] <- pixelData$totalBiomass)

  return(biomassMap)
}

#' Create `minRelativeB` table
#'
#' The table contains expert-based values for minimum relative biomass of each shade tolerance
#' class (the minimum relative biomass a cohort with a given shade tolerance should have to be able
#' to germinate), in each unique ecoregion group.
#' All ecoregion groups currently have the same values.
#'
#' @template pixelCohortData
#'
#' @return a data.frame of min relative biomass values per ecoregion group.
#'
#' @export
makeMinRelativeB <- function(pixelCohortData) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  ## D. Cyr's values result in too many cohorts in more moisture-limited forests of Western Canada.
  ## https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/blob/master/LandisInputs/BSW/biomass-succession-main-inputs_BSW_Baseline.txt
  ##
  ## Adjusted values for western forests:
  minRelativeB <- data.frame(
    ecoregionGroup = as.factor(levels(pixelData$ecoregionGroup)),
    minRelativeBDefaults()
    # X1 = 0.15, ## 0.2
    # X2 = 0.25, ## 0.4
    # X3 = 0.50, ## 0.5
    # X4 = 0.75, ## 0.7
    # X5 = 0.85  ## 0.9
  )

  return(minRelativeB)
}

#' minRelativeB defaults for Western Boreal Forest Canada
#'
#' @export
minRelativeBDefaults <- function() {
  data.frame(X1 = 0.15, X2 = 0.25, X3 = 0.35, X4 = 0.45, X5 = 0.55)
}

#' Create `makePixelGroupMap`
#'
#' Create the `makePixelGroupMap` raster containing `pixelGroups` in `pixelCohortData`.
#'
#' @template pixelCohortData
#' @template rasterToMatch
#'
#' @return a raster with pixel groups
#'
#' @export
makePixelGroupMap <- function(pixelCohortData, rasterToMatch) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  pixelGroupMap <- rasterRead(rasterToMatch)

  ## suppress this message call no non-missing arguments to min;
  ## returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(pixelGroupMap[pixelData$pixelIndex] <- as.integer(pixelData$pixelGroup))

  return(pixelGroupMap)
}

#' Create `standAgeMap`
#'
#' Create the `standAgeMap` raster containing age estimates for `pixelCohortData`.
#' A separate [reproducible::prepInputs()] call will source Canadian National Fire Data Base
#' data to update ages of recently burned pixels. To suppress this, pass NULL/NA `fireURL`

#' @template rasterToMatch
#' @param dataSource Character. One of KNN, NTEMS, or SCANFI.
#'   Defaults to SCANFI for `dataYear` 2020.
#'   Also available:
#'   - KNN for `dataYear` 2011; 2001
#' When `dataSource = "SCANFI"`:
#' - For `dataYear = 2020`, the SCANFI stand age map for 2020 is used directly.
#' - For `dataYear != 2020`, the 2020 SCANFI map is **back-adjusted** to the requested
#'   `dataYear` using NTEMS fire and harvest disturbance layers (1985–2020) and the
#'   2001 kNN stand age map.
#'   Disturbed pixels are assigned ages based on year of disturbance and kNN-based
#'   estimates where applicable, while undisturbed areas are aged by subtracting
#'   the difference between `dataYear` and 2020. Negative ages are set to zero.
#'   This process fills missing or inconsistent SCANFI values and records pixel IDs
#'   that were imputed or adjusted.
#' A separate [reproducible::prepInputs()] call can be used to source Canadian
#' National Fire Database (NFDB) fire polygons, allowing further stand age correction
#' for burned areas. To suppress this, set `fireURL = NULL` or `fireURL = NA`.
#' @param dataYear Numeric. Year for which data is obtained. Can be 2001 or 2011 for KNN, 2020 for SCANFI V1,
#'   or 1985-2025 (5 year intervals) for SCANFI V2.
#' @param dataVersion Character. SCANFI product version for data. Default is currently V2. V1 also available.
#' @param ageURL URL for age map download. Will be supplied based on `dataSource` and `dataYear`
#' @param ageFun passed to 'fun' arg of [reproducible::prepInputs()] of stand age map
#' @param maskWithRTM passed to [reproducible::prepInputs()] of stand age map
#' @param method passed to [reproducible::prepInputs()] for reprojecting the stand age map
#' @param datatype passed to [reproducible::prepInputs()] of stand age map
#' @param writeTo passed to [reproducible::prepInputs()] of stand age map
#' @param firePerimeters fire raster layer fire year values.
#' @param fireURL url to download fire polygons used to update age map. If NULL or NA age
#'   imputation is bypassed. Requires passing `rasterToMatch`. Only used if `firePerimeters`
#'   is missing.
#' @param fireFun passed to [reproducible::prepInputs()] of fire data. Only used if `firePerimeters`
#'   is missing.
#' @param fireField field used to rasterize fire polys. Only used if `firePerimeters`
#'   is missing.
#' @template destinationPath
#' @param ... additional arguments passed to [reproducible::prepInputs()]
#'
#' @return a raster layer stand age map corrected for fires, with an attribute vector of pixel IDs
#'  for which ages were corrected. If no corrections were applied the attribute vector is `integer(0)`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(SpaDES.tools)
#' library(terra)
#' library(reproducible)
#' randomPoly <- randomStudyArea(size = 1e7)
#' randomPoly
#' ras2match <- rast(res = 250, ext = ext(randomPoly), crs = crs(randomPoly))
#' ras2match <- rasterize(randomPoly, ras2match)
#' tempDir <- file.path(tempdir(), "ex_prepInputsStandAgeMap")
#'
#' ## NOT USING FIRE PERIMETERS TO CORRECT STAND AGE
#' ## rasterToMatch does not need to be provided, but can be for masking/cropping.
#' standAge <- prepInputsStandAgeMap(
#'   destinationPath = tempDir,
#'   rasterToMatch = ras2match,
#'   fireURL = NA
#' ) ## or NULL
#' attr(standAge, "imputedPixID")
#'
#' ## USING FIRE PERIMETERS TO CORRECT STAND AGE
#' ## ideally, get the firePerimenters layer first
#' firePerimeters <- Cache(prepInputsFireYear,
#'   url = paste0(
#'     "https://cwfis.cfs.nrcan.gc.ca/downloads",
#'     "/nfdb/fire_poly/current_version/NFDB_poly.zip"
#'   ),
#'   destinationPath = tempDir,
#'   rasterToMatch = ras2match
#' )
#'
#' ## Example adjusting SCANFI to 2000 using NTEMS and kNN
#' standAge2000 <- prepInputsStandAgeMap(
#'   destinationPath = tempdir(),
#'   dataSource = "SCANFI",
#'   dataYear = 2000,
#'   rasterToMatch = rast(res = 250, ext = ext(vect(randomStudyArea(size = 1e7))))
#' )
#' attr(standAge2000, "imputedPixID")
#' }
prepInputsStandAgeMap <- function(
    rasterToMatch = NULL,
    dataSource = "SCANFI",
    dataYear = 2020,
    dataVersion = "V2",
    ageURL = NULL,
    ageFun = "terra::rast",
    maskWithRTM = TRUE,
    method = "bilinear",
    datatype = "INT2U",
    destinationPath = NULL,
    writeTo = NULL,
    firePerimeters = NULL,
    fireURL = paste0(
      "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/",
      "fire_poly/current_version/NFDB_poly.zip"
    ),
    fireFun = "terra::vect",
    fireField = "YEAR",
    ...
) {

  digRTM <- .robustDigest(rasterToMatch)
  dots <- list(...)
  if (is.null(writeTo) && !is.null(dots$filename2)) {
    writeTo <- dots$filename2
  }

  #track pixels that are imputed
  allImputedPixels <- integer(0)

  if (dataSource == "SCANFI") {
    if (dataVersion == "V1") {
      #SCANFI V1 has only one dataYear so age is adjusted using NTEMS disturbance layers
      ageURL <- paste0("https://drive.google.com/file/d/1OdZ7Tznk53KceEyt9dFOBOkxDHEX5X0U")

      standAgeMap <- prepInputs(
        url = ageURL,
        destinationPath = destinationPath,
        datatype = datatype,
        method = method,
        fun = ageFun,
        datatype = datatype,
        to = rasterToMatch,
        ...
      ) |> Cache(.functionName = "prepInputs_ageMapFromSCANFI", omitArgs = "to",
                 .cacheExtra = digRTM)

      if (dataYear != "2020") {
        #use NTEMS to identify disturbances (harvest and fire) that occurred between dataYear and 2020
        #example pixel disturbed in 2002 with dataYear 2000 - use kNN 2001 estimate minus one - if negative, set to 0
        #example pixel disturbed in 1995 with dataYear 2000 -> set stand age to dataYear - YOD = 5
        #example pixel undisturbed with dataYear 2000 -> subtract 20 from SCANFI estimate - if negative, set to 0
        #example pixel disturbed in 2001 with dataYear 2000 - kNN will not have correct age, so set standAge to 16
        #As the largest observed disturbance occurred in 2001 and the time series begins in 1985,
        # the minimum age in 2000 would be 2000 - 1985 + 1
        message(
          "SCANFI data is currently available for 2020 only - age will be adjusted to ",
          dataYear,
          " using various data sources"
        )
        # Download and align NTEMS fire and harvest disturbance layers
        fire_NTEMS <- prepInputs(
          url = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip",
          destinationPath = destinationPath,
          to = standAgeMap,
          method = "near"
        )  |> Cache(.functionName = "prepInputs_CA_ForestFire1985to2020")
        harvest_NTEMS <- prepInputs(
          url = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Harvest_1985-2020.zip",
          destinationPath = destinationPath,
          to = standAgeMap,
          method = "near"
        ) |> Cache(.functionName = "prepInputs_CA_Harvest1985to2020")
        NAflag(fire_NTEMS) <- 0
        NAflag(harvest_NTEMS) <- 0

        baseKNN <- prepInputs(
          url = paste0(
            "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-",
            "attributes_attributs-forests-canada/2001-attributes_attributs-2001/",
            "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif"
          ),
          destinationPath = destinationPath,
          to = standAgeMap,
          method = "near"
        )
        newVals <- data.table::data.table(
          fireYear = as.vector(fire_NTEMS),
          harvestYear = as.vector(harvest_NTEMS),
          pixelID = 1:ncell(standAgeMap)
        )
        newVals <- newVals[
          !is.na(fireYear) | !is.na(harvestYear),
          .(distYear = min(fireYear, harvestYear, na.rm = TRUE)),
          .(pixelID)
        ]
        kNN_AgeAdj <- dataYear - 2001 #if dataYear is 2000, subtract one - if 2010, add nine

        newVals[distYear >= dataYear, newAge := baseKNN[pixelID] + kNN_AgeAdj]
        newVals[distYear < dataYear, newAge := dataYear - distYear]
        newVals[, SCANFIage := standAgeMap[pixelID]]
        if (dataYear == 2000) {
          newVals[distYear == 2001, newAge := 16] #assume these stands were at least 16 (1985 start date of TS)
        }
        #final safety catches - likely disagreement over what is forest
        newVals <- newVals[c(!is.na(newAge) & !is.na(SCANFIage))] #ie kNN and SCANFI agree non-forest
        SCANFI_AgeAdj <- dataYear - 2020
        newStandAgeMap <- standAgeMap + SCANFI_AgeAdj
        newStandAgeMap[newVals$pixelID] <- newVals$newAge
        #some zeroes remain -
        newStandAgeMap[newStandAgeMap < 0] <- 0
        standAgeMap <- newStandAgeMap
        allImputedPixels <- c(newVals$pixelID)
        rm(newStandAgeMap, baseKNN, harvest_NTEMS, fire_NTEMS)
      }
      #if baseYear is 2020, proceed with 2020 standAgeMap
    }
    if (dataVersion == "V2") {
      if (dataYear == "1985") {
        ageURL <- paste0("https://drive.google.com/file/d/1KoGbzQZB-wR8LkngAKTZGShDR0ILBK5r")
      } else if (dataYear == "1990") {
        ageURL <- paste0("https://drive.google.com/file/d/1083izY1R-uBRdnnJL6KTYTPg3RFKCU6V")
      } else if (dataYear == "1995") {
        ageURL <- paste0("https://drive.google.com/file/d/12aqGD4YQKrZXUrzSPIOz2DcVY00mQZAT")
      } else if (dataYear == "2000") {
        ageURL <- paste0("https://drive.google.com/file/d/1v-9sx1a_-WqKfuULT80QZ1b0CAI8lLGV")
      } else if (dataYear == "2005") {
        ageURL <- paste0("https://drive.google.com/file/d/1HSWZ7aHW9GxfYTI3aAZRz5i2rD9-pymq")
      } else if (dataYear == "2010") {
        ageURL <- paste0("https://drive.google.com/file/d/1PIC0pvDUZFx7DfauJvB2WsJjgBgkCaXb")
      } else if (dataYear == "2015") {
        ageURL <- paste0("https://drive.google.com/file/d/1_ZGdjepqS3tGHKykAG5SCecic4Qn-Nx-")
      } else if (dataYear == "2020") {
        ageURL <- paste0("https://drive.google.com/file/d/1nXPS3bpFUESYieNfXO25OKlZJEgqtRnD")
      } else if (dataYear == "2025") {
        ageURL <- paste0("https://drive.google.com/file/d/1mM2z-_sjt9HGv1JQIt2VZkhW_HrldhgT")
      } else {
        stop("SCANFI V2 data is currently available for 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, and 2025 only")
      }

      standAgeMap <- prepInputs(
        url = ageURL,
        destinationPath = destinationPath,
        datatype = datatype,
        method = method,
        fun = ageFun,
        # datatype = datatype,
        to = rasterToMatch,
        ...
      ) |> Cache(.functionName = paste0("prepInputs_ageMapFromSCANFI", "_", dataYear), 
                 omitArgs = "to",
                 .cacheExtra = digRTM)

    }
  }
  else {
    if (is.null(ageURL)) {
      if (dataSource == "KNN") {
        if (dataYear == "2011") {
          ageURL <- paste0(
            "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
            "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
            "NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif"
          )
        } else if (dataYear == "2001") {
          ageURL <- paste0(
            "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
            "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/",
            "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif"
          )
        } else {
          stop("KNN data is available for 2001 or 2011 only")
        }
      } else {
        stop("unrecognized dataSource")
      }
    }

    if (is.null(rasterToMatch)) {
      maskWithRTM <- FALSE
    }

    standAgeMap <- Cache(
      prepInputs,
      ...,
      maskWithRTM = maskWithRTM,
      method = method,
      datatype = datatype,
      destinationPath = destinationPath,
      url = ageURL,
      rasterToMatch = rasterToMatch
    )
  }

  #get NFDB fires
  getFires <- if (
    is.null(firePerimeters) && (isFALSE(is.null(fireURL)) && isFALSE(is.na(fireURL)))
  ) {
    TRUE
  } else {
    FALSE
  }

  if (is(standAgeMap, "SpatRaster")) {
    standAgeMap <- as.int(standAgeMap + 0.5)
    # vals <- as.vector(standAgeMap[])
  } else {
    vals <- standAgeMap[]
    standAgeMap[] <- asInteger(vals)
  }


  if (getFires) {
    if (isFALSE(is.null(rasterToMatch))) {
      firePerimeters <- prepInputsFireYear(
        ...,
        url = fireURL,
        fun = fireFun,
        fireField = fireField,
        destinationPath = destinationPath,
        rasterToMatch = rasterToMatch
      ) |> Cache(omitArgs = "rasterToMatch",
                 .cacheExtra = digRTM)
    } else {
      message(
        "No 'rasterToMatch' or 'firePerimeters' supplied; ages will NOT be adjusted using fire data."
      )
    }
  }

  if (isFALSE(is.null(firePerimeters))) {
    standAgeMap <- replaceAgeInFires(standAgeMap, firePerimeters, startTime = dataYear)
    imputedPixID <- attr(standAgeMap, "imputedPixID")
    #From SCANFI
    allImputedPixels <- unique(c(allImputedPixels, imputedPixID))
  }

  if (!is.null(writeTo)) {
    standAgeMap <- writeTo(standAgeMap, writeTo, ...)
  }

  attr(standAgeMap, "imputedPixID") <- allImputedPixels
  return(standAgeMap)
}

#' Create `rawBiomassMap`
#'
#' Create the `rawBiomassMap` raster containing biomass estimates for `pixelCohortData`.
#'
#'
#' @param dataSource Character. One of KNN, NTEMS, or SCANFI.
#'   Defaults to KNN for `dataYear` 2001.
#'   Also available:
#'   - KNN for `dataYear` 2011;
#'   - NTEMS for `dataYear` 2015;
#'   - SCANFI for `dataYear` 2000, 2010, or 2020.
#'
#' @param dataYear Numeric. Year for which data is obtained. Can be 2001 or 2011 for KNN, 2015 for NTEMS,
#'    2000, 2010, or 2020 for SCANFI V1, or 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020 (default), 2025 possible for V2.
#'
#' @param dataVersion Character. SCANFI product version for data. Default is currently V2. V1 also available.
#'
#' @param ... arguments passed to [reproducible::prepInputs()] and [reproducible::Cache()].
#' If the following arguments are not provided, the following values will be used:
#'   \itemize{
#'     \item{`url`: by default, the 2020 SCANFI stand biomass map is downloaded from
#'       a private google drive }
#'     \item{`useSAcrs` and `projectTo`: `FALSE` and `NA`}
#'     \item{`method`: `"bilinear"`}
#'     \item{`datatype`: `"INT2U"`}
#'     \item{`overwrite`: `TRUE`}
#'     \item{`omitArgs`: `c("destinationPath", "targetFile", "userTags", "stable")`}
#'   }
#'
#' @return a `rawBiomassMap` raster
#'
#' @export
prepRawBiomassMap <- function(dataSource = "SCANFI", dataYear = "2020", dataVersion = "V2", ...) {
  Args <- list(...)

  if (!(dataSource %in% c("KNN", "NTEMS", "SCANFI"))) {
    stop("Data Source must be either KNN, NTEMS, or SCANFI")
  }
  if (is.null(Args$url)) {
    if (dataSource == "KNN") {
      if (dataYear == "2011") {
        Args$url <- paste0(
          "http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
          "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
          "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
        )
      } else if (dataYear == "2001") {
        Args$url <- paste0(
          "http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
          "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/",
          "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
        )
      } else {
        stop("KNN data is available for 2001 or 2011 only")
      }
    } else if (dataSource == "NTEMS") {
      if (dataYear == "2015") {
        Args$url <- paste0(
          "https://drive.google.com/file/d/19R4IXxByGvG3V3oE6VjhYnwqTQjQGVC-/view?usp=drive_link"
        )
      } else {
        stop("NTEMS data is currently available for 2015 only")
      }
    } else if (dataSource == "SCANFI") {
      if (dataVersion == "V1") {
        if (dataYear == "2000") {
          Args$url <- paste0("https://drive.google.com/file/d/1lubpotPt-Tr_x1PHnLP6YL36fGg5Ic6h")
        } else if (dataYear == "2010") {
          Args$url <- paste0("https://drive.google.com/file/d/1J3izr9d0IaUs0H4GWJNbn6rCan-Or7Jf")
        } else if (dataYear == "2020") {
          Args$url <- paste0("https://drive.google.com/file/d/1lexPzmm4zeY_5nljoNmsIlzrZYd1TpG_")
        } else {
          stop("SCANFI V1 data is currently available for 2000, 2010, and 2020 only")
        }
      } else if (dataVersion == "V2") {
        if (dataYear == "1985") {
          Args$url <- paste0("https://drive.google.com/file/d/1rSbEs9PpjS5D5n4R2QL0UYVP1xGepRg5")
        } else if (dataYear == "1990") {
          Args$url <- paste0("https://drive.google.com/file/d/15psX9X_ElAxg3oZqfp9QOpy1b_i8OO-4")
        } else if (dataYear == "1995") {
          Args$url <- paste0("https://drive.google.com/file/d/1KdR4k9Bb95-Y2pFCTBCkuC9yIepE8bK6")
        } else if (dataYear == "2000") {
          Args$url <- paste0("https://drive.google.com/file/d/1GJMLSZweBW3dngDf3RRs4lTX-eMMr7_6")
        } else if (dataYear == "2005") {
          Args$url <- paste0("https://drive.google.com/file/d/1NxgZXKPiFTWRTHo7b40jCJqLvvvMp1tY")
        } else if (dataYear == "2010") {
          Args$url <- paste0("https://drive.google.com/file/d/1HFbXmH6o_2zXezEC6wQKRlPWcwjWl-VX")
        } else if (dataYear == "2015") {
          Args$url <- paste0("https://drive.google.com/file/d/1aYzXALVkOvW18CgRXBnmhqdvWtqOo-7Q")
        } else if (dataYear == "2020") {
          Args$url <- paste0("https://drive.google.com/file/d/13-atqi_7ogRPIFxOoJZoUDYdQCJ5-a_u")
        } else if (dataYear == "2025") {
          Args$url <- paste0("https://drive.google.com/file/d/12MFxY0F9go8zDXNpNsx8UvQUo4cnTdE0")
        } else {
          stop("SCANFI V2 data is currently available for 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, and 2025 only")
        }
      }
    }
  }
  ## NOTE: only calling httr2::request here because listed in Imports, to satisfy R CMD check;
  ##       httr is actually needed for reproducible::prepInputs() but it's only a Suggests there;
  ##       see LandR#113 and discussion therein
  url <- httr2::request(Args$url)$url

  if (is.null(Args$writeTo)) {
    if (!is.null(Args$filename2)) {
      Args$writeTo <- Args$filename2
      Args$filename2 <- NULL
    }
  }

  if (is.null(Args$overwrite)) {
    ## when prepInputs below fails for some reason, often the file gets downloaded, but it is corrupted
    ##   If it didn't fail, then the `Cache` will work and not trigger a new prepInputs, so it is safe
    ##   and won't redownload
    Args$overwrite <- TRUE
  }

  Args2 <- list()
  if (is.null(Args$omitArgs)) {
    Args2$omitArgs <- c("destinationPath", "targetFile", "stable")
  }

  if (is.null(Args2$quick)) {
    Args2$quick <- c("writeTo")
  }

  rawBiomassMap <- do.call(prepInputs, args = Args) |>
    Cache(quick = Args2$quick, .functionName = "prepInputsRawBiomassMap", omitArgs = Args2$omitArgs)

  return(rawBiomassMap)
}

#' Create a raster of fire perimeters
#'
#' @param ... Additional arguments passed to [reproducible::prepInputs()]
#' @template rasterToMatch
#' @param fireField field used to rasterize fire polys
#' @param earliestYear the earliest fire date to allow
#'
#' @return a `SpatRaster` layer of fire perimeters with fire year values.
#'
#' @export
#'
#' @examplesIf !isTRUE(as.logical(Sys.getenv("CI")))
#' ## NOTE: runs locally, skipped on CI -- downloads NFDB_poly.zip from
#' ## cwfis.cfs.nrcan.gc.ca (same predicate as testthat::skip_on_ci())
#' withr::local_options(list(
#'   reproducible.useTerra = TRUE,
#'   reproducible.rasterRead = "terra::rast"
#' ))
#'
#' randomPoly <- LandR::randomStudyArea(
#'   size = 1e+8, seed = 5
#' )
#'
#' ras2match <- terra::rast(
#'   randomPoly,
#'   vals = 1,
#'   res = 100,
#' )
#' ras2match <- terra::mask(ras2match, randomPoly)
#'
#' firePerimeters <- prepInputsFireYear(
#'   url = paste0(
#'     "https://cwfis.cfs.nrcan.gc.ca/downloads/",
#'     "nfdb/fire_poly/current_version/NFDB_poly.zip"
#'   ),
#'   destinationPath = file.path(tempdir(), "ex_prepInputsFireYear"),
#'   rasterToMatch = ras2match,
#'   earliestYear = 1950
#' )
#'
#' if (interactive()) {
#'   terra::plot(firePerimeters)
#'   terra::plot(randomPoly, add = TRUE)
#' }
#'
#' withr::deferred_run()
#'
prepInputsFireYear <- function(..., rasterToMatch, fireField = "YEAR", earliestYear = 1950) {
  dots <- list(...)
  fun <- if (is.null(dots$fun)) "terra::vect" else dots$fun
  dots$fun <- NULL # need to do this or else it will pass double to the prepInputs

  ## invalid NFDB polygons will cause Rstudio to crash during postProcess as of 8/21/2024
  ## removing invalid polygons is far faster than fixing the 0.1% of data
  ## projectTo must be rasterToMatch due to terra rasterize, but don't project yet because of NFDB
  postProcessArgs <- dots[names(dots) %in% c("to", "projectTo", "studyArea", "maskTo")]
  if (length(postProcessArgs) == 0) {
    postProcessArgs$cropTo <- rasterToMatch
    postProcessArgs$projectTo <- rasterToMatch
    postProcessArgs$maskTo <- rasterToMatch
  }
  postProcessArgs$projectTo <- rasterToMatch

  preProcessArgs <- dots[!names(dots) %in% names(postProcessArgs)]
  ## you can crop without worrying about geometry
  preProcessArgs$cropTo <- rasterToMatch

  ## Load polygons
  files <- do.call(preProcess, append(list(fun = fun), preProcessArgs))
  files2 <- files$checkSums[result %in% "OK"]$actualFile
  shpFiles <- grep(files2, pattern = ".shp$", value = TRUE)

  preProcessArgs2 <- preProcessArgs
  preProcessArgs2$url <- NULL

  ## There are at least 2 .shp files now (as of Dec 2, 2025)
  vv <- Map(shp = shpFiles, function(shp) terra::vect(file.path(preProcessArgs2$destinationPath, shp)))
  vvv <- terra::vect(vv)
  # 
  # unique(vvv$SRC_AGENCY)
  # sa <- setupStudyArea(list(NAME_1 = c("Alberta", "British Columbia", "Yukon", "Northwest Territories") |> 
  #                             paste(collapse = "|"))) |> 
  #   terra::project(terra::crs(rasterToMatch))
  # unique(vvv$SRC_AGENCY)
  # 
  
  allFires <- do.call(postProcess, append(list(vvv), preProcessArgs2))
  # lots <- Map(shp = shpFiles, function(shp) {
  #   preProcessArgs2$targetFile = file.path(preProcessArgs2$destinationPath, shp)
  #   shp1 <- terra::vect(preProcessArgs2$targetFile)
  #   preProcessArgs2$targetFile <- NULL
  #   do.call(postProcess, append(list(x = shp1), preProcessArgs2))
  # })
  # 
  # allFires <- lots[[1]]
  # 
  # if (length(lots) > 1) {
  #   for (i in 2:length(lots)) {
  #     allFires <- rbind(allFires, lots[[i]])
  #   }
  # }

  # allFires <- do.call(prepInputs, append(list(fun = fun), preProcessArgs))

  ## the reason this isn't combined into one function is due to geometry issues in NFDB
  allFires <- allFires[terra::is.valid(allFires), ] ## drop invalid geometries

  ## If no valid polygons, return empty raster
  if (nrow(allFires) == 0) {
    if (inherits(rasterToMatch, "SpatRaster")) {
      fireRas <- rast(rasterToMatch, vals = NA)
    } else {
      fireRas <- raster::raster(rasterToMatch)
      fireRas[] <- NA
    }
    return(fireRas)
  }

  ## Transform to raster CRS if needed
  if (!identical(crs(allFires), crs(rasterToMatch))) {
    allFires <- terra::project(allFires, crs(rasterToMatch))
  }

  if (isTRUE(grepl("vect", fun))) {
    allFires <- st_as_sf(allFires)
  }

  allFires <- st_zm(allFires)

  allFires <- st_cast(allFires, "MULTIPOLYGON") ## collapse them into a single multipolygon
  allFires <- st_transform(allFires, crs(rasterToMatch))
  if (!is(allFires[[fireField]], "numeric")) {
    warning("Chosen fireField will be coerced to numeric")
    d[[fireField]] <- as.numeric(as.factor(d[[fireField]]))
  }
  if (is(rasterToMatch, "SpatRaster")) {
    if (!is(allFires, "SpatVector")) {
      allFires <- vect(allFires)
    }

    ## fun = max to take the most recent fire year
    fireRas <- terra::rasterize(allFires, rasterToMatch, field = fireField, fun = max)
    fireRas[
      !is.na(terra::values(fireRas, mat = FALSE)) &
        terra::values(fireRas, mat = FALSE) < earliestYear
    ] <- NA
  } else {
    .requireNamespace("fasterize", stopOnFALSE = TRUE)
    fireRas <- fasterize::fasterize(d, raster = rasterToMatch, field = fireField)
    fireRas[!is.na(as.vector(fireRas[])) & as.vector(fireRas[]) < earliestYear] <- NA
  }

  ## This may potentially result in dots intended for postProcess being lost.
  fireRas <- do.call(postProcess, append(list(x = fireRas), postProcessArgs))

  return(fireRas)
}

#' Replace stand age with time since last fire
#'
#' @param standAgeMap a raster layer stand age map
#' @param firePerimeters the earliest fire date to allow
#' @template startTime
#'
#' @return a raster layer stand age map corrected for fires, with an attribute vector of pixel IDs
#'  for which ages were corrected. If no corrections were applied the attribute vector is `integer(0)`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' randomPoly <- terra::vect(SpaDES.tools::randomStudyArea(size = 1e7))
#' randomPoly
#' ras2match <- terra::rast(resolution = 250,
#'                          ext = terra::ext(randomPoly),
#'                          crs = terra::crs(randomPoly))
#' ras2match <- terra::rasterize(randomPoly, ras2match)
#' tempDir <- file.path(tempdir(), "ex_replaceAgeInFires")
#'
#' standAge <- reproducible::prepInputsStandAgeMap(
#'   destinationPath = tempDir,
#'   rasterToMatch = ras2match,
#'   fireURL = NA
#' ) ## or NULL
#' attr(standAge, "imputedPixID")
#'
#' firePerimeters <- reproducible::Cache(prepInputsFireYear,
#'   url = paste0(
#'     "https://cwfis.cfs.nrcan.gc.ca/downloads",
#'     "/nfdb/fire_poly/current_version/NFDB_poly.zip"
#'   ),
#'   destinationPath = tempDir,
#'   rasterToMatch = ras2match
#' )
#' standAge <- replaceAgeInFires(standAge, firePerimeters)
#' attr(standAge, "imputedPixID")
#' }
replaceAgeInFires <- function(standAgeMap, firePerimeters, startTime) {
  if (missing(startTime)) {
    startTime <- 0
  }
  if (startTime < 1950 || startTime > 2023) {
    message("'startTime' is missing or is not within a reasonable range of 1950 to 2023, ")
    message("  --> The most recent fire year will be used.")
    startTime <- max(firePerimeters[], na.rm = TRUE)
  }

  toChange <- !is.na(as.vector(firePerimeters[])) &
    as.vector(firePerimeters[]) <= asInteger(startTime)
  standAgeMap[] <- asInteger(as.vector(standAgeMap[]))
  standAgeMap[toChange] <- asInteger(startTime) - asInteger(firePerimeters[][toChange])
  imputedPixID <- which(toChange)

  attr(standAgeMap, "imputedPixID") <- imputedPixID
  return(standAgeMap)
}

#' Create `rasterToMatch` and `rasterToMatchLarge`
#'
#' `rasterToMatch` and `rasterToMatchLarge` raster layers are created
#'   from `studyArea` and `studyAreaLarge` polygons (respectively)
#'   using a template raster (often `rawBiomassMap`)
#'
#' @template studyArea
#' @param studyAreaLarge same as `studyArea`, but larger and completely
#'   covering it.
#' @template rasterToMatch
#' @template rasterToMatchLarge
#' @template destinationPath
#' @param templateRas a template raster used to make `rasterToMatch`
#'   and/or `rasterToMatchLarge`. Must match `studyAreaLarge`.
#' @template studyAreaName
#' @template cacheTags
#'
#' @export
prepRasterToMatch <- function(
    studyArea,
    studyAreaLarge,
    rasterToMatch,
    rasterToMatchLarge,
    destinationPath,
    templateRas,
    studyAreaName,
    cacheTags = NULL
) {
  if (is.null(rasterToMatch) || is.null(rasterToMatchLarge)) {
    ## if we need rasterToMatch/rasterToMatchLarge, that means a) we don't have it,
    ## but b) we will have templateRas

    if (is.null(rasterToMatchLarge) && !is.null(rasterToMatch)) {
      rasterToMatchLarge <- rasterToMatch
    } else if (is.null(rasterToMatchLarge) && is.null(rasterToMatch)) {
      warning(paste0(
        "rasterToMatch and rasterToMatchLarge are missing. Both will be created \n",
        "from templateRas and studyArea/studyAreaLarge.\n
                     If this is wrong, provide both rasters"
      ))

      if (is.null(templateRas)) {
        stop(paste(
          "Please provide a template raster to make rasterToMatch(Large).",
          "An option is to use 'rawBiomassMap'"
        ))
      }
      if (!.compareRas(templateRas, studyAreaLarge, stopOnError = FALSE)) {
        ## note that extents/origin may never align if the resolution and projection do not allow for it
        templateRas <- Cache(
          postProcessTo,
          templateRas,
          cropTo = studyAreaLarge,
          maskTo = studyAreaLarge,
          # studyArea = studyAreaLarge,
          # useSAcrs = FALSE,
          overwrite = TRUE,
          userTags = c("postRTMtemplate")
        )
        templateRas <- fixErrors(templateRas)
      }
      rasterToMatchLarge <- templateRas
    }

    if (!anyNA(as.vector(rasterToMatchLarge[]))) {
      whZeros <- as.vector(rasterToMatchLarge[]) == 0
      if (sum(whZeros) > 0) {
        ## means there are zeros instead of NAs for RTML --> change
        rasterToMatchLarge[whZeros] <- NA
        message(
          "There were no NAs on the rasterToMatchLarge, but there were zeros;",
          " converting these zeros to NA."
        )
      }
    }

    RTMvals <- as.vector(rasterToMatchLarge[])
    rasterToMatchLarge[!is.na(RTMvals)] <- 1 # converts to RAM object

    ## use try -- overwrite can fail with terra + windows if raster was loaded
    ## previously by another module; unlinking/deleting the file does not work in this case.
    rasterToMatchLargeTmp <- try(Cache(
      writeOutputs,
      rasterToMatchLarge,
      datatype = "INT2U",
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatchLarge"),
      omitArgs = c("userTags")
    ))
    if (!is(rasterToMatchLargeTmp, "try-error")) {
      rasterToMatchLarge <- rasterToMatchLargeTmp
    }
    if (is.null(rasterToMatch)) {
      rtmFilename <- .suffix(
        file.path(destinationPath, "rasterToMatch.tif"),
        paste0("_", studyAreaName)
      )
      rasterToMatch <- Cache(
        postProcessTo(
          from = rasterToMatchLarge,
          cropTo = studyArea, # needs to keep crs of original; can't use `to`
          maskTo = studyArea, # needs to keep crs of original; can't use `to`
          method = "bilinear",
          datatype = "INT2U",
          # writeTo = rtmFilename, # can't save w/ terra b/c same filename as RTML
          overwrite = TRUE
        ),
        userTags = c(cacheTags, "rasterToMatch"),
        omitArgs = c("destinationPath", "targetFile", "userTags", "stable", "writeTo", "overwrite")
      )
    }
    ## covert to 'mask'
    if (!anyNA(rasterToMatch[])) {
      whZeros <- as.vector(rasterToMatch[]) == 0
      if (sum(whZeros) > 0) {
        # means there are zeros instead of NAs for RTML --> change
        rasterToMatch[whZeros] <- NA
        message("There were no NAs on the RTM, but there were zeros; converting these zeros to NA.")
      }
    }

    RTMvals <- as.vector(rasterToMatch[])
    rasterToMatch[!is.na(RTMvals)] <- 1
  }

  return(list(rasterToMatch = rasterToMatch, rasterToMatchLarge = rasterToMatchLarge))
}
