utils::globalVariables(c(
  "currentLCC", "destinationPath", "endLCC", "pixelID", "writeTo"
))

#' Obtain an LCC layer for a given year from NTEMS, with forest matching the FAO definition
#'
#' @param year stack of species layers rasters
#' @param disturbedCode value assigned to pixels that are forest per FAO definition but not in LCC year
#' @param resampleMethod method used when resampling LCC layers to match `rasterToMatch`
#' @param ... passed to `prepInputs`
#'
#' @return a `SpatRaster` with corrected forest pixels
#'
#' @importFrom terra lapp ncell inMemory
#' @export
prepInputs_NTEMS_LCC_FAO <- function(year = 2010, disturbedCode = 240, resampleMethod = "near", ...) {

  if (year > 2019 || year < 1984) {
    stop("LCC for this year is unavailable")
  }
  newFilename <- NULL
  writeToFN <- NULL
  dots <- list(...)
  if (!is.null(dots$writeTo)) {
    #must pass a different file name to prepInputs as the object that is Cached
    #will inevitably be modified later in this function
    writeToFN <- dots$writeTo
    #assign a temporary filename for the raw LCC
    newFilename <- paste0("raw_", basename(dots$writeTo))
    dots$writeTo <- NULL
  }

  if (is.null(dots$rasterToMatch) && is.null(dots$cropTo) && is.null(dots$to)) {
    warning("the NTEMS raster file is too large to process without cropping via `rasterToMatch` or `cropTo`")
  }

  if (isTRUE(getOption("reproducible.gdalwarp"))) {
    message("temporarily setting reproducible.usegdalwarp to FALSE to avoid error")
    opts <- options(reproducible.gdalwarp = FALSE)
    on.exit(options(opts), add = TRUE)
  }
  ## Data codes:
  ## 0 = no change; 20 = water; 31 = snow_ice; 32 = rock_rubble; 33 = exposed_barren_land;
  ## 40 = bryoids; 50 = shrubs; 80 = wetland; 81 = wetland-treed; 100 = herbs; 210 = coniferous;
  ## 220 = broadleaf; 230 = mixedwood
  lccURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", year, ".zip")
  lccTF <- paste0("CA_forest_VLCE2_", year, ".tif")

  #fix dots
  dots$url <- lccURL
  dots$targetFile <- lccTF
  dots$method <- resampleMethod
  # dots$writeTo <- newFilename
  lcc <- do.call(prepInputs, dots)
  lcc <- buildVRT(lcc, writeTo = newFilename,
                         destinationPath = dots$destinationPath)

  if (!inMemory(lcc)) {
    faoFilename <- paste0("FAO_", dots$writeTo)
  } else {
    faoFilename <- NULL
  }

  # dots$writeTo <- writeToFN

  ## 2024-12: see #110; don't delete CA_forest_VLCE2 raster even though it's 24GB
  ## deleting it results in redownload every time and breaks parallel sims (race condition)
  # toUnlink <- ifelse(is.null(dots$destinationPath), lccTF,
  #                    file.path(dots$destinationPath, lccTF))
  # unlink(toUnlink)

  #restore dots$writeTo - it will be NULL if it wasn't passed

  ## 1 is forest, 2 is land that can meet the FAO definition of forest
  ## do not pass dots, or the filename is passed and is overwritten
  url <- "https://opendata.nfis.org/downloads/forest_change/CA_FAO_forest_2019.zip"
  #let terra options dictate whether fao is on disk or not

  fao <- prepInputs(
    url = url,
    method = resampleMethod, destinationPath = dots$destinationPath, cropTo = lcc,
    maskTo = lcc, projectTo = lcc
  )
  fao <- buildVRT(fao, writeTo = faoFilename, destinationPath = dots$destinationPath)

  ## pixels may not be disturbed yet if year is prior to 2019 (FAO year)
  ## adjust non-forest LCC that are disturbed forest to disturbedCode

  DisturbedAdjust <- function(LCC, FAO, newVal = disturbedCode) {
    LCC[FAO == 2 & !LCC %in% c(210, 81, 220, 230)] <- newVal
    return(LCC)
  }

  input <- c(lcc, fao)
  out <- terra::lapp(input, fun = DisturbedAdjust, usenames = FALSE)
  # lcc <- terra::init(lcc, as.vector(out))

  #assign it to itself or it stays in memory
  out <- buildVRT(out, writeTo = writeToFN,
                  destinationPath = dots$destinationPath) #overwrite lcc

  gc()
  return(out)
}

#' Reclassify non-flammable pixels that become flammable - herbaceous or shrubby - vegetation
#'
#' @param rstLCC input lcc layer with bare soil class that may become vegetated
#' @param endYear NTEMS LCC year to use for correcting transition from bare to non-forest
#' @param lccToAdjust lcc values of the bare class
#' @param nonforestLCC allowable lcc values for bare to become
#' @param ... non-spatial arguments passed to `prepInputs` e.g. `destinationPath`
#' @return a `SpatRaster` with non-flammable pixels corrected if they become flammable non-forest
#'
#' @export
prepInputs_NTEMS_Nonforest <- function(rstLCC, endYear = 2019, lccToAdjust = 33,
                                       nonforestLCC = c(50, 100), ...) {
  if (is.null(rstLCC)) {
    ## allow a more graceful fail than imploding upon reading the NTEMS dataset
    stop("rstLCC should not be NULL")
  }

  if (isTRUE(getOption("reproducible.gdalwarp"))) {
    message("temporarily setting reproducible.usegdalwarp to FALSE to avoid error")
    opts <- options(reproducible.gdalwarp = FALSE)
    on.exit(options(opts), add = TRUE)
  }

  lccURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", endYear, ".zip")
  lccTF <- paste0("CA_forest_VLCE2_", endYear, ".tif")
  finalLCC <- prepInputs(
    url = lccURL, targetFile = lccTF, method = "near",
    cropTo = rstLCC, projectTo = rstLCC, maskTo = rstLCC, ...
  )

  toFix <- data.table(currentLCC = values(rstLCC, mat = FALSE), pixelID = seq_len(ncell(finalLCC)))
  toFix <- toFix[currentLCC %in% lccToAdjust]

  toFix[, endLCC := values(finalLCC, mat = FALSE)[pixelID]]
  toFix <- toFix[endLCC %in% nonforestLCC, newLCC := endLCC]
  toFix <- toFix[!is.na(newLCC)]

  ## adjust pixels
  rstLCC[toFix$pixelID] <- toFix$newLCC

  return(rstLCC)
}

#' Obtain an Dominant species layer for a given year from NTEMS
#'
#' @template destinationPath
#' @param year stack of species layers rasters
#' @param ... passed to `prepInputs`
#' @template sppEquiv
#' @template sppEquivCol
#' @return a `SpatRaster` with dominant species
#' @importFrom terra unique
#'
#' @export
prepInputs_NTEMS_DominantSpecies <- function(year = 2011, destinationPath, sppEquiv = LandR::sppEquivalencies_CA,
                                             sppEquivCol = "LandR", ...) {
  if (year > 2022 || year < 1984) {
    stop("Dominant species for this year is unavailable")
  }

  dots <- list(...)

  if (is.null(dots$rasterToMatch) && is.null(dots$cropTo)) {
    warning("the NTEMS raster file is large and will take significant time to prepare without `rasterToMatch` or `cropTo` defined")
  }

  if (!is.null(dots$cropTo)) {
    cropTo <- dots$cropTo
  } else {
    cropTo <- dots$rasterToMatch
  }

  if (!is.null(dots$maskTo)) {
    maskTo <- dots$maskTo
  } else {
    maskTo <- dots$rasterToMatch
  }


  if (isTRUE(getOption("reproducible.gdalwarp"))) {
    message("temporarily setting reproducible.usegdalwarp to FALSE to avoid error")
    opts <- options(reproducible.gdalwarp = FALSE)
    on.exit(options(opts), add = TRUE)
  }

  domSppURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_Tree_Species_Classification_", year, ".zip")
  domSppTF <- paste0("Canada_Tree_Species_Classification_HMM_", year, ".tif")
  domSpp <- prepInputs(
    url = domSppURL, targetFile = domSppTF, # cropping and masking nation wide raster to study area but NOT reprojecting
    destinationPath = destinationPath,
    cropTo = cropTo, maskTo = maskTo
  )

  sppEquiv <- sppEquiv[, .SD, .SDcol = c("NTEMS_Species_Code", sppEquivCol)] # matching NTEMS spp code to sppEquivCol
  uniqueVals <- as.data.table(terra::unique(domSpp))
  setnames(uniqueVals, new = "NTEMS_Species_Code")
  uniqueVals <- sppEquiv[uniqueVals, on = c("NTEMS_Species_Code")] # pulling all species from NTEMS layer
  uniqueVals <- na.omit(uniqueVals) # removing non-treed areas

  if (!is.null(dots$projectTo)) {
    projectTo <- dots$projectTo
  } else if (!is.null(dots$rasterToMatch)) {
    projectTo <- dots$rasterToMatch
  } else {
    projectTo <- NULL
  }

  domSpp <- lapply(uniqueVals[["NTEMS_Species_Code"]], FUN = function(spp, ras = domSpp,
                                                                      template = projectTo) { # converting each species to binary layers and reprojecting to save user computation time
    newMap <- domSpp
    newMap[!domSpp[] == spp] <- 0
    newMap[domSpp[] == spp] <- 1

    if (!is.null(template)) {
      newMap <- postProcess(newMap, to = template, method = "average")
    }
    return(newMap)
  })

  domSpp <- rast(domSpp)
  names(domSpp) <- uniqueVals[[sppEquivCol]]

  return(domSpp)
}
