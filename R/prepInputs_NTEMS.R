utils::globalVariables(c(
  "currentLCC", "endLCC", "pixelID"
))

#' Obtain an LCC layer for a given year from NTEMS, with forest matching the FAO definition
#'
#' @param year stack of species layers rasters
#' @param disturbedCode value assigned to pixels that are forest per FAO definition but not in LCC year
#' @param ... passed to `prepInputs`
#'
#' @return a `SpatRaster` with corrected forest pixels
#'
#' @export
prepInputs_NTEMS_LCC_FAO <- function(year = 2010, disturbedCode = 1, ...) {
  if (year > 2019 | year < 1984) {
    stop("LCC for this year is unavailable")
  }

  dots <- list(...)

  if (is.null(dots$rasterToMatch) && is.null(dots$cropTo)) {
    stop("the NTEMS raster file is too large to process without cropping via `rasterToMatch` or `cropTo`")
  }

  if (isTRUE(getOption("reproducible.gdalwarp"))) {
    message("temporarily setting reproducible.usegdalwarp to FALSE to avoid error")
    options("reproducible.gdalwarp" = FALSE)
    on.exit(options("reproducible.gdalwarp" = TRUE), add = TRUE)
  }

  ## Data codes:
  ## 0 = no change; 20 = water; 31 = snow_ice; 32 = rock_rubble; 33 = exposed_barren_land;
  ## 40 = bryoids; 50 = shrubs; 80 = wetland; 81 = wetland-treed; 100 = herbs; 210 = coniferous;
  ## 220 = broadleaf; 230 = mixedwood
  lccURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", year, ".zip")
  lccTF <- paste0("CA_forest_VLCE2_", year, ".tif")
  lcc <- prepInputs(url = lccURL, targetFile = lccTF, method = "near", ...)

  ## 1 is forest, 2 is disturbed forest
  ## do not pass dots, or the filename is passed and is overwritten
  fao <- prepInputs(url = "https://opendata.nfis.org/downloads/forest_change/CA_FAO_forest_2019.zip",
                    method = "near", destinationPath = dots$destinationPath, cropTo = lcc,
                    maskTo = lcc, projectTo = lcc)

  ## 10 is not a class in use - make it disturbed forest
  ## pixels may not be disturbed yet if year is prior to 2019 (FAO year)
  ## adjust non-forest LCC that are disturbed forest to 10
  lccDat <- data.table(pixelID = 1:ncell(lcc), lcc = values(lcc, mat = FALSE))
  lccDat <- lccDat[!is.na(lcc) & lcc %in% c(210, 81, 220, 230)]
  lccDat[, fao := values(fao, mat = FALSE)[pixelID]]
  lccDat <- lccDat[!is.na(fao) & fao == 2]

  lcc[lccDat$pixelID] <- disturbedCode
  rm(lccDat)
  gc()

  return(lcc)
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
    options("reproducible.gdalwarp" = FALSE)
    on.exit(options("reproducible.gdalwarp" = TRUE), add = FALSE)
  }

  lccURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", endYear, ".zip")
  lccTF <- paste0("CA_forest_VLCE2_", endYear, ".tif")
  finalLCC <- prepInputs(url = lccURL, targetFile = lccTF, method = "near",
                         cropTo = rstLCC, projectTo = rstLCC, maskTo = rstLCC, ...)

  toFix <- data.table(currentLCC = values(rstLCC, mat = FALSE), pixelID = seq_len(ncell(finalLCC)))
  toFix <- toFix[currentLCC %in% lccToAdjust]

  toFix[, endLCC := values(finalLCC, mat = FALSE)[pixelID]]
  toFix <- toFix[endLCC %in% nonforestLCC, newLCC := endLCC]
  toFix <- toFix[!is.na(newLCC)]

  #adjust pixels
  rstLCC[toFix$pixelID] <- toFix$newLCC

  return(rstLCC)
}
