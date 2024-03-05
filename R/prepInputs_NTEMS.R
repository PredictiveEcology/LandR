# utils::globalVariables(c(
#  "foo"
# ))

#' Obtain an LCC layer for a given year from NTEMS, with forest matching the FAO definition
#'
#' @inheritParams reproducible::postProcessTo
#' @param year stack of species layers rasters
#' @param rasterToMatch template raster to crop to
#' @param destinationPath destination folder for FAO and LCC data
#' @param disturbedCode value assigned to pixels that are forest per FAO definition but not in LCC year
#'
#' @return
#' A \code{SpatRaster} with corrected forest pixels
#'
#' @export
#' @importFrom reproducible prepInputs
#'
prepInputs_NTEMS_LCC_FAO <- function(year = 2010, disturbedCode = 0, rasterToMatch, ...) {
  if (year > 2019 | year < 1984) {
    stop("LCC for this year is unavailable")
  }

  #TODO: what is the cache option to avoid cachign these enormous 23 GB rasters

  # Data codes: 0 = no change; 20 = water; 31 = snow_ice; 32 = rock_rubble; 33 = exposed_barren_land;
  # 40 = bryoids; 50 = shrubs; 80 = wetland; 81 = wetland-treed; 100 = herbs; 210 = coniferous;
  # 220 = broadleaf; 230 = mixedwood

  lccURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", year, ".zip")
  lccTF <- paste0("CA_forest_VLCE2_", year, ".tif")
  lcc <- prepInputs(url = lccURL, targetFile = lccTF, method = "near", ...)

  #1 is forest, #2 is disturbed forest
  fao <- prepInputs(url = "https://opendata.nfis.org/downloads/forest_change/CA_FAO_forest_2019.zip",
                    method = "near", ...)
  #10 is not a class in use - make it disturbed forest
  #pixels may not be disturbed yet if year is prior to 2019 (FAO year)
  #adjust non-forest LCC that are disturbed forest to 10

  lcc[!lcc[] %in% c(210, 81, 220, 230) & fao[] == 2] <- disturbedCode

  return(lcc)
}

#'reclassify non-flammable pixels that become flammable - herbaceous or shrubby - vegetation
#'
#' @param rstLCC input lcc layer with bare soil class that may become vegetated
#' @param ... passed to prepInputs - crop, mask, and project will be to rstLCC
#' @param endYear NTEMS LCC year to use for correcting transition from bare to non-forest
#' @param lccToAdjust lcc values of the bare class
#' @param nonforestLCC allowable lcc values for bare to become
#' @return
#' A \code{SpatRaster} with non-flammable pixels corrected if they become flammable non-forest
#'
#' @export
#' @importFrom data.table data.table setnames
prepInputs_NTEMS_Nonforest <- function(rstLCC, endYear = 2019, lccToAdjust = 33,
                                       nonforestLCC = c(50, 100), ...) {

  lccURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", endYear, ".zip")
  lccTF <- paste0("CA_forest_VLCE2_", endYear, ".tif")
  endLCC <- prepInputs(url = lccURL, targetFile = lccTF, method = "near",
                       cropTo = rstLCC, projectTo = rstLCC, maskTo = rstLCC, ...)

  toFix <- data.table(currentLCC = rstLCC[], endLCC = endLCC[], id = 1:ncell(endLCC))
  setnames(toFix, new = c("currentLCC", "endLCC", "id"))
  toFix[currentLCC == lccToAdjust & endLCC %in% nonforestLCC, newLCC := endLCC]
  toFix <- toFix[!is.na(newLCC), .(id, newLCC)]

  #adjust pixels
  rstLCC[toFix$id] <- toFix$newLCC

  return(rstLCC)
}




