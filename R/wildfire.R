#' Download and prepare raster fire data from NFI
#'
#' @param dPath destination path for archive
#' @template rasterToMatch
#' @param url location from which to download the wildfire raster layer(s)
#'
#' @return a raster with values representing fire year 1985-2015
#'
#' @export
getWildfire_NFI <- function(dPath, rasterToMatch, url = NULL) {
  if (is.null(url)) {
    url <- paste0(
      "https://opendata.nfis.org/downloads/forest_change/",
      "CA_forest_wildfire_year_DNBR_Magnitude_1985_2015.zip"
    )
  }

  ## this is an enormous raster - we want the second raster in the stack, wildfire year
  fireDL <- preProcess(
    url = url,
    destinationPath = dPath,
    targetFile = raster::extension(basename(url), "tif"),
    alsoExtract = "similar"
  )
  wildfireYear <- terra::rast(fireDL$targetFilePath)
  wildfireYear2 <- wildfireYear$CA_forest_wildfire_year_DNBR_Magnitude_1985_2015_2 ## only need year

  ## TODO: Caching of postProcessTerra() directly is still experimental/inconsistent; avoid for now
  wildfire_SA <- postProcessTerra(wildfireYear2, to = rasterToMatch, method = "near")
  wildfire_SA[wildfire_SA == 0L] <- NA_integer_
  wildfire_SA <- wildfire_SA + 1900
  wildfire_SA <- terra::as.int(wildfire_SA)

  return(wildfire_SA)
}
