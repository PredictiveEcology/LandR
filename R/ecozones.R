#' Prepare ecozones raster
#'
#' Download, rasterize, crop, mask, and reproject Canadian national ecozones shapefile.
#'
#' @param url character. ecozones shapefile url.
#'            Default: <http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip>.
#'
#' @template destinationPath
#'
#' @template studyArea
#'
#' @template rasterToMatch
#'
#' @return `SpatRaster`
#'
#' @export
prepEcozonesRst <- function(url = NULL, destinationPath, studyArea = NULL, rasterToMatch = NULL) {
  if (is.null(url)) {
    url <- "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
  }

  targetCRS <- NULL

  if (!is.null(rasterToMatch)) {
    targetCRS <- sp::proj4string(rasterToMatch)
  }

  if (!is.null(studyArea)) {
    studyArea <- sf::as_Spatial(studyArea)
    targetCRS <- sp::proj4string(studyArea)
  }

  ecozone_shp <- prepInputs(
    url = url,
    targetFile = "ecozones.shp",
    alsoExtract = "similar",
    fun = "sf::st_read",
    destinationPath = destinationPath,
    studyArea = studyArea,
    targetCRS = targetCRS
  )

  ecozone_shp[["ZONE_NAME"]] <- as.factor(ecozone_shp[["ZONE_NAME"]])
  ecozone <- terra::rasterize(
    ecozone_shp,
    terra::rast(rasterToMatch),
    field = "ZONE_NAME",
    fun = "sum"
  ) |>
    terra::rast() |>
    terra::as.int()

  if (is(rasterToMatch, "Raster")) {
    rasterToMatch <- terra::rast(rasterToMatch)
  }

  ecozone <- terra::project(ecozone, rasterToMatch, method = "near")
  ecozone <- terra::mask(ecozone, rasterToMatch)

  return(ecozone)
}
