#' Get digital elevation map of Canada
#'
#' Defaults to using 3 arcmin DEM of Canada in lonlat.
#'
#' @template studyArea
#' @template rasterToMatch
#' @template destinationPath
#'
#' @return `RasterLayer`
#'
#' @export
#' @importFrom raster raster
#' @importFrom reproducible postProcess prepInputs projectInputs
prepInputsCanDEM <- function(studyArea, rasterToMatch, destinationPath) {
  dem_url <- "https://drive.google.com/file/d/1Cd4J5I2LHGUJosDQgfKrb_T93ZVSWCkt/"

  lonlat <- CRS("+proj=longlat +datum=WGS84")
  studyArea_lonlat <- projectInputs(studyArea, lonlat)

  dem <- prepInputs(
    targetFile = "Canada_3ArcMinuteDEM.asc",
    url = dem_url,
    destinationPath = destinationPath,
    fun = "raster::raster", method = "bilinear"
  )
  crs(dem) <- lonlat
  dem <- postProcess(dem, studyArea = studyArea_lonlat, rasterToMatch = rasterToMatch)
}
