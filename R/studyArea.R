#' Create default study areas for use with LandR modules
#'
#' This simply re-exports \code{SpaDES.tools::randomStudyArea}
#'
#' @export
#' @importFrom utils getFromNamespace
randomStudyArea <- getFromNamespace("randomStudyArea", "SpaDES.tools")

#' Get and prep GADM maps
#'
#' A wrapper around \code{\link[raster]{getData}} that does additional cleanup and reprojects.
#'
#' @param country Character string giving the country code (default \code{"CAN"}).
#' @param level   Numeric giving the ADM level to use (default \code{1}).
#' @param proj     Character string giving the projection to use.
#' @param dPath   The destination path in which to save the data file.
#'
#' @export
#' @importFrom raster getData
#' @importFrom reproducible fixErrors
#' @importFrom sp spTransform
prepGADM <- function(country = "CAN", level = 1, proj, dPath) {
  getData("GADM", country = country, level = level, path = dPath) %>%
    spTransform(prj) %>%
    fixErrors(objectName = paste0("GDAM_", country, level))
}
