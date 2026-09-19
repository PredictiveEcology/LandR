#' Create default study areas for use with LandR modules
#'
#' This simply re-exports `SpaDES.tools::randomStudyArea`
#'
#' @inheritParams SpaDES.tools::randomStudyArea
#'
#' @export
randomStudyArea <- utils::getFromNamespace("randomStudyArea", "SpaDES.tools")

#' Ecological boundaries within a study area
#'
#' Extend a study area boundary to include the full ecological boundary polygons intersecting it.
#'
#' @template studyArea
#'
#' @template destinationPath
#'
#' @param type character. The polygon type to use:
#'             one of "ECODISTRICT", "ECOREGION", "ECOPROVINCE", or "ECOZONE".
#'
#' @returns spatial polygons object of the same class as `studyArea`
#'
#' @examplesIf !isTRUE(as.logical(Sys.getenv("CI")))
#' ## NOTE: runs locally, skipped on CI -- downloads from sis.agr.gc.ca, which throttles
#' ## connections from CI runners (same predicate as testthat::skip_on_ci())
#' ## using SpatVector objects
#' sa <- randomStudyArea(size = 1e9)
#' sa_eco <- studyAreaEco(studyArea = sa)
#'
#' ## using sf objects
#' sa_sf <- sf::st_as_sf(sa)
#' sa_eco_sf <- studyAreaEco(studyArea = sa_sf)
#'
#' if (interactive()) {
#'   ggplot() +
#'     geom_sf(data = sa_eco_sf, fill = "gray") +
#'     geom_sf(data = sa_sf, fill = "violet", alpha = 0.3)
#' }
#'
#' @export
studyAreaEco <- function(studyArea = NULL, destinationPath = tempdir(),
                         type = c("ECOZONE", "ECOPROVINCE", "ECOREGION", "ECODISTRICT")) {
  stopifnot(is.null(studyArea) || inherits(studyArea, "sf") || inherits(studyArea, "SpatVector"))

  is_sf <- inherits(studyArea, "sf")

  if (!is.null(studyArea) && !is_sf) {
    studyArea <- sf::st_as_sf(studyArea)
  }

  url <- switch(
    tolower(type[1]),
    ecodistrict = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    ecoregion = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
    ecoprovince = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
    ecozone = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
  )
  eco <- prepInputs(
    url = url,
    destinationPath = destinationPath,
    projectTo = studyArea,
    fun = "sf::st_read"
  )

  if (!is.null(studyArea)) {
    eco <- eco[which(sapply(sf::st_intersects(eco, studyArea), length) > 0), ]

    if (!is_sf) {
      eco <- terra::vect(eco)
    }
  }

  return(eco)
}
