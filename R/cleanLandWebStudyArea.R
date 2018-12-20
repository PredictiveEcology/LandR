#' Clean up the LandWeb study area from David Andison
#'
#' @param poly A polygon or character string identifying the path to polygon
#' @param minFRI Numeric or integer, indicating the minimum fire return interval
#'               that will be part of the cleanup of polygon. Anything below
#'               this will be \code{NA}.
#' @export
#' @importFrom raster shapefile
#'
.cleanLandWebStudyArea <- function(poly, minFRI = 40) {
  if (is.character(poly)) {
    createPrjFile(poly)
    poly <- raster::shapefile(poly)
  }

  stopifnot(any(c("LTHFC", "LTHRC") %in% names(poly)))
  if (!isTRUE("LTHRC" %in% names(poly))) {
    # Apparently, sometimes it is LTHFC, sometimes LTHRC; get rid of LTHFC
    poly$LTHRC <- poly$LTHFC #nolint
    poly$LTHFC <- NULL #nolint

    # The fires of Fire Return Interval 30 years are not correctly simulated
    # by LandMine, so they are removed.
    poly$fireReturnInterval <- poly$LTHRC
    poly$LTHRC <- NULL #nolint
  }
  poly$fireReturnInterval[poly$fireReturnInterval <= minFRI] <- NA
  poly@data <- poly@data[, !(names(poly) %in% "ECODISTRIC")]

  poly
}
