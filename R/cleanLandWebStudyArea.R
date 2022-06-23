#' Clean up the LandWeb study area from David Andison
#'
#' @param poly A polygon or character string identifying the path to polygon
#' @param minFRI Numeric or integer, indicating the minimum fire return interval
#'               that will be part of the cleanup of polygon. Anything below
#'               this will be `NA`.
#' @export
#' @importFrom pemisc createPrjFile
#' @importFrom raster shapefile
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

#' Do an arbitrary set of operations on a polygon
#'
#' @param poly A polygon object, or a character string identifying the shapefile
#'             path to load, and clean.
#' @param fn   A function identifying the type of cleaning to do.
#' @param type If `fn` is not known, an character string can be specified to
#'             identify which `fn` to use.
#'             This MUST be a known type for this function.
#' @param ...  Passed to `fn`
#'
#' @export
#' @importFrom stats na.omit
polygonClean <- function(poly, fn = NULL, type = NULL, ...) {
  if (is.null(fn)) {
    if (is.null(type)) {
      stop("Either fn or type must be specified")
    } else {
      if (type == "LandWeb")
        fn <- .cleanLandWebStudyArea
      else
        stop("Unknown type")
    }
  }
  poly <- fn(poly, ...)
}
