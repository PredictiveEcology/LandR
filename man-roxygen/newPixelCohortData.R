#' @param newPixelCohortData must be a complete \code{cohortData} object with newly created
#'        cohorts. They do not have to have \code{pixelGroup} values yet;
#'        they can be overlapping with \code{cohortData}, (i.e., they can be regenerated on empty
#'        pixels or on already occupied pixels).
#'        Must contain the columns: \code{pixelIndex}, \code{speciesCode}, \code{ecoregionGroup}.
#'        The remaining 4 (see \code{cohortData}) will be created with \code{0}s.
