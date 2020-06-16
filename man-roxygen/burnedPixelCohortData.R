#' @param burnedPixelCohortData An expanded \code{cohortData} \code{data.table} with pixel-level
#'   cohort information on burnt pixels and columns:
#'   \code{severity} - fire severity in that pixel calculated based on fire behaviour properties;
#'   \code{firetolerance} - species-level fire tolerance;
#'   \code{severityToleranceDif} - the difference between \code{severity} and \code{firetolerance}.
