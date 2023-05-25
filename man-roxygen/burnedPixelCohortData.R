#' @param burnedPixelCohortData An expanded `cohortData` `data.table` with pixel-level
#'   cohort information on burnt pixels and the following (optional) columns:
#'   `severity` - fire severity in that pixel calculated based on fire behaviour properties;
#'   `firetolerance` - species-level fire tolerance;
#'   `severityToleranceDif` - the difference between `severity` and `firetolerance`.
