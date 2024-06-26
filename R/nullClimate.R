#' Null climate effect
#'
#' Default climate effects function in the case where no climate effect is simulated
#'
#' @template cohortData
#' @param ... additional arguments that are passed to LandR.CS
#'
#' @details the `cohortData` object is used to calculate the \%
#'    reduction/increase of mortality and growth biomasses per cohort.
#'
#' @return `data.table` with `pixelGroup`, `age` and `speciesCode`, as well as
#'   `mortPred` and `growthPred` columns with \% reduction/increase of mortality and
#'   growth biomasses resulting from a climate effect.
#'   These percentages are later multiplied by the by baseline biomasses of mortality and growth
#'   (e.g. 0\% meaning total mortality or growth reduction, and 100\% meaning no reduction).
#'   This default, no climate effect, function outputs 100\% for all cohorts for both
#'   `mortPred` and `growthPred`.
#'
#' @export
calculateClimateEffect <- function(cohortData, ...) {
  predObj <- unique(cohortData[, .(pixelGroup, age, speciesCode)])
  predObj[, `:=`(mortPred = 100, growthPred = 100)]
  return(predObj)
}
