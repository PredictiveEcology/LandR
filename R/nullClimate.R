#' \code{calculateClimateEffect}
#'
#' @param ... additional arguments that are passed to LandR.CS
#' @return NULL in place of a model object
#' @importFrom data.table data.table
#' @export
calculateClimateEffect <- function(cohortData, ...) {
  predObj <- unique(cohortData[, .(pixelGroup, age, speciesCode)])
  predObj[, `:=`(mortPred = 100, growthPred = 100)]
  return(predObj)
}
