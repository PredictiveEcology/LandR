#' \code{calculateClimateEffect}
#'
#' @param ... additional arguments that are passed to LandR.CS
#' @return NULL in place of a model object
#' @importFrom data.table data.table
#' @export
calculateClimateEffect <- function(...) {
  data.table('mortPred' = 100, 'growthPred' = 100)
}
