if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", "lightProb", "shadetolerance", "siteShade", "year"))
}

#' Assign light probability
#'
#' @param sufficientLight TODO: description needed
#' @param newCohortData  TODO: description needed
#'
#' @return  TODO: description needed
#'
#' @export
assignLightProb <- function(sufficientLight, newCohortData) {
  ## for each line, get the survival probability from sufficientLight table note
  ## that sufficentLight is a table of survival probs for each tolerance level
  ## (row) by and shade level (column) siteShade + 2 is necessary to skip the
  ## first column
  newCohortData[, lightProb := sufficientLight[cbind(shadetolerance, siteShade + 2)]]
}

#' Convert numeric values to rounded integers
#'
#' Essentially a wrapper around \code{round}, rather than \code{truncate}, which is what \code{as.integer}
#' does. Internally, this is simply \code{as.integer(floor(x + 0.5))}.
#'
#' @note Values ending in \code{.5} will be rounded up, whether positive or negative.
#' This is different than \code{round}.
#'
#' @param x A numeric vector
#'
#' @return An integer vector of length \code{x}, rounded to zero decimal places
#'   prior to \code{as.integer}
#'
#' @export
#' @examples
#' x <- seq(-2, 2, 0.25)
#' data.frame(dbl = x, int = asInteger(x))
asInteger <- function(x)
  as.integer(floor(x + 0.5))

#' Resample
#'
#' Imports the non-exported function \code{SpaDES.tools:::resample}.
#'
#' @importFrom utils getFromNamespace
#' @keywords internal
#' @rdname resample
#' @seealso \code{\link[SpaDES.tools]{resample}}
.resample <- getFromNamespace("resample", "SpaDES.tools")

#' Test whether disturbance should be scheduled
#'
#' @param disturbanceLayer a \code{RasterLayer} object
#' @param currentYear time of simulation
#' @param disturbanceType the type of disturbance (e.g. Burn, Harvest)
#' @return logical vector signifying whether to schedule event
#' @export
#'
#' @examples
#' doEvent <- scheduleDisturbance(sim$rstCurrentBurn, time(sim), disturbanceType = "Burn")
scheduleDisturbance <- function(disturbanceLayer, currentYear, disturbanceType){
  if (!disturbanceType %in% c("Burn", "Harvest")) {
    stop("Please ensure disturbance type is either 'Burn' or 'Harvest'")
  }
  if (is.null(disturbanceLayer) || names(disturbanceLayer) != paste0(disturbanceType, currentYear)) {
    TRUE
  } else {
    FALSE
  }
}
