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
#' Simple wrapper to round, rather than truncate, values.
#'
#' @note Values ending in .5 will be rounded up, which is reasonable for positive values of \code{x},
#' but may not be desired for negative values.
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
