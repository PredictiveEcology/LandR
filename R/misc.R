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
  newCohortData[ , lightProb := sufficientLight[cbind(shadetolerance, siteShade + 2)]]
}

#' Convert numeric values to rounded integers
#'
#' Simple wrapper around \code{as.integer} to round, rather than truncate, values.
#'
#' @param x A numeric vector
#'
#' @return An integer vector of length \code{x}, rounded to zero decimal places
#'   prior to \code{as.integer}
#'
#' @export
asInteger <- function(x)
  as.integer(round(x, 0))

#' Resample
#'
#' Imports the non-exported function \code{SpaDES.tools:::resample}.
#'
#' @keywords internal
#' @rdname resample
#' @seealso \code{\link[SpaDES.tools]{resample}}
.resample <- getFromNamespace("resample", "SpaDES.tools")

