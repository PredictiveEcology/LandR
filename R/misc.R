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

#' Pull out the values from speciesEcoregion table for current time
#'
#' @param speciesEcoregion A \code{data.table} with \code{speciesEcoregion} values
#' @param currentTime The current time e.g., \code{time(sim)}
#'
#' @note
#' TODO
#'
#' @export
speciesEcoregionLatestYear <- function(speciesEcoregion, currentTime) {
  spEco <- speciesEcoregion[year <= currentTime]
  spEco[year == max(spEco$year)]
}
