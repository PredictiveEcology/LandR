if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", "age", "ecoregion", "ecoregionGroup",
                           "lightProb", "maxANPP", "maxB", "maxB_eco", "pixelIndex",
                           "shadetolerance", "siteShade", "speciesposition", "sumB",
                           "temppixelGroup", "year"))
}

#' Add new cohorts
#'
#' Does the following:
#' \enumerate{
#'   \item add new cohort data into \code{cohortdata};
#'   \item assign initial biomass and age for new cohort;
#'   \item assign the new pixelgroup to the pixels that have new cohort;
#'   \item update the pixelgroup map.
#' }
#'
#' @param newCohortData must be a complete cohortData object with entirely new
#'   \code{pixelGroup} values (i.e., already non-overlapping with
#'   \code{cohortData}. Columns it must have: \code{pixelGroup}, \code{speciesCode},
#'   \code{ecoregionGroup}
#' @param cohortData a \code{data.table} with columns:
#'   "pixelGroup", "ecoregionGroup", "speciesCode", "age", "B", "mortality", "aNPPAct", "sumB"
#' @param pixelGroupMap Raster layer with pixel values equal to a pixel group number that
#'   correspondsd exactly to ]\code{pixelGroup} column in \code{cohortData}
#' @param time Current time e.g., time(sim). This is used to extract the correct
#'   parameters in \code{speciesEcoregion} table if there are different values over time
#' @param speciesEcoregion A speciesEcoregion table.
#'
#' @return
#' A \code{data.table} with a new, \code{rbindlist}ed cohortData
#'
#' @export
#' @importFrom data.table rbindlist set setkey
#' @importFrom raster getValues
#'
addNewCohorts <- function(newCohortData, cohortData, pixelGroupMap, time, speciesEcoregion) {
  ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum biomass per ecoregion, join to new cohort data
  specieseco_current <- speciesEcoregion[year <= time]
  specieseco_current <- na.omit(specieseco_current)
  specieseco_current <- setkey(specieseco_current[year == max(specieseco_current$year),
                                                  .(speciesCode, maxANPP, maxB, ecoregionGroup)],
                               speciesCode, ecoregionGroup)
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]
  newCohortData <- newCohortData[specieseco_current, on = c("speciesCode", "ecoregionGroup"), nomatch = 0]
  #newCohortData <- setkey(newCohortData, speciesCode, ecoregionGroup)[specieseco_current, nomatch = 0]
  set(newCohortData, NULL, "age", 1)  ## set age to 1
  set(newCohortData, NULL, "sumB", 0)
  ## set biomass - if B=0, it's getting maxANPP ???
  set(newCohortData, NULL, "B",
      as.integer(pmax(1, newCohortData$maxANPP *
                        exp(-1.6 * newCohortData$sumB / newCohortData$maxB_eco))))
  set(newCohortData, NULL, "B", as.integer(pmin(newCohortData$maxANPP, newCohortData$B)))

  newCohortData <- newCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, age, B,
                                     mortality = 0, aNPPAct = 0, sumB = 0)]
  cohortData <- rbindlist(list(cohortData, newCohortData))
  pixelGroupsLostToDisturbance <- setdiff(cohortData$pixelGroup, getValues(pixelGroupMap))
  cohortData <- cohortData[!pixelGroup %in% pixelGroupsLostToDisturbance]
  return(cohortData)
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
  ## for each line, get the survival probability from sufficientLight table
  ## note that sufficentLight is a table of survival probs for each tolerance level (row) by and shade level (column)
  ## siteShade + 2 is necessary to skip the first column
  newCohortData[ , lightProb := sufficientLight[cbind(shadetolerance, siteShade + 2)]]
}
