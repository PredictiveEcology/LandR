if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", "age", "ecoregion", "ecoregionGroup",
                           "lightProb", "maxANPP", "maxB", "maxB_eco", "pixelIndex",
                           "shadetolerance", "siteShade", "speciesposition",
                           "speciesGroup", "speciesInt", "sumB",
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
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#'
addNewCohorts <- function(newCohortData, cohortData, pixelGroupMap, time, speciesEcoregion) {
  ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum biomass per ecoregion, join to new cohort data
  if (isTRUE(getOption("LandR.assertions"))) {
    ncd <- copy(newCohortData)
    CD <- copy(cohortData)
  }

  # This removes the duplicated pixels within pixelGroup, i.e., the reason we want pixelGroups
  newCohortData <- unique(newCohortData, by = c("pixelGroup", "speciesCode"))

  specieseco_current <- speciesEcoregion[year <= time]
  specieseco_current <- na.omit(specieseco_current)
  specieseco_current <- setkey(specieseco_current[year == max(specieseco_current$year),
                                                  .(speciesCode, maxANPP, maxB, ecoregionGroup)],
                               speciesCode, ecoregionGroup)

  # Note that after the following join, some cohorts will be lost due to lack of parameters in speciesEcoregion
  #  These need to be modified in pixelGroupMap
  missingNewCohortData <- newCohortData[!specieseco_current, on = c("speciesCode", "ecoregionGroup")]
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]
  newCohortData <- specieseco_current[newCohortData, on = c("speciesCode", "ecoregionGroup")]
  # newCohortData <- newCohortData[specieseco_current, on = c("speciesCode", "ecoregionGroup"),
  #                                nomatch = 0]
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

  # unique(cohortData, by = c("pixelGroup", "speciesCode", "age"))
  cohortData <- rbindlist(list(cohortData, newCohortData))
  return(cohortData)
}



#' Remove missing cohorts from cohortData based on pixelGroupMap
#'
#'
#' @param cohortData a \code{data.table} with columns:
#'   "pixelGroup", "ecoregionGroup", "speciesCode", "age", "B", "mortality", "aNPPAct", "sumB"
#' @param pixelGroupMap Raster layer with pixel values equal to a pixel group number that
#'   correspondsd exactly to ]\code{pixelGroup} column in \code{cohortData}
#' @param firePixelTable A data.table with 2 columns, \code{pixelIndex} and \code{pixelGroup}.
#'   This will be used in conjunction
#'   with cohortData and pixelGroupMap to ensure that everything matches correctly.
#'
#' @return
#' A \code{data.table} with a new, cohortData
#'
#' @export
#' @importFrom data.table rbindlist set setkey
#' @importFrom raster getValues
#'
rmMissingCohorts <- function(cohortData, pixelGroupMap, firePixelTable) {
  pgmVals <- na.omit(getValues(pixelGroupMap))
  pgsStillInCDGoneFromPGM <- setdiff(cohortData$pixelGroup, pgmVals)

  # REMOVE lines in cohortData that are no longer in the pixelGroupMap
  cohortData <- cohortData[!pixelGroup %in% pgsStillInCDGoneFromPGM]

  if (isTRUE(getOption("LandR.assertions"))) {
    # All pixels in pgsStillInCDGoneFromPGM should have been touched by a fire
    test1 <- isTRUE(all(pgsStillInCDGoneFromPGM %in% na.omit(firePixelTable$pixelGroup)))

    # There should still be some burned pixel groups that are still on the map, i.e., only some pixels from a PG got burned
    burnedByPGStillOnMap <- setdiff(na.omit(firePixelTable$pixelGroup), pgsStillInCDGoneFromPGM)
    test2 <- isTRUE(all(burnedByPGStillOnMap %in% pgmVals))

    test3 <- length(setdiff(cohortData$pixelGroup, pgmVals)) == 0

    if (!isTRUE(all(test1, test2, test3)))
      stop("cohortData and pixelGroupMap don't match")
  }

  cohortData
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


#' Create the correct string for pixelGroups
#'
#' @param maxPixelGroup A length 1 numeric/integer indicating the current maximum pixelGroup value
#' @param ecoregionGroup  A vector of ecoregionGroup strings
#' @param speciesGroup  A vector of speciesGroup strings
#'
#' @return  TODO: description needed
#'
#' @export
makePixelGroups <- function(maxPixelGroup, ecoregionGroup, speciesGroup) {
  as.integer(maxPixelGroup) +
    as.integer(factor(paste(ecoregionGroup, speciesGroup, sep = "_")))
}

#' Add the correct \code{pixelGroups} to a \code{pixelCohortData} object
#'
#' @inheritParams makePixelGroups
#' @param pixelCohortData  # pixel groups are groups of identical pixels based
#'   on \code{speciesGroup} x \code{Age} and \code{ecoregionGroup}.
#'
#' @note
#' This should not (yet) be used where age is an issue.
#'
#' @export
#' @importFrom data.table setkey
#' @importFrom SpaDES.core paddedFloatToChar
addPixelGroup <- function(pixelCohortData, maxPixelGroup) {

  pixelCohortData[, speciesInt := as.integer(speciesCode)]
  pixelCohortData[, speciesGroup := sum(2^(unique(speciesInt)-1)),  by = "pixelIndex"]
  pixelCohortData[, speciesGroup := paddedFloatToChar(speciesGroup, padL = max(nchar(as.character(speciesGroup))))]
  setkey(pixelCohortData, ecoregionGroup, speciesGroup)
  pixelCohortData[ , pixelGroup := makePixelGroups(maxPixelGroup, ecoregionGroup, speciesGroup)]
  pixelCohortData[, c("speciesInt", "speciesGroup") := NULL]
  pixelCohortData
}
