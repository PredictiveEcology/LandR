if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", "age", "aNPPAct", "ecoregion", "ecoregionGroup",
                           "maxANPP", "maxB", "maxB_eco", "mortality", "pixelIndex",
                           "speciesposition", "speciesGroup", "speciesInt", "sumB",
                           "temppixelGroup", "uniqueCombo", "uniqueComboByRow", "uniqueComboByPixelIndex",
                           "year"))
}

#' Add cohorts to \code{cohortData} and \code{pixelGroupMap}
#'
#' This is a wrapper for  \code{generatePixelGroups}, \code{initiateNewCohort} and
#' updates to \code{pixelGroupMap} via assignment to new \code{pixelIndex}
#' values in \code{newCohortData}. By running these all together,
#' there is less chance that they will diverge. There are some checks
#' internally for consistency.
#'
#' Does the following:
#' \enumerate{
#'   \item add new cohort data into \code{cohortdata};
#'   \item assign initial biomass and age for new cohort;
#'   \item assign the new pixelgroup to the pixels that have new cohort;
#'   \item update the pixelgroup map.
#' }
#'
#' @param newCohortData must be a complete cohortData object with newly
#' created cohorts. They do not have to have \code{pixelGroup} values yet; they
#' can be overlapping with \code{cohortData}, (i.e., they can be regenerated on
#' empty pixels or on alreayd occupied pixels). Columns it must have:
#' \code{pixelGroup}, \code{speciesCode}, \code{age}, \code{ecoregionGroup}. The
#' remaining 4 (see \code{cohortData}) will be created with \code{0}s
#' @param cohortData a \code{data.table} with columns: \code{pixelGroup},
#'   \code{ecoregionGroup}, \code{speciesCode}, \code{age}, \code{B},
#'   \code{mortality}, \code{aNPPAct}, and \code{sumB}.
#' @param pixelGroupMap Raster layer with pixel values equal to a pixel group
#'   number that correspondsd exactly to \code{pixelGroup} column in
#'   \code{cohortData}
#' @param time Current time e.g., time(sim). This is used to extract the correct
#'   parameters in \code{speciesEcoregion} table if there are different values
#'   over time
#' @param speciesEcoregion A \code{speciesEcoregion} table.
#' @param firePixelTable A data.table with 2 columns, \code{pixelIndex} and
#'   \code{pixelGroup}
#' @param successionTimestep The time between successive seed dispersal events.
#'   In LANDIS-II, this is called "Succession Timestep". This is used here
#'
#' @return
#' A list of length 2, \code{cohortData} and \code{pixelGroupMap}, with
#' \code{newCohortData} inserted.
#'
#' @export
#' @rdname updateCohortData
#' @importFrom crayon green magenta
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
updateCohortData <- function(newCohortData, cohortData, pixelGroupMap, time, speciesEcoregion,
                             firePixelTable = NULL, successionTimestep) {

  maxPixelGroup <- as.integer(maxValue(pixelGroupMap))

  if (!is.null(firePixelTable)) {
    pixelGroupMap[firePixelTable$pixelIndex] <- 0L
  }
  relevantPixels <- pixelGroupMap[][newCohortData$pixelIndex]
  zeroOnPixelGroupMap <- relevantPixels == 0

  if (!"age" %in% colnames(newCohortData))
    newCohortData[, age := 1L]

  allNewPixelGroups <- all(zeroOnPixelGroupMap)
  if (all(zeroOnPixelGroupMap)) {
    # Deal with pixels on the map that have no pixelGroup -- these are burned
    # pixels --> the entirely newly regenerated pixels DOes not require a
    # re-pixelGroupMaping  -- can just add to existing pixelGroup values
    message(crayon::green("  Regenerating only open pixels (e.g., likely resprouting & serotiny only)"))
    columnsForPG <- c("ecoregionGroup", "speciesCode", "age") # no Biomass because they all have zero
    cd <- newCohortData[,c("pixelIndex", columnsForPG), with = FALSE]
    newCohortData[, pixelGroup :=
                    generatePixelGroups(cd, maxPixelGroup = maxPixelGroup,
                                        columns = columnsForPG)]#,
    #successionTimestep = successionTimestep)

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(newCohortData[, c("pixelIndex", "pixelGroup")],
                             by = c("pixelIndex"))
  } else {
    # This is for situations where there are some empty pixels being filled,
    #   and some occupied pixels getting infilling. This requires a wholesale
    #   re-pixelGroup
    message(crayon::green("  Regenerating open and pixels with biomass (likely after seed dispersal)"))

    allNewPixelGroups <- FALSE

    pixelIndex <- which(pixelGroupMap[] %in% cohortData$pixelGroup)
    cohortDataPixelIndex <- data.table(pixelIndex = pixelIndex,
                                       pixelGroup = pixelGroupMap[][pixelIndex])
    cdLong <- cohortDataPixelIndex[cohortData, on = "pixelGroup", allow.cartesian = TRUE]
    cohorts <- rbindlist(list(cdLong, newCohortData), use.names = TRUE, fill = TRUE)

    columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B")
    cd <- cohorts[,c("pixelIndex", columnsForPG), with = FALSE]
    cohorts[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L,
                                                columns = columnForPG)]

    if (isTRUE(getOption("LandR.assertions"))) {

      uniquePixelsInCohorts <- pixelGroupMap[][unique(cohorts$pixelIndex)]
      pixelsOnMap <- sum(!is.na(pixelGroupMap[]), na.rm = TRUE)
      lenUniquePixelsInCohorts <- length(unique(cohorts$pixelIndex))
      lenBurnedPixels <- sum(pixelGroupMap[] == 0, na.rm = TRUE) # 30927
      pixelsRegeneratedOnZeros <- sum(uniquePixelsInCohorts == 0)
      allPixelsNotInCohortData <- pixelGroupMap[][-unique(cohorts$pixelIndex)]
      numPixelsNoRegen <- sum(allPixelsNotInCohortData == 0, na.rm = TRUE)
      tableB <- table(allPixelsNotInCohortData) # 25166

      test2 <- identical(as.integer(pixelsOnMap - tableB), lenUniquePixelsInCohorts)
      test3 <- identical(as.integer(pixelsOnMap - (lenBurnedPixels - pixelsRegeneratedOnZeros)), lenUniquePixelsInCohorts)

      uniqueAllPixelsNotInCohortData <- unique(allPixelsNotInCohortData)
      test1 <- all(uniqueAllPixelsNotInCohortData %in% c(NA, 0L))
      if (!test1 | !test2 | !test3) {
        stop("Every value on pixelGroupMap greater than 0 must have a pixelIndex in cohorts.",
             " This test is failing, i.e., there are some pixelGroupMaps have pixelGroups, and aren't in cohorts.")
      }
    }

    # Bring to pixelGroup level -- this will squash the data.table
    allCohortData <- cohorts[ , .(ecoregionGroup = ecoregionGroup[1],
                                  mortality = mortality[1],
                                  aNPPAct = aNPPAct[1],
                                  sumB = sumB[1]),
                              by = uniqueCohortDefinition]

    theNewOnes <- is.na(allCohortData$B)
    cohortData <- allCohortData[!theNewOnes]
    newCohortData <- allCohortData[theNewOnes]

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(cohorts[, c("pixelIndex", "pixelGroup")],
                             by = c("pixelIndex"))

  }

  # update pixelGroupMap
  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (isTRUE(getOption("LandR.assertions"))) {
    if (!isTRUE(all(pixelsToChange$pixelGroup ==
                    pixelGroupMap[][pixelsToChange$pixelIndex])))
      stop("pixelGroupMap and newCohortData$pixelGroupMap don't match in updateCohortData fn")
  }

  ## give biomass in pixels that have serotiny/resprouting
  # newCohortData[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]

  ##########################################################
  # Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone)
  ##########################################################
  cohortData <- .initiateNewCohorts(newCohortData, cohortData,
                                    pixelGroupMap,
                                    time = time, speciesEcoregion = speciesEcoregion)

  outs <- rmMissingCohorts(cohortData, pixelGroupMap, firePixelTable)

  if (isTRUE(getOption("LandR.assertions"))) {
    maxPixelGroupFromCohortData <- max(outs$cohortData$pixelGroup)
    maxPixelGroup <- as.integer(maxValue(outs$pixelGroupMap))
    test1 <- (!identical(maxPixelGroup, maxPixelGroupFromCohortData))
    if (test1) {
      stop("The sim$pixelGroupMap and cohortData have unmatching pixelGroup.",
           " They must be matching.",
           " If this occurs, please contact the module developers")
    }
  }

  message(crayon::magenta("NUMBER OF UNIQUE PIXELGROUPS: ", length(unique(outs$cohortData$pixelGroup)),
                          " AND FORESTED PIXELS: ", sum(!is.na(outs$pixelGroupMap[]))))

  return(list(cohortData = outs$cohortData,
              pixelGroupMap = outs$pixelGroupMap))
}


#' \code{.initiateNewCohorts} will calculate new values for \code{B}, add
#' \code{age}, then \code{rbindlist} this with \code{cohortData}
#'
#' @inheritParams updateCohortData
#' @return
#' \code{.initiateNewCohorts} returns A \code{data.table} with a new,
#' \code{rbindlist}ed cohortData
#'
#' @rdname updateCohortData
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
.initiateNewCohorts <- function(newCohortData, cohortData, pixelGroupMap, time, speciesEcoregion) {
  ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum biomass per ecoregion, join to new cohort data
  namesNCD <- names(newCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newCohortData[, pixelGroup := getValues(pixelGroupMap)[pixelIndex]]
    } else {
      stop("newCohortData must have either pixelIndex or pixelGroup")
    }
  }

  specieseco_current <- speciesEcoregionLatestYear(speciesEcoregion, time)
  specieseco_current <- setkey(specieseco_current[, .(speciesCode, maxANPP, maxB, ecoregionGroup)],
                               speciesCode, ecoregionGroup)

  # Note that after the following join, some cohorts will be lost due to lack of
  #  parameters in speciesEcoregion. These need to be modified in pixelGroupMap.
  missingNewCohortData <- newCohortData[!specieseco_current, on = uniqueSpeciesEcoregionDefinition]
  specieseco_current <- specieseco_current[!is.na(maxB)]
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]
  newCohortData <- specieseco_current[newCohortData, on = uniqueSpeciesEcoregionDefinition]
  newCohortData <- newCohortData[!is.na(maxB)]
  # newCohortData <- newCohortData[specieseco_current, on = uniqueSpeciesEcoregionDefinition,
  #                                nomatch = 0]
  #newCohortData <- setkey(newCohortData, speciesCode, ecoregionGroup)[specieseco_current, nomatch = 0]
  set(newCohortData, NULL, "age", 1)  ## set age to 1
  set(newCohortData, NULL, "sumB", 0)
  ## set biomass - if B=0, it's getting maxANPP ???
  if ("B" %in% names(newCohortData))
    newCohortData[, B := NULL]
  set(newCohortData, NULL, "B",
      as.integer(pmax(1, newCohortData$maxANPP *
                        exp(-1.6 * newCohortData$sumB / newCohortData$maxB_eco))))
  set(newCohortData, NULL, "B", as.integer(pmin(newCohortData$maxANPP, newCohortData$B)))

  newCohortData <- newCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, age, B,
                                     mortality = 0L, aNPPAct = 0L, sumB = 0)]

  # This removes the duplicated pixels within pixelGroup, i.e., the reason we want pixelGroups
  newCohortData <- unique(newCohortData, by = uniqueCohortDefinition)

  cohortData <- rbindlist(list(cohortData, newCohortData), fill = TRUE, use.names = TRUE)
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
#' A \code{list} with 2 \code{data.table} objects, \code{cohortData} and \code{pixelGroupMap},
#' each updated based on missing pixelGroups in the other.
#'
#' @export
#' @importFrom data.table rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
rmMissingCohorts <- function(cohortData, pixelGroupMap, firePixelTable) {
  pgmValues <- data.table(pixelGroup = getValues(pixelGroupMap),
                          pixelIndex = seq(ncell(pixelGroupMap)))

  pgmVals <- na.omit(pgmValues)
  pgmVals <- pgmVals[pixelGroup > 0]
  whPgsStillInCDGoneFromPGM <- !cohortData$pixelGroup %in% pgmVals$pixelGroup
  pgsStillInCDGoneFromPGM <- cohortData[whPgsStillInCDGoneFromPGM,] #setdiff(unique(cohortData$pixelGroup), unique(pgmVals$pixelGroup))
  #pgsStillInPGMGoneFromCD <- setdiff(unique(pgmVals$pixelGroup), unique(cohortData$pixelGroup))
  whPgsStillInPGMGoneFromCD <- !pgmVals$pixelGroup %in% cohortData$pixelGroup
  pgsStillInPGMGoneFromCD <- pgmVals[whPgsStillInPGMGoneFromCD,]

  # REMOVE lines in cohortData that are no longer in the pixelGroupMap
  cohortData <- cohortData[!pixelGroup %in% pgsStillInCDGoneFromPGM$pixelGroup]
  # REMOVE pixels in pixelGroupMap that are no longer in the cohortData
  pixelGroupMap[pgsStillInPGMGoneFromCD$pixelIndex] <- NA

  if (isTRUE(getOption("LandR.assertions"))) {
    testCohortData(cohortData, pixelGroupMap, message = "rmMissingCohorts")
  }

  if (NROW(unique(cohortData[pixelGroup == 67724]$ecoregionGroup)) > 1) stop()

  return(list(cohortData = cohortData,
              pixelGroupMap = pixelGroupMap))
}


#' Create the correct string for pixelGroups
#'
#' @inheritParams generatePixelGroups
#' @param ecoregionGroup  A vector of ecoregionGroup strings
#' @param speciesGroup  A vector of speciesGroup strings
#'
#' @return  TODO: description needed
#'
.makePixelGroups <- function(maxPixelGroup, ecoregionGroup, speciesGroup,
                             columns = c("ecoregionGroup", "speciesGroup", "age")) {
  as.integer(maxPixelGroup) +
    as.integer(factor(paste(ecoregionGroup, speciesGroup, sep = "_")))
}

#' Add the correct \code{pixelGroups} to a \code{pixelCohortData} object
#'
#' @param pixelCohortData  # pixel groups are groups of identical pixels based
#'   on \code{speciesGroup} x \code{Age} and \code{ecoregionGroup}.
#' @param maxPixelGroup A length 1 numeric/integer indicating the current maximum pixelGroup value
#' @param columns A character vector of column names to use as part of the generation of unique
#'   combinations of features. Default is \code{c("ecoregionGroup", "speciesCode", "age", "biomass")}
#'
#' @return
#' Returns a vector of pixelGroup in the original order of the input \code{pixelCohortData}.
#' This should likely be added to the \code{pixelCohortData} object immediately.
#'
#'
#' @export
#' @importFrom data.table setkey
#' @importFrom SpaDES.core paddedFloatToChar
generatePixelGroups <- function(pixelCohortData, maxPixelGroup,
                                columns = c("ecoregionGroup", "speciesCode", "age", "biomass")) {

  columnsOrig <- columns
  columns <- columns[columns %in% names(pixelCohortData)]
  columns2 <- paste0(columns, "2")
  if (!all(columns == columnsOrig))
    message("Creating pixelGroup values, but not using all columns requested. Only using, ",
            paste(columns, collapse = ", "), " instead of ", paste(columnsOrig, collapse = ", "))

  pcd <- pixelCohortData # no copy -- just for simpler name

  # Convert to unique numeric
  pcd[ , c(columns2) := lapply(.SD, function(x) {
    a <- as.integer(factor(x))
  }), .SDcols = columns]

  # concatenate within rows -- e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
  pcd[, uniqueComboByRow :=
        as.integer(factor(do.call(paste, append(list(sep = "_"), as.list(.SD))))),
      .SDcols = columns2]

  # concatenate within pixelIndex
  pcd[ , c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"), by = "pixelIndex"]
  pcd[ , c("pixelGroup") := as.integer(maxPixelGroup) + as.integer(factor(uniqueComboByPixelIndex))]

  # clean up
  pcd[, c(columns2, "uniqueComboByPixelIndex", "uniqueComboByRow") := (NULL)]
  return(pcd$pixelGroup)
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

#' A test that pixelGroupMap and cohortData match
#'
#' @inheritParams updateCohortData
#' @param sim If the simList is included, then the browser() call will be more useful
#' @param maxExpectedNumDiverge A numeric, length 1, indicating by how many they
#'   can diverge. Default 1.
#' @param message An optional message to print. This may help identify where this function
#'   was called.
#' @note
#' TODO
#'
#' @export
#' @importFrom crayon green
#' @importFrom stats na.omit
testCohortData <- function(cohortData, pixelGroupMap, sim, maxExpectedNumDiverge = 1,
                           message = "") {
  browser(expr = exists("aaaa"))

  a <- sort(unique(na.omit(pixelGroupMap[])))
  b <- sort(unique(na.omit(cohortData$pixelGroup)))
  test1 <- sum(!a %in% b)  # can be 1 because there could be pixelGroup of 0, which is OK to not match
  test2 <- sum(!b %in% a)  # can be 1 because there could be pixelGroup of 0, which is OK to not match
  if (test1 > maxExpectedNumDiverge || test2 > maxExpectedNumDiverge) {
    if (nchar(message) > 0) message(message)
    if (test1 > maxExpectedNumDiverge) message("test1 is ", test1, " -- too many pixelGroups on pixelGroupMap")
    if (test2 > maxExpectedNumDiverge) message("test2 is ", test2, " -- too many pixelGroups in cohortData")
    stop("The sim$pixelGroupMap and cohortData have unmatching pixelGroup. They must be matching.",
         " Please contact the module developers")
  } else {
    message(crayon::green("  -- assertion passed using testCohortData --"))
  }
}

.ageRndUpSuccessionTimestep <- function(age, successionTimestep) {
  as.integer(ceiling(as.numeric(age) / successionTimestep) * successionTimestep)
}

#' The columns in a cohortData that define "unique"
#'
#' If 2 pixels have identical values in all of these columns, then they are the
#' same \code{pixelGroup}
#' @export
#' @rdname uniqueDefinitions
uniqueCohortDefinition <- c("pixelGroup", "speciesCode", "age", "B")

#' The columns in a cohortData that define "unique"
#'
#' If 2 pixels have identical values in all of these columns, then they are the
#' same \code{pixelGroup}
#' @export
#' @rdname uniqueDefinitions
uniqueSpeciesEcoregionDefinition <- c("speciesCode", "ecoregionGroup")


#' Summary for cohortData
#' @param cohortData A cohortData object
#' @export
describeCohortData <- function(cohortData) {
  vals <- c("biomass", "totalBiomass", "age", "cover")
  names(vals) <- vals
  out <- lapply(vals, function(val)
    .cohortMessages(cohortData, val)
  )
  message(magenta("Pixels with non-NA cover:, ", cohortData[!is.na(cover), length(unique(pixelIndex))]))
}

.cohortMessages <- function(cohortData, val) {
  out <- list()
  if (val %in% colnames(cohortData)) {
    pixelsNA <- NROW(cohortData[is.na(get(val)), unique("pixelIndex"), with = FALSE])
    message(magenta("Pixels with missing", val, ":", format(pixelsNA, big.mark = ",")))
    pixelsZero <- NROW(cohortData[, all(get(val) == 0), by = "pixelIndex"][get("V1") ==TRUE])
    message(magenta("Pixels with all(",val," == 0): ", format(pixelsZero, big.mark = ",")))
    pixelsBiomassNonZero <- NROW(cohortData[, any(get(val) > 0), by = "pixelIndex"][get("V1") ==TRUE])
    message(magenta("Pixels with all(",val," > 0): ", format(pixelsBiomassNonZero, big.mark = ",")))
    out <- list(pixelsNA = pixelsNA, pixelsZero = pixelsZero, pixelsBiomassNonZero = pixelsBiomassNonZero)
  }
  return(invisible(out))
}
