if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(
    ".", ".I", ":=", "age", "aNPPAct", "columnsForPG", "cover", "coverOrig",
    "ecoregion", "ecoregionGroup", "hasBadAge",
    "imputedAge", "initialEcoregion", "initialEcoregionCode", "initialPixels",
    "lcc", "maxANPP", "maxB", "maxB_eco", "mortality",
    "newPossLCC", "outBiomass", "pixelIndex", "pixels",
    "speciesposition", "speciesGroup", "speciesInt", "sumB",
    "temppixelGroup", "totalBiomass",
    "uniqueCombo", "uniqueComboByRow", "uniqueComboByPixelIndex", "V1", "year"
  ))
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
#'   \item assign initial B and age for new cohort;
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
updateCohortData <- function(newCohortData, cohortData, pixelGroupMap, time,
                             speciesEcoregion, firePixelTable = NULL,
                             successionTimestep) {

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
    message(crayon::green("  Regenerating open and pixels with B (likely after seed dispersal)"))

    allNewPixelGroups <- FALSE

    pixelIndex <- which(pixelGroupMap[] %in% cohortData$pixelGroup)
    cohortDataPixelIndex <- data.table(pixelIndex = pixelIndex,
                                       pixelGroup = pixelGroupMap[][pixelIndex])
    cdLong <- cohortDataPixelIndex[cohortData, on = "pixelGroup", allow.cartesian = TRUE]
    cohorts <- rbindlist(list(cdLong, newCohortData), use.names = TRUE, fill = TRUE)

    columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B")
    cd <- cohorts[,c("pixelIndex", columnsForPG), with = FALSE]
    cohorts[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L,
                                                columns = columnsForPG)]

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

  ## give B in pixels that have serotiny/resprouting
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

  message(crayon::magenta("NUMBER OF UNIQUE PIXELGROUPS:", length(unique(outs$cohortData$pixelGroup)),
                          ", FORESTED PIXELS:", sum(!is.na(outs$pixelGroupMap[])),
                          ", PIXELS WITH NO PIXEL GROUP:", sum(outs$pixelGroupMap[]==0, na.rm = TRUE)))

  return(list(cohortData = outs$cohortData,
              pixelGroupMap = outs$pixelGroupMap))
}

#' Initiate new cohorts
#'
#' Calculate new values for \code{B}, add \code{age}, then \code{rbindlist} this
#' with \code{cohortData}.
#'
#' @inheritParams updateCohortData
#' @return
#' \code{.initiateNewCohorts} returns A \code{data.table} with a new,
#' \code{rbindlist}ed cohortData
#'
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
#' @rdname updateCohortData
.initiateNewCohorts <- function(newCohortData, cohortData, pixelGroupMap, time,
                                speciesEcoregion) {
  ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum B per ecoregion, join to new cohort data
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
  ## set B - if B=0, it's getting maxANPP ???
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
#' @param cohortData A \code{data.table} with columns:
#'   \code{pixelGroup}, \code{ecoregionGroup}, \code{speciesCode}, \code{age},
#'   \code{B}, \code{mortality}, \code{aNPPAct}, ond \code{sumB}.
#' @param pixelGroupMap Raster layer with pixel values equal to a pixel group number
#'   that correspondsd exactly to ]\code{pixelGroup} column in \code{cohortData}.
#' @param firePixelTable A data.table with 2 columns, \code{pixelIndex} and \code{pixelGroup}.
#'   This will be used in conjunction with \code{cohortData} and \code{pixelGroupMap}
#'   to ensure that everything matches correctly.
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
#'   combinations of features. Default is \code{c("ecoregionGroup", "speciesCode", "age", "B")}
#'
#' @return
#' Returns a vector of pixelGroup in the original order of the input \code{pixelCohortData}.
#' This should likely be added to the \code{pixelCohortData} object immediately.
#'
#' @export
#' @importFrom data.table setkey
#' @importFrom SpaDES.core paddedFloatToChar
generatePixelGroups <- function(pixelCohortData, maxPixelGroup,
                                columns = c("ecoregionGroup", "speciesCode",
                                            "age", "B")) {

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
  vals <- c("B", "totalBiomass", "age", "cover")
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


#' Convert Land Cover Classes to another value in is neighbourhood
#'
#' This will search around the pixels on \code{rstLCC} that have
#' \code{pixelClassesToReplace}, and search in iteratively increasing
#' radii outwards for other Land Cover Classes than \code{pixelClassesToReplace}.
#' It will then take the cohorts that were in pixels with \code{pixelClassesToReplace}
#' and assign them new values in the output object. This function will
#' also check that it must be an \code{ecoregionCode} that already exists in
#' \code{cohortData}, i.e., not create new \code{ecoregionCode} values.
#'
#' @param pixelClassesToReplace Integer vector of classes that are are to be replaced, e.g.,
#'      34, 35, 36 on LCC2005, which are burned young, burned 10yr, and cities
#'
#' @param rstLCC LCC raster, e.g., LCC2005
#'
#' @param ecoregionGroupVec The vector of ecoregionGroup codes
#'
#' @param pixelCohortData A \code{data.table} with individual cohorts, with data
#'     for every pixel, columns: \code{initialEcoregionCode}, \code{speciesCode}.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom data.table rbindlist setnames
#' @importFrom raster raster
#' @importFrom SpaDES.core paddedFloatToChar
#' @importFrom SpaDES.tools spread2
convertUnwantedLCC <- function(pixelClassesToReplace = 34:36, rstLCC,
                               ecoregionGroupVec, speciesEcoregion,
                               availableERC_by_Sp) {
  rstUnwantedLCC <- integer(length(ecoregionGroupVec))
  rstUnwantedLCC[] <- NA;
  rstUnwantedLCC[gsub(".*_", "", ecoregionGroupVec) %in% pixelClassesToReplace] <- 1
  theUnwantedPixels <- which(rstUnwantedLCC == 1)
  theUnwantedPixels <- theUnwantedPixels[theUnwantedPixels %in% availableERC_by_Sp$pixelIndex]
  #cdEcoregionCodes <- as.character(unique(pixelCohortData$initialEcoregionCode))
  #availableERC_by_Sp <- unique(pixelCohortData, by = c("initialEcoregionCode", "speciesCode", "B"))
  #availableERC_by_Sp <- availableERC_by_Sp[eval(rowsInPCDToKeep)]
  #availableERC_by_Sp <- unique(availableERC_by_Sp, by = c("initialEcoregionCode", "speciesCode"))
  #availableERC_by_Sp[, `:=`(age = NULL, logAge = NULL, totalBiomass = NULL,
  #                          cover = NULL, coverOrig = NULL, B = NULL, lcc = NULL)]
  if (getOption("LandR.assertions")) {
    #theUnwantedCellsFromCD <- unique(pixelCohortData[lcc %in% pixelClassesToReplace]$pixelIndex)
    #iden <- identical(sort(unique(theUnwantedPixels)), sort(theUnwantedCellsFromCD))
    #if (!iden)
    #  stop("values of 34 and 35 on pixelCohortData and sim$LCC2005 don't match")
  }
  iterations <- 1
  # cd <- pixelCohortData[, .(pixelIndex, initialEcoregionCode, speciesCode)]
  numCharEcoregion <- nchar(gsub("_.*", "", availableERC_by_Sp$initialEcoregionCode[1]))
  while (length(theUnwantedPixels) > 0) {
    out <- spread2(rstLCC, start = theUnwantedPixels, asRaster = FALSE,
                   iterations = iterations, allowOverlap = TRUE, spreadProb = 1)
    out <- out[initialPixels != pixels] # rm pixels which are same as initialPixels --> these are known wrong
    iterations <- iterations + 1
    out[, lcc := rstLCC[][pixels]]
    # out <- unique(availableERC_by_Sp)[out, on = c("pixelIndex" = "pixels"), nomatch = 0] # join the availableERC_by_Sp which has initialEcoregionCode
    out[lcc %in% c(pixelClassesToReplace), lcc:=NA]
    out <- na.omit(out)
    # out[, `:=`(potentialNewERG = ecoregionGroupVec[pixels])]
    # attach speciesCode -- need this so that we know which speciesCodes to test against speciesEcoregion
    a <- out[availableERC_by_Sp, allow.cartesian = TRUE, on = c("initialPixels" = "pixelIndex"), nomatch = NA]
    out5 <- availableERC_by_Sp[out, allow.cartesian = TRUE,
                               on = c("pixelIndex" = "initialPixels"), nomatch = NA] # join the availableERC_by_Sp which has initialEcoregionCode

    # Any that don't have all species, should be removed
    out5[, toDelete := any(is.na(speciesCode)), by = "pixelIndex"]
    out6 <- out5[toDelete == FALSE]
    # out7 <- out5[toDelete == TRUE]

    # browser(expr = out6[pixelIndex == 194455])
    rowsToKeep <- out6[, list(keep = .resample(.I, 1)), by = c("pixelIndex")] # random sample of available, weighted by abundance
    out2 <- out6[rowsToKeep$keep,
                 list(newPossLCC = lcc,
                      initialEcoregion = substr(initialEcoregionCode, 1, numCharEcoregion),
                      pixelIndex)]
    #out2 <- out2[rowsToKeep$keep]
    out2[, ecoregionGroup := paste0(initialEcoregion, "_",
                                    paddedFloatToChar(as.numeric(newPossLCC), padL = 2, padR = 0))]
    out2[, initialEcoregion := NULL]

    # remove combinations of ecoregionGroup and speciesCode that don't exist -- Now this excludes B = 0
    #out2 <- availableERC_by_Sp[out2, on = c("pixelIndex" = "initialPixels", "initialEcoregionCode" = "ecoregionGroup", "speciesCode"),
    #                           nomatch = NA]
    #hasMatch <- out2[, all(!is.na(rasterToMatch)), by = c("pixelIndex", "initialEcoregionCode", "speciesCode")][V1 == TRUE]
    keepPixels <- unique(out2$pixelIndex)
    # out2 <- out2[ecoregionGroup %in% cdEcoregionCodes] # remove codes that don't exist in pixelCohortData

    theUnwantedPixels <- theUnwantedPixels[!theUnwantedPixels %in% keepPixels]
    out2 <- unique(out2)
    if (!exists("out3")) {
      out3 <- out2
    } else {
      out3 <- rbindlist(list(out2, out3))
    }
  }

  # setnames(out3, c("initialPixels", "initialEcoregionCode"), c("pixelIndex", "ecoregionGroup"))
  out3[, `:=`(newPossLCC = NULL)]
  # out3 <- unique(out3, by = c("pixelIndex", "ecoregionGroup"))

  out3
}

#' Generate initial \code{cohortData} table
#'
#' Takes a single \code{data.table} input, which has the following columns in addition to
#' others that will be labelled with species name, and contain percent cover of each:
#' \itemize{
#'   \item age
#'   \item logAge
#'   \item initialEcoregionCode
#'   \item totalBiomass
#'   \item pixelIndex
#'   \item lcc
#' }
#'
#' @param inputDataTable A \code{data.table} with columns described above.
#' @param sppColumns A vector of the names of the columns in \code{inputDataTable} that
#'   represent percent cover by species
#' @param pixelGroupBiomassClass Round B to the nearest \code{pixelGroupBiomassClass}
#'   to establish unique pixelGroups
#'
#' @author Eliot McIntire
#' @export
#' @importFrom crayon blue
#' @importFrom data.table melt setnames
#' @importFrom reproducible Cache
makeAndCleanInitialCohortData <- function(inputDataTable, sppColumns, pixelGroupBiomassClass) {
  ### Create groupings
  if (isTRUE(getOption("LandR.assertions"))) {
    expectedColNames <- c("age", "logAge", "initialEcoregionCode", "totalBiomass",
                          "lcc", "pixelIndex")
    if (!all(expectedColNames %in% colnames(inputDataTable)))
      stop("Column names for inputDataTable must include ", expectedColNames)
    if (!all(sppColumns %in% colnames(inputDataTable)))
      stop("Species names are incorrect")
    if (!all(unlist(lapply(inputDataTable[, sppColumns, with = FALSE],
                           function(x) all(x >= 0 & x <= 100)  ))))
      stop("Species columns are not percent cover between 0 and 100. This may",
           " be because they more NA values than the Land Cover raster")

  }
  coverColNames <- grep(colnames(inputDataTable), pattern = "cover", value = TRUE)
  newCoverColNames <- gsub("cover\\.", "", coverColNames)
  setnames(inputDataTable, old = coverColNames, new = newCoverColNames)
  message(crayon::blue("Create initial cohortData object, with no pixelGroups yet"))

  cohortData <- data.table::melt(inputDataTable,
                                 value.name = "cover",
                                 measure.vars = newCoverColNames,
                                 variable.name = "speciesCode")
  cohortData[, coverOrig := cover]

  if (getOption("LandR.assertions"))
    #describeCohortData(cohortData)
    message(blue("assign B = 0 and age = 0 for pixels where cover = 0, ",
                 "\n  because cover is most reliable dataset"))
  cohortData[cover == 0, `:=`(age = 0L, B = 0L)]
  message(blue("assign totalBiomass = 0 sum(cover) = 0 in a pixel, ",
               "\n  because cover is most reliable dataset"))
  cohortData <- cohortData[, sum(cover) == 0, by = "pixelIndex"][V1 == TRUE][
    cohortData, on = "pixelIndex"][V1 == TRUE, totalBiomass := 0L]
  cohortData[, V1 := NULL]

  ######################
  # message(crayon::blue("POSSIBLE ALERT -- assume deciduous cover is 1/2 the conversion to B as conifer"))
  # cohortData[speciesCode == "Popu_sp", cover := asInteger(cover / 2)] # CRAZY TODO -- DIVIDE THE COVER BY 2 for DECIDUOUS -- will only affect mixed stands
  cohortData[ , cover := {
    sumCover <- sum(cover)
    if (sumCover > 100) {
      cover <- asInteger(cover/(sumCover + 0.0001) * 100L)
    }
    cover
  }, by = "pixelIndex"]

  # Biomass -- by cohort
  message(crayon::blue("Divide total B of each pixel by the relative cover of the cohorts"))
  cohortData[ , B := asInteger(mean(totalBiomass) * cover / 100), by = "pixelIndex"] # /100 because cover is percent
  message(crayon::blue("Round B to nearest P(sim)$pixelGroupBiomassClass"))
  cohortData[ , B := asInteger(ceiling(B / pixelGroupBiomassClass) *
                                       pixelGroupBiomassClass)] # /100 because cover is percent
  message(blue("Set B to 0 where cover > 0 and age = 0, because B is least quality dataset"))
  cohortData[ , totalBiomass := asInteger(totalBiomass)]

  ######################################################
  # Impute missing ages on poor age dataset
  ######################################################
  cohortDataMissingAge <- cohortData[, hasBadAge := all(age == 0 & cover > 0) | any(is.na(age)), by = "pixelIndex"][
    hasBadAge == TRUE]
  cohortDataMissingAgeUnique <- unique(cohortDataMissingAge,
                                       by = c("initialEcoregionCode", "speciesCode"))[
                                         , .(initialEcoregionCode, speciesCode)]
  cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[cohortData, on = c("initialEcoregionCode", "speciesCode"), nomatch = 0]
  ageQuotedFormula <- quote(age ~ B * speciesCode + (1 | initialEcoregionCode) + cover)
  cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[, .(B, age, speciesCode, initialEcoregionCode, cover)]
  message(blue("Impute missing age values"))
  system.time(outAge <- Cache(statsModel, form = ageQuotedFormula, .specialData = cohortDataMissingAgeUnique))

  print(outAge$rsq)

  cohortDataMissingAge[, imputedAge := pmax(0L, asInteger(predict(outAge$mod, newdata = cohortDataMissingAge)))]
  cohortData <- cohortDataMissingAge[, .(pixelIndex, imputedAge, speciesCode)][cohortData, on = c("pixelIndex", "speciesCode")]
  cohortData[!is.na(imputedAge), `:=`(age = imputedAge, logAge = log(imputedAge))]
  cohortData[, `:=`(imputedAge = NULL, hasBadAge = NULL)]

  #######################################################
  # set B to zero if age is zero because B is lowest quality dataset
  #######################################################
  message(blue("Set recalculate totalBiomass as sum(B); many biomasses will have been set to 0 in previous steps"))
  cohortData[cover > 0 & age == 0, B := 0L]
  cohortData[, totalBiomass := sum(B), by = "pixelIndex"]

  # This was unused, but for beta regression, this can allow 0s and 1s without needing a separate
  #   model for zeros and ones
  # https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
  # cohortData[ , coverProp := (cover/100 * (NROW(cohortData) - 1) + 0.5) / NROW(cohortData)]

  cohortData
}

#' The generic statistical model -- to run lmer or glmer
#'
#' This does a few things including R squared, gets the fitted values.
#' It appears that running the models "asis" without this wrapper
#' does not work with \code{Cache}. The return of the model in
#' a list solves this problem.
#'
#' @param form A quoted formula to test
#' @param .specialData The custom dataset required for the model
#' @param ... Anything passed to args for the model
#'
#' @export
#' @importFrom lme4 glmer lmer
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom stats fitted predict
statsModel <- function(form, .specialData, ...) {
  if ("family" %in% names(list(...))) {
    modelFn <- lme4::glmer
  } else {
    modelFn <- lme4::lmer
  }
  mod <- modelFn(
    formula = eval(form),
    data = .specialData,
    ...)

  list(mod = mod, pred = fitted(mod), rsq = MuMIn::r.squaredGLMM(mod))
}

#' Default columns that define pixel groups
#'
#' @export
columnsForPixelGroups <- c("ecoregionGroup", "speciesCode", "age", "B")
