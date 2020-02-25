if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(
    ".", ".I", ":=", "..groupVar", "age", "aNPPAct", "cover", "coverOrig",
    "ecoregion", "ecoregionGroup", "hasBadAge",
    "imputedAge", "initialEcoregion", "initialEcoregionCode", "initialPixels",
    "lcc", "maxANPP", "maxB", "maxB_eco", "mortality",
    "newPossLCC", "noPixels", "oldSumB", "ord", "outBiomass", "oldEcoregionGroup",
    "pixelGroup2", "pixelIndex", "pixels", "possERC",
    "speciesposition", "speciesGroup", "speciesInt", "state", "sumB",
    "temppixelGroup", "toDelete", "totalBiomass",
    "uniqueCombo", "uniqueComboByRow", "uniqueComboByPixelIndex", "V1", "year"
  ))
}

#' Add cohorts to \code{cohortData} and \code{pixelGroupMap}
#'
#' This is a wrapper for  \code{generatePixelGroups}, \code{initiateNewCohort} and updates to
#' \code{pixelGroupMap} via assignment to new \code{pixelIndex} values in \code{newPixelCohortData}.
#' By running these all together, there is less chance that they will diverge.
#' There are some checks internally for consistency.
#'
#' Does the following:
#' \enumerate{
#'   \item add new cohort data into \code{cohortData};
#'   \item assign initial \code{B} and \code{age} for new cohort;
#'   \item assign the new \code{pixelGroup} to the pixels that have new cohort;
#'   \item update the \code{pixelGroup} map.
#' }
#'
#' @template newPixelCohortData
#' @template cohortData
#' @template pixelGroupMap
#' @template currentTime
#' @template speciesEcoregion
#'
#' @param treedFirePixelTableSinceLastDisp A data.table with at least 2 columns, \code{pixelIndex} and \code{pixelGroup}.
#'   This will be used in conjunction with \code{cohortData} and \code{pixelGroupMap}
#'   to ensure that everything matches correctly.
#' @param provenanceTable table with pixel-based provenances for reforestation
#' @param successionTimestep The time between successive seed dispersal events.
#'   In LANDIS-II, this is called "Succession Timestep". This is used here
#'
#' @param verbose Integer, where increasing number is increasing verbosity. Currently,
#'    only level 1 exists; but this may change.
#'
#' @template doAssertion
#'
#' @return
#' A list of length 2, \code{cohortData} and \code{pixelGroupMap}, with
#' \code{newPixelCohortData} inserted.
#'
#' @export
#' @rdname updateCohortData
#' @importFrom crayon green magenta
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom SpaDES.core paddedFloatToChar
#' @importFrom stats na.omit
updateCohortData <- function(newPixelCohortData, cohortData, pixelGroupMap, currentTime,
                             speciesEcoregion, treedFirePixelTableSinceLastDisp = NULL,
                             provenanceTable = NULL,
                             successionTimestep,
                             verbose = getOption("LandR.verbose", TRUE),
                             doAssertion = getOption("LandR.assertions", TRUE)) {

  maxPixelGroup <- as.integer(maxValue(pixelGroupMap))

  if (!is.null(treedFirePixelTableSinceLastDisp)) {
    ## only set to 0 the pixels that became empty (have no survivors)
    emptyPix <- setdiff(treedFirePixelTableSinceLastDisp$pixelIndex,
                        newPixelCohortData[type == "survivor", pixelIndex])
    pixelGroupMap[emptyPix] <- 0L
  }
  relevantPixels <- pixelGroupMap[][newPixelCohortData$pixelIndex]
  zeroOnPixelGroupMap <- relevantPixels == 0

  if (!"age" %in% colnames(newPixelCohortData))
    newPixelCohortData[, age := 1L]

  if (all(zeroOnPixelGroupMap)) {
    # Deal with pixels on the map that have no pixelGroup -- these are burned
    # pixels --> the entirely newly regenerated pixels do not require a
    # re-pixelGroupMaping  -- can just add to existing pixelGroup values
    if (verbose > 0)
      message(crayon::green("  Regenerating only disturbed pixels (i.e. resprouting & serotiny, in addition to surviving cohorts)"))
    columnsForPG <- c("ecoregionGroup", "speciesCode", "age") # no Biomass because they all have zero
    cd <- newPixelCohortData[,c("pixelIndex", columnsForPG), with = FALSE]
    newPixelCohortData[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = maxPixelGroup,
                                                           columns = columnsForPG)]#,
    #successionTimestep = successionTimestep)

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(newPixelCohortData[, c("pixelIndex", "pixelGroup")],
                             by = c("pixelIndex"))
  } else {
    # This is for situations where there are some empty pixels being filled,
    #   and some occupied pixels getting infilling. This requires a wholesale
    #   re-pixelGroup
    if (verbose > 0)
      message(crayon::green("  Regenerating open and pixels with B (likely after seed dispersal, or partial mortality following disturbance)"))

    pixelIndex <- which(pixelGroupMap[] %in% cohortData$pixelGroup)
    cohortDataPixelIndex <- data.table(pixelIndex = pixelIndex,
                                       pixelGroup = pixelGroupMap[][pixelIndex])
    cdLong <- cohortDataPixelIndex[cohortData, on = "pixelGroup", allow.cartesian = TRUE]
    cohorts <- rbindlist(list(cdLong, newPixelCohortData), use.names = TRUE, fill = TRUE)

    columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B")
    cd <- cohorts[, c("pixelIndex", columnsForPG), with = FALSE]
    cohorts[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L, columns = columnsForPG)]

    # Bring to pixelGroup level -- this will squash the data.table
    if (is.null(cohorts[["sumB"]])) {
      cohorts[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]
    }
    allCohortData <- cohorts[ , .(ecoregionGroup = ecoregionGroup[1],
                                  mortality = mortality[1],
                                  aNPPAct = aNPPAct[1],
                                  sumB = sumB[1]),
                              by = uniqueCohortDefinition]

    theNewOnes <- is.na(allCohortData$B)
    cohortData <- allCohortData[!theNewOnes]
    newPixelCohortData <- allCohortData[theNewOnes]

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(cohorts[, c("pixelIndex", "pixelGroup")], by = c("pixelIndex"))
  }

  # update pixelGroupMap
  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (doAssertion) {
    if (!isTRUE(all(pixelsToChange$pixelGroup == pixelGroupMap[][pixelsToChange$pixelIndex])))
      stop("pixelGroupMap and newPixelCohortData$pixelGroupMap don't match in updateCohortData fn")
  }

  ## give B in pixels that have serotiny/resprouting
  # newPixelCohortData[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]

  ##########################################################
  # Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone)
  ##########################################################
  cohortData <- .initiateNewCohorts(newPixelCohortData, cohortData,
                                      pixelGroupMap, currentTime = currentTime,
                                      speciesEcoregion = speciesEcoregion,
                                      successionTimestep = successionTimestep)

  outs <- rmMissingCohorts(cohortData, pixelGroupMap)

  assertCohortData(outs$cohortData, outs$pixelGroupMap,
                   doAssertion = doAssertion, verbose = verbose)

  if (doAssertion) {
    maxPixelGroupFromCohortData <- max(outs$cohortData$pixelGroup)
    maxPixelGroup <- as.integer(maxValue(outs$pixelGroupMap))
    test1 <- (!identical(maxPixelGroup, maxPixelGroupFromCohortData))
    if (test1) {
      stop("The sim$pixelGroupMap and cohortData have unmatching pixelGroup.",
           " They must be matching.",
           " If this occurs, please contact the module developers")
    }
  }

  if (verbose > 0) {
    nPixForest <- sum(!is.na(outs$pixelGroupMap[]))
    nPixGrps <- length(unique(outs$cohortData$pixelGroup))
    nPixNoPixGrp <- sum(outs$pixelGroupMap[] == 0, na.rm = TRUE)
    nPixTreed <- sum(outs$pixelGroupMap[] != 0, na.rm = TRUE)

    nDigits <- max(nchar(c(nPixForest, nPixGrps, nPixNoPixGrp))) + 3
    message(crayon::magenta("NUMBER OF FORESTED PIXELS          :",
                            paddedFloatToChar(nPixForest, padL = nDigits, pad = " ")))
    message(crayon::magenta("NUMBER OF PIXELS WITH TREES        :",
                            paddedFloatToChar(nPixTreed, padL = nDigits, pad = " ")))
    message(crayon::magenta("NUMBER OF UNIQUE PIXELGROUPS       :",
                            paddedFloatToChar(nPixGrps, padL = nDigits, pad = " ")))
    message(crayon::magenta("NUMBER OF PIXELS WITH NO PIXELGROUP:",
                            paddedFloatToChar(nPixNoPixGrp, padL = nDigits, pad = " ")))
  }

  return(list(cohortData = outs$cohortData,
              pixelGroupMap = outs$pixelGroupMap))
}

#' Initiate new cohorts
#'
#' Calculate new values for \code{B}, add \code{age}, then \code{rbindlist} this
#' with \code{cohortData}.
#'
#' @template newPixelCohortData
#' @template cohortData
#'
#' @return A \code{data.table} with a new \code{rbindlist}ed \code{cohortData}
#'
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
#' @rdname updateCohortData
.initiateNewCohorts <- function(newPixelCohortData, cohortData, pixelGroupMap, currentTime,
                                speciesEcoregion, successionTimestep) {
  ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum B per ecoregion, join to new cohort data
  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, pixelGroup := getValues(pixelGroupMap)[pixelIndex]]
    } else {
      stop("newPixelCohortData must have either pixelIndex or pixelGroup")
    }
  }

  ## simplifying to pixelGroup again
  if (!is.null(newPixelCohortData[["pixelIndex"]]))
    set(newPixelCohortData, NULL, "pixelIndex", NULL)
  newPixelCohortData <- newPixelCohortData[!duplicated(newPixelCohortData), ]   ## faster than unique

  specieseco_current <- speciesEcoregionLatestYear(speciesEcoregion, currentTime)
  specieseco_current <- setkey(specieseco_current[, .(speciesCode, maxANPP, maxB, ecoregionGroup)],
                               speciesCode, ecoregionGroup)

  ## Note that after the following join, some cohorts will be lost due to lack of
  ##  parameters in speciesEcoregion. These need to be modified in pixelGroupMap.
  #missingNewPixelCohortData <- newPixelCohortData[!specieseco_current, on = uniqueSpeciesEcoregionDefinition]
  specieseco_current <- specieseco_current[!is.na(maxB)]
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]
  newPixelCohortData <- specieseco_current[newPixelCohortData, on = uniqueSpeciesEcoregionDefinition]
  newPixelCohortData <- newPixelCohortData[!is.na(maxB)]

  if (any(newPixelCohortData$age > 1))
    stop("newPixelCohortData should only have new cohorts aged 1")
  set(newPixelCohortData, NULL, "age", 1L)  ## set age to 1

  ## Ceres: this was causing new cohorts to be initialized with maxANPP.
  ## instead, calculate total biomass of older cohorts
  # set(newPixelCohortData, NULL, "sumB", 0L)
  if (!is.null(newPixelCohortData[["sumB"]]))
    set(newPixelCohortData, NULL, "sumB", NULL)
  cohortData[age >= successionTimestep, oldSumB := sum(B, na.rm = TRUE), by = "pixelGroup"]

  newPixelCohortData <- unique(cohortData[, .(pixelGroup, oldSumB)],
                               by = "pixelGroup")[newPixelCohortData, on = "pixelGroup"]
  set(newPixelCohortData, which(is.na(newPixelCohortData$oldSumB)), "oldSumB", 0)   ## faster than [:=]
  setnames(newPixelCohortData, "oldSumB", "sumB")
  set(cohortData, NULL, "oldSumB", NULL)

  if ("B" %in% names(newPixelCohortData))
    newPixelCohortData[, B := NULL]
  set(newPixelCohortData, NULL, "B",
      asInteger(pmax(1, newPixelCohortData$maxANPP *
                       exp(-1.6 * newPixelCohortData$sumB / newPixelCohortData$maxB_eco))))
  set(newPixelCohortData, NULL, "B", asInteger(pmin(newPixelCohortData$maxANPP, newPixelCohortData$B)))

  newPixelCohortData <- newPixelCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, age, B,
                                               mortality = 0L, aNPPAct = 0L)]

  if (getOption("LandR.assertions")) {
    if (isTRUE(NROW(unique(newPixelCohortData, by = uniqueCohortDefinition)) != NROW(newPixelCohortData)))
      stop("Duplicated new cohorts in a pixelGroup. Please debug LandR:::.initiateNewCohorts")
  }

  cohortData <- rbindlist(list(cohortData, newPixelCohortData), fill = TRUE, use.names = TRUE)
  # cohortData[, sumB := sum(B, na.rm = TRUE), by = "pixelGroup"]  ## recalculate sumB
  # if (!is.integer(cohortData[["sumB"]]))
  #   set(cohortData, NULL, "sumB", asInteger(cohortData[["sumB"]]))

  return(cohortData)
}

#' Remove missing cohorts from \code{cohortData} based on \code{pixelGroupMap}
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template doAssertion
#'
#' @return
#' A \code{list} with 2 \code{data.table} objects, \code{cohortData} and \code{pixelGroupMap},
#' each updated based on missing \code{pixelGroups} in the other.
#'
#' @export
#' @importFrom data.table rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
rmMissingCohorts <- function(cohortData, pixelGroupMap,
                             doAssertion = getOption("LandR.assertions", TRUE)) {

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

  assertCohortData(cohortData, pixelGroupMap, message = "rmMissingCohorts",
                   doAssertion = doAssertion)

  if (NROW(unique(cohortData[pixelGroup == 67724]$ecoregionGroup)) > 1) stop()

  return(list(cohortData = cohortData,
              pixelGroupMap = pixelGroupMap))
}

#' Add the correct \code{pixelGroups} to a \code{pixelDataTable} object
#'
#' Generates unique groupings of a \code{data.table} object where one or more rows can
#' all belong to the same \code{pixelIndex}. Pixel groups will be identical pixels based
#' on unique combinations of \code{columns}.
#'
#' @param pixelDataTable  A \code{data.table} with column-based descriptions.
#'   Must have a column called \code{pixelIndex}, which allows for multiple rows to be associated
#'   with a single pixel.
#' @param maxPixelGroup A length 1 numeric indicating the current maximum \code{pixelGroup} value;
#'    the \code{pixelGroup} numbers returned will start at \code{maxPixelGroup + 1}.
#' @param columns A character vector of column names to use as part of the generation of unique
#'   combinations of features. Default is \code{c("ecoregionGroup", "speciesCode", "age", "B")}
#'
#' @return
#' Returns a vector of \code{pixelGroup} in the original order of the input \code{pixelDataTable}.
#' This should likely be added to the \code{pixelDataTable} object immediately.
#'
#' @export
#' @importFrom data.table setkey setorderv
#' @importFrom plyr mapvalues
#' @importFrom SpaDES.core paddedFloatToChar
generatePixelGroups <- function(pixelDataTable, maxPixelGroup,
                                columns = c("ecoregionGroup", "speciesCode", "age", "B")) {
  columnsOrig <- columns
  columns <- columns[columns %in% names(pixelDataTable)]
  columns2 <- paste0(columns, "2")
  if (!all(columns == columnsOrig))
    message("Creating pixelGroup values, but not using all columns requested. Only using, ",
            paste(columns, collapse = ", "), " instead of ", paste(columnsOrig, collapse = ", "))

  pcd <- pixelDataTable # no copy -- just for simpler name

  if (getOption("LandR.assertions")) {
    pcdOrig <- data.table::copy(pcd)
  }
  # concatenate within rows -- e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
  pcd[, uniqueComboByRow := do.call(paste, as.list(.SD)), .SDcols = columns]

  # concatenate within pixelIndex
  pcd[ , c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"), by = "pixelIndex"]
  pcd[ , c("pixelGroup") := as.integer(maxPixelGroup) + as.integer(factor(uniqueComboByPixelIndex))]

  if (getOption("LandR.assertions")) { # old algorithm
    # prepare object 1 (pcd) for checking below
    pcd[, ord := 1:.N]
    setorderv(pcd, c("pixelIndex"))
    pcd[, pixelGroup2 := mapvalues(pixelGroup, from = unique(pixelGroup), to = as.character(seq_along(unique(pixelGroup))))]
    setorderv(pcd, "ord")

    pcdOld <- data.table::copy(pcdOrig)

    # Convert to unique numeric
    pcdOld[ , c(columns2) := lapply(.SD, function(x) {
      a <- as.integer(factor(x))
    }), .SDcols = columns]

    # concatenate within rows -- e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
    pcdOld[, uniqueComboByRow := as.integer(factor( do.call(paste, as.list(.SD)))),
           .SDcols = columns2]

    # concatenate within pixelIndex
    pcdOld[ , c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"),
            by = "pixelIndex"]
    pcdOld[ , c("pixelGroup") := as.integer(maxPixelGroup) +
              as.integer(factor(uniqueComboByPixelIndex))]
    # prepare object 2 (pcdOld) for checking below
    pcdOld[, ord := 1:.N]
    setorderv(pcdOld, c("pixelIndex"))
    pcdOld[, pixelGroup2 := mapvalues(pixelGroup, from = unique(pixelGroup),
                                      to = as.character(seq_along(unique(pixelGroup))))]
    setorderv(pcdOld, "ord")

    # The check
    if (!identical(pcdOld$pixelGroup2, pcd$pixelGroup2))
      stop("new generatePixelGroups algorithm failing")
  }

  return(pcd$pixelGroup)
}

#' Pull out the values from \code{speciesEcoregion} table for current time
#'
#' @template speciesEcoregion
#' @template currentTime
#'
#' @return
#' The \code{speciesEcoregion} input object, but with data from only one year, the year
#' that is less than or equal to the \code{currentTime}
#'
#' @export
speciesEcoregionLatestYear <- function(speciesEcoregion, currentTime) {
  spEco <- speciesEcoregion[year <= currentTime]
  spEco[year == max(year)]
}

#' @keywords internal
.ageRndUpSuccessionTimestep <- function(age, successionTimestep) {
  as.integer(ceiling(as.numeric(age) / successionTimestep) * successionTimestep)
}

#' The columns in a \code{cohortData} that define "unique"
#'
#' If two pixels have identical values in all of these columns, they are the same \code{pixelGroup}.
#'
#' @export
#' @rdname uniqueDefinitions
uniqueCohortDefinition <- c("pixelGroup", "speciesCode", "age", "B")

#' @export
#' @rdname uniqueDefinitions
uniqueSpeciesEcoregionDefinition <- c("speciesCode", "ecoregionGroup")

#' Summary for \code{cohortData}
#'
#' @template cohortData
#'
#' @export
describeCohortData <- function(cohortData) {
  vals <- c("B", "totalBiomass", "age", "cover")
  names(vals) <- vals
  out <- lapply(vals, function(val)
    .cohortMessages(cohortData, val)
  )
  message(magenta("Pixels with non-NA cover:, ", cohortData[!is.na(cover), length(unique(pixelIndex))]))
}

#' @keywords internal
.cohortMessages <- function(cohortData, val) {
  out <- list()
  if (val %in% colnames(cohortData)) {
    pixelsNA <- NROW(cohortData[is.na(get(val)), unique("pixelIndex"), with = FALSE])
    message(magenta("Pixels with missing", val, ":", format(pixelsNA, big.mark = ",")))
    pixelsZero <- NROW(cohortData[, all(get(val) == 0), by = "pixelIndex"][get("V1") == TRUE])
    message(magenta("Pixels with all(",val," == 0): ", format(pixelsZero, big.mark = ",")))
    pixelsBiomassNonZero <- NROW(cohortData[, any(get(val) > 0), by = "pixelIndex"][get("V1") == TRUE])
    message(magenta("Pixels with all(",val," > 0): ", format(pixelsBiomassNonZero, big.mark = ",")))
    out <- list(pixelsNA = pixelsNA, pixelsZero = pixelsZero, pixelsBiomassNonZero = pixelsBiomassNonZero)
  }
  return(invisible(out))
}

#' Convert Land Cover Classes (LCC) to another value in its neighbourhood
#'
#' This will search around the pixels on \code{rstLCC} that have
#' \code{classesToReplace}, and search in iteratively increasing
#' radii outwards for other Land Cover Classes than the those indicated in
#' \code{classesToReplace}. This will constrain
#' It will then take the cohorts that were in pixels with \code{classesToReplace}
#' and assign them new values in the output object. This function will
#' also check that it must be an \code{ecoregionCode} that already exists in
#' \code{cohortData}, i.e., not create new \code{ecoregionCode} values. See Details.
#'
#' @details
#' This function is designed to be used in highly constrained situations, where it is not
#' just replacing a Land Cover Class by a neighbouring Land Cover Class. But it can
#' be used for the simpler cases of simply replacing a Land Cover Class.
#'
#' @param pixelClassesToReplace Deprecated. Use \code{classesToReplace}
#'
#' @param classesToReplace Integer vector of classes that are are to be replaced,
#'     e.g., 34, 35, 36 on LCC2005, which are burned young, burned 10 year, and cities.
#'
#' @param rstLCC LCC raster, e.g., LCC2005
#'
#' @param theUnwantedPixels An optional vector of pixel IDs that need to be changed.
#'   If not provided, then pixels to change will be taken from the match between
#'   \code{availableERC_by_Sp} and \code{classesToReplace}. Supplying this allows
#'   the user to only replace some of the pixels with a given class.
#'
#' @param ecoregionGroupVec Deprecated. Use \code{availableERC_by_Sp}
#'
#' @param speciesEcoregion Deprecated. Use \code{availableERC_by_Sp}
#'
#' @param availableERC_by_Sp A \code{data.table} or \code{data.frame} with 3 columns:
#'   \code{speciesCode}, \code{initialEcoregionCode} and \code{pixelIndex}.
#'   \code{pixelIndex} is the pixel id for each line in the \code{data.table};
#'   \code{speciesCode} is the species name in the pixel (can have more than one
#'     species per pixel, so multiple rows per pixel); and,
#'   \code{initialEcoregionCode} is the unique codes that are "available" to be
#'     used as a replacement for \code{classesToReplace}. \code{initialEcoregionCode}
#'     must be a character vector, with one or no "_" used as a separator, with the last
#'     component being the Land Cover Class that matches \code{classesToReplace}, e.g.,
#'     \code{"242_18"}. If there is no "_" in this code, then the codes must match the
#'     \code{classesToReplace} exactly, e.g., \code{"11"}.
#'     If \code{pixelIndex} is missing, the function will fill it
#'     with \code{seq(ncell(rstLCC))}. If \code{speciesCode} is missing, the function
#'     will replace it with a dummy value (\code{"allSpecies"}).
#'
#' @template doAssertion
#'
#' @return
#' A \code{data.table} with two columns, \code{pixelIndex} and \code{ecoregionGroup}.
#' This represents the new codes to used in the \code{pixelIndex} locations.
#' These should have no values overlapping with \code{classesToReplace}.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom data.table as.data.table is.data.table rbindlist setnames
#' @importFrom raster raster
#' @importFrom SpaDES.core paddedFloatToChar
#' @importFrom SpaDES.tools spread2
convertUnwantedLCC <- function(classesToReplace = 34:36, rstLCC,
                               availableERC_by_Sp, theUnwantedPixels,
                               ecoregionGroupVec, speciesEcoregion, pixelClassesToReplace,
                               doAssertion = getOption("LandR.assertions", TRUE)) {
  if (!missing(pixelClassesToReplace))
    stop("pixelClassesToReplace is deprecated. Please use classesToReplace")
  if (!missing(ecoregionGroupVec))
    stop("ecoregionGroupVec is deprecated. Please use availableERC_by_Sp")
  if (!missing(speciesEcoregion))
    stop("speciesEcoregion is deprecated. Please use availableERC_by_Sp")

  if (!is.data.table(availableERC_by_Sp))
    if (is.data.frame(availableERC_by_Sp)) {
      availableERC_by_Sp <- as.data.table(availableERC_by_Sp)
    } else {
      stop("availableERC_by_Sp must be a data.table or data.frame")
    }

  if (all(grepl("_", availableERC_by_Sp$initialEcoregionCode))) {
    hasPreDash <- TRUE
    availableERC_by_Sp[, ecoregion := gsub("_.*", "", initialEcoregionCode)]
  } else {
    hasPreDash <- FALSE
  }

  genericSpeciesName <- "allSpecies"
  if (!"speciesCode" %in% names(availableERC_by_Sp)) {
    availableERC_by_Sp[, speciesCode := genericSpeciesName]
  }
  if (!"pixelIndex" %in% names(availableERC_by_Sp)) {
    availableERC_by_Sp[, pixelIndex := seq(ncell(rstLCC))]
  }
  if (!is.character(availableERC_by_Sp$initialEcoregionCode) && hasPreDash) {
    availableERC_by_Sp[, initialEcoregionCode := as.character(initialEcoregionCode)]
  }

  # if (!all(grepl("_", availableERC_by_Sp$initialEcoregionCode))) {
  #   numChar <- max(nchar(availableERC_by_Sp$initialEcoregionCode), na.rm = TRUE)
  #   availableERC_by_Sp[!is.na(initialEcoregionCode),
  #                      initialEcoregionCode := paste0("1_", paddedFloatToChar(as.integer(initialEcoregionCode), numChar))]
  # }

  if (missing(theUnwantedPixels)) {
    if (FALSE) { # This is old section... can be deleted soon (April 10, 2019, Eliot)
      rstUnwantedLCC <- integer(length(ecoregionGroupVec))
      rstUnwantedLCC[] <- NA
      rstUnwantedLCC[gsub(".*_", "", ecoregionGroupVec) %in% classesToReplace] <- 1
      theUnwantedPixels1 <- which(rstUnwantedLCC == 1)
      theUnwantedPixels1 <- theUnwantedPixels1[theUnwantedPixels1 %in%
                                                 unique(availableERC_by_Sp$pixelIndex)]
    } else {
      theUnwantedRows <- gsub(".*_", "", availableERC_by_Sp$initialEcoregionCode) %in%
        as.character(classesToReplace)
      theUnwantedPixels <- sort(unique(availableERC_by_Sp[theUnwantedRows, "pixelIndex"])[[1]])
    }
  }

  if (doAssertion) {
    #  stop("values of 34 and 35 on pixelCohortData and sim$LCC2005 don't match")
  }
  iterations <- 1
  # remove the lines that have the code "classesToReplace"
  availableERG2 <- if (hasPreDash) {
    availableERC_by_Sp[-which(gsub(".*_", "", initialEcoregionCode) %in%
                                classesToReplace)]
  } else {
    availableERC_by_Sp[-which(initialEcoregionCode %in% classesToReplace)]
  }

  availableERG2 <- unique(availableERG2, by = c("speciesCode", "initialEcoregionCode"))
  availableERG2[, `:=`(pixelIndex = NULL)]

  numCharIEC <- max(nchar(availableERC_by_Sp$initialEcoregionCode), na.rm = TRUE)
  if (hasPreDash) {
    numCharEcoregion <- max(nchar(availableERG2$ecoregion), na.rm = TRUE)
    numCharLCCCodes <- numCharIEC - numCharEcoregion - 1 # minus 1 for the dash
    availableERG2[, `:=`(ecoregion = NULL)]
  } else {
    numCharLCCCodes <- numCharIEC
  }

  currentLenUnwantedPixels <- length(theUnwantedPixels)
  repeatsOnSameUnwanted <- 0

  while (length(theUnwantedPixels) > 0) {
    message("Converting unwanted LCCs: ", length(theUnwantedPixels), " pixels remaining.")
    out <- spread2(rstLCC, start = theUnwantedPixels, asRaster = FALSE,
                   iterations = iterations, allowOverlap = TRUE, spreadProb = 1)
    out <- out[initialPixels != pixels] # rm pixels which are same as initialPixels --> these are known wrong
    iterations <- iterations + 1
    out[, lcc := rstLCC[][pixels]]
    out[lcc %in% c(classesToReplace), lcc := NA]
    out <- na.omit(out)
    out5 <- availableERC_by_Sp[out[, state := NULL], allow.cartesian = TRUE,
                               on = c("pixelIndex" = "initialPixels"), nomatch = NA] # join the availableERC_by_Sp which has initialEcoregionCode

    if (hasPreDash) {
      out5[, possERC := paste0(ecoregion, "_",
                               paddedFloatToChar(as.integer(lcc), padL = numCharLCCCodes, padR = 0))]
    } else {
      out5[, possERC := lcc]
    }
    out7 <- out5[availableERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode"), nomatch = NA]
    out6 <- na.omit(out7)

    # sanity check -- don't let an infinite loop
    if (currentLenUnwantedPixels == length(theUnwantedPixels)) {
      repeatsOnSameUnwanted <- repeatsOnSameUnwanted + 1
    } else {
      currentLenUnwantedPixels <- length(theUnwantedPixels)
      repeatsOnSameUnwanted <- 0
    }

    if (repeatsOnSameUnwanted > 5) {
      out2 <- data.table(newPossLCC = NA, pixelIndex = theUnwantedPixels, ecoregionGroup = NA)
      message("  removing ", NROW(theUnwantedPixels), " pixel of class ",
              paste(rstLCC[theUnwantedPixels], collapse = ", "), " because couldn't",
              " find a suitable replacement")
      theUnwantedPixels <- integer()
    }

    if (NROW(out6) > 0) {
      ## take random sample of available, weighted by abundance
      rowsToKeep <- out6[, list(keep = .resample(.I, 1)), by = c("pixelIndex")]
      out8 <- out6[rowsToKeep$keep]
      out2 <- out8[, list(newPossLCC = lcc, pixelIndex)]
      if (hasPreDash) {
        out2[, initialEcoregion := substr(out8[, initialEcoregionCode], 1, numCharEcoregion)]
        out2[, ecoregionGroup := paste0(initialEcoregion, "_",
                                        paddedFloatToChar(as.integer(newPossLCC), padL = 2, padR = 0))] #nolint
        out2[, initialEcoregion := NULL]
      } else {
        out2[, ecoregionGroup := as.integer(newPossLCC)] #nolint
      }

      ## remove combinations of ecoregionGroup and speciesCode that don't exist
      ## -- Now this excludes B = 0
      keepPixels <- unique(out2$pixelIndex)
      theUnwantedPixels <- theUnwantedPixels[!theUnwantedPixels %in% keepPixels]
      out2 <- unique(out2)

      if (!exists("out3")) {
        out3 <- out2
      } else {
        out3 <- rbindlist(list(out2, out3))
      }
    }
  }

  if (!exists("out3")) {
    out3 <- data.table(pixelIndex = NA, ecoregionGroup = NA)[!is.na(pixelIndex)]
  } else {
    # setnames(out3, c("initialPixels", "initialEcoregionCode"), c("pixelIndex", "ecoregionGroup"))
    out3[, `:=`(newPossLCC = NULL)]
    # out3 <- unique(out3, by = c("pixelIndex", "ecoregionGroup"))
    out3 <- unique(out3)
  }
  out3
}

#' Generate template \code{cohortData} table
#'
#' Internal function used by \code{\link{makeAndCleanInitialCohortData}}.
#'
#' @param inputDataTable A \code{data.table} with columns described above.
#'
#' @param pixelGroupBiomassClass Round B to the nearest \code{pixelGroupBiomassClass}
#'   to establish unique \code{pixelGroups}.
#'
#' @template doAssertion
#'
#' @param rescale Logical. If \code{TRUE}, the default, cover for each species will be rescaled
#'   so all cover in \code{pixelGroup} or pixel sums to 100.
#'
#' @importFrom crayon blue
#' @importFrom data.table melt setnames
#' @keywords internal
.createCohortData <- function(inputDataTable, pixelGroupBiomassClass,
                              doAssertion = getOption("LandR.assertions", TRUE), rescale = TRUE) {
  coverColNames <- grep(colnames(inputDataTable), pattern = "cover", value = TRUE)
  newCoverColNames <- gsub("cover\\.", "", coverColNames)
  setnames(inputDataTable, old = coverColNames, new = newCoverColNames)
  message(blue("Create initial cohortData object, with no pixelGroups yet"))
  cohortData <- data.table::melt(inputDataTable,
                                 value.name = "cover",
                                 measure.vars = newCoverColNames,
                                 variable.name = "speciesCode")
  cohortData[, coverOrig := cover]
  if (isTRUE(doAssertion))
    if (any(duplicated(cohortData)))
      warning(".createCohortData: cohortData contains duplicate rows.")

  if (doAssertion)
    #describeCohortData(cohortData)
    message(blue("assign B = 0 and age = 0 for pixels where cover = 0,\n",
                 "because cover is most reliable dataset"))

  hasCover0 <- which(cohortData[["cover"]] == 0)
  cncd <- colnames(cohortData)
  if (any(c("age", "logAge") %in% cncd)) {
    set(cohortData, hasCover0, "age", 0L)
    set(cohortData, hasCover0, "logAge", -Inf)
  }
  if (any(c("B", "totalBiomass") %in% cncd)) {
    set(cohortData, hasCover0, "B", 0)
    message(blue("assign totalBiomass = 0 if sum(cover) = 0 in a pixel, ",
                 "\n  because cover is most reliable dataset"))
    cohortData <- cohortData[, sum(cover) == 0, by = "pixelIndex"][V1 == TRUE][
      cohortData, on = "pixelIndex"][V1 == TRUE, totalBiomass := 0L]
    cohortData[, V1 := NULL]
  }

  ## CRAZY TODO: DIVIDE THE COVER BY 2 for DECIDUOUS -- will only affect mixed stands
  # message(crayon::blue(paste("POSSIBLE ALERT:",
  #                            "assume deciduous cover is 1/2 the conversion to B as conifer")))
  # cohortData[speciesCode == "Popu_sp", cover := asInteger(cover / 2)]
  set(cohortData, NULL, "cover", as.numeric(cohortData[["cover"]]))
  if (isTRUE(rescale))
    cohortData[ , cover := {
      sumCover <- sum(cover)
      if (sumCover > 100) {
        cover <- cover / (sumCover + 0.0001) * 100L
      }
      cover
    }, by = "pixelIndex"]

  set(cohortData, NULL, "cover", asInteger(cohortData[["cover"]]))

  if (any(c("B", "totalBiomass") %in% cncd)) {
    # Biomass -- by cohort (NOTE: divide by 100 because cover is percent)
    set(cohortData, NULL, "B", as.numeric(cohortData[["B"]]))
    message(blue("Divide total B of each pixel by the relative cover of the cohorts"))
    cohortData[ , B := mean(totalBiomass) * cover / 100, by = "pixelIndex"]
    message(blue("Round B to nearest P(sim)$pixelGroupBiomassClass"))
    cohortData[ , B := ceiling(B / pixelGroupBiomassClass) * pixelGroupBiomassClass]

    message(blue("Set B to 0 where cover > 0 and age = 0, because B is least quality dataset"))
    cohortData[cover > 0 & age == 0, B := 0L]
    cohortData[ , totalBiomass := asInteger(totalBiomass)]
    set(cohortData, NULL, "B", asInteger(cohortData[["B"]]))
  }

  return(cohortData)
}

#' Generate initial \code{cohortData} table
#'
#' Takes a single \code{data.table} input, which has the following columns in addition to
#'    others that will be labelled with species name, and contain percent cover of each:
#' \itemize{
#'   \item \code{pixelIndex} (integer)
#'   \item \code{age} (integer)
#'   \item \code{logAge} (numeric)
#'   \item \code{initialEcoregionCode} (factor)
#'   \item \code{totalBiomass} (integer)
#'   \item \code{lcc} (integer)
#'   \item \code{rasterToMatch} (integer)
#'   \item \code{speciesCode} (factor)
#'   \item \code{cover} (integer)
#'   \item \code{coverOrig} (integer)
#'   \item \code{B} (integer)
#' }
#'
#' @param inputDataTable A \code{data.table} with columns described above.
#'
#' @param sppColumns A vector of the names of the columns in \code{inputDataTable} that
#'   represent percent cover by species, rescaled to sum up to 100\%.
#'
#' @param pixelGroupBiomassClass Round B to the nearest \code{pixelGroupBiomassClass}
#'   to establish unique \code{pixelGroups}.
#'
#' @template doAssertion
#'
#' @param doSubset Turns on/off subsetting. Defaults to \code{TRUE}.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom crayon blue
#' @importFrom data.table melt setnames
#' @importFrom reproducible Cache
#' @rdname makeAndCleanInitialCohortData
makeAndCleanInitialCohortData <- function(inputDataTable, sppColumns, pixelGroupBiomassClass,
                                          doAssertion = getOption("LandR.assertions", TRUE),
                                          doSubset = TRUE) {
  ### Create groupings
  if (doAssertion) {
    expectedColNames <- c("age", "logAge", "initialEcoregionCode", "totalBiomass",
                          "lcc", "pixelIndex")
    if (!all(expectedColNames %in% colnames(inputDataTable)))
      stop("Column names for inputDataTable must include ",
           paste(expectedColNames, collapse = " "))
    if (!all(sppColumns %in% colnames(inputDataTable)))
      stop("Species names are incorrect")
    if (!all(unlist(lapply(inputDataTable[, sppColumns, with = FALSE],
                           function(x) all(x >= 0 & x <= 100)))))
      stop("Species columns are not percent cover between 0 and 100. This may",
           " be because they more NA values than the Land Cover raster")
  }

  cohortData <- Cache(.createCohortData, inputDataTable = inputDataTable,
                      pixelGroupBiomassClass = pixelGroupBiomassClass,
                      doAssertion = doAssertion)

  ######################################################
  # Impute missing ages on poor age dataset
  ######################################################
  cohortDataMissingAge <- cohortData[, hasBadAge := all(age == 0 & cover > 0) |
                                       any(is.na(age)), by = "pixelIndex"][hasBadAge == TRUE]

  if (NROW(cohortDataMissingAge) > 0) {
    cohortDataMissingAgeUnique <- unique(cohortDataMissingAge,
                                         by = c("initialEcoregionCode", "speciesCode"))[
                                           , .(initialEcoregionCode, speciesCode)]
    cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[
      cohortData, on = c("initialEcoregionCode", "speciesCode"), nomatch = 0]
    ageModel <- quote(lme4::lmer(age ~ B * speciesCode + (1 | initialEcoregionCode) + cover))
    cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[, .(B, age, speciesCode,
                                                                 initialEcoregionCode, cover)]
    cohortDataMissingAgeUnique <- subsetDT(cohortDataMissingAgeUnique,
                                           by = c("initialEcoregionCode", "speciesCode"),
                                           doSubset = doSubset)
    message(blue("Impute missing age values: started", Sys.time()))
    outAge <- Cache(statsModel, modelFn = ageModel,
                    uniqueEcoregionGroups = unique(cohortDataMissingAgeUnique$initialEcoregionCode),
                    .specialData = cohortDataMissingAgeUnique,
                    omitArgs = ".specialData")
    message(blue("                           completed", Sys.time()))
    message(outAge$rsq)

    ## allow.new.levels = TRUE because some groups will have only NA for age for all species
    cohortDataMissingAge[
      , imputedAge := pmax(0L, asInteger(predict(outAge$mod,
                                                 newdata = cohortDataMissingAge,
                                                 allow.new.levels = TRUE)))]
    cohortData <- cohortDataMissingAge[, .(pixelIndex, imputedAge, speciesCode)][
      cohortData, on = c("pixelIndex", "speciesCode")]
    cohortData[!is.na(imputedAge), `:=`(age = imputedAge, logAge = log(imputedAge))]
    cohortData[, `:=`(imputedAge = NULL)]
  }
  cohortData[, `:=`(hasBadAge = NULL)]

  #######################################################
  # set B to zero if age is zero because B is lowest quality dataset
  #######################################################
  message(blue("Set recalculate totalBiomass as sum(B);",
               "many biomasses will have been set to 0 in previous steps"))
  cohortData[cover > 0 & age == 0, B := 0L]
  cohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]

  # This was unused, but for beta regression, this can allow 0s and 1s without
  #   needing a separate model for zeros and ones
  # https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
  # cohortData[ , coverProp := (cover/100 * (NROW(cohortData) - 1) + 0.5) / NROW(cohortData)]

  return(cohortData)
}

#' Subset a \code{data.table} with random subsampling within \code{by} groups
#'
#' @param DT A \code{data.table}
#' @param by Character vector of column names to use for groups
#' @param doSubset Logical or numeric indicating the number of subsamples to use
#'
#' @export
subsetDT <- function(DT, by, doSubset = TRUE) {
  if (!is.null(doSubset)) {
    if (!isFALSE(doSubset)) {
      sam <- if (is.numeric(doSubset)) doSubset else 50
      message("subsampling initial dataset for faster estimation of maxBiomass parameter: ",
              "using maximum of ", sam, " samples per combination of ecoregionGroup and speciesCode. ",
              "Change 'doSubset' to a different number if this is not enough")
      # subset -- add line numbers of those that were sampled
      a <- DT[, list(lineNum = .I[sample(.N, size = min(.N, sam))]), by = by]
      # Select only those row numbers from whole dataset
      DT <- DT[a$lineNum]
    }
  }
  return(DT)
}

#' The generic statistical model to run (\code{lmer} or \code{glmer})
#'
#' This does a few things including R squared, gets the fitted values.
#' It appears that running the models "as is" without this wrapper does not work with \code{Cache}.
#' The return of the model in a list solves this problem.
#' For Caching, the \code{.specialData} should be "omitted" via \code{omitArgs}, and
#' \code{uniqueEcoregionGroups} should not be omitted.
#'
#' @param modelFn A quoted expression of type \code{package::model(Y ~ X, ...)}, omitting
#'   the \code{data} argument. E.g. \code{lme4::glmer(Y ~ X + (X|G), family = poisson)}
#' @param uniqueEcoregionGroups Unique values of \code{ecoregionGroups}.
#'   This is the basis for the statistics, and can be used to optimize caching,
#'   e.g. ignore \code{.specialData} in \code{.omitArgs}.
#' @param sumResponse a sum of all the response variable values
#'   Also to be used to optimize caching, e.g. ignore \code{.specialData}
#'   in \code{.omitArgs}.
#' @param .specialData The custom dataset required for the model.
#'
#' @export
#' @importFrom crayon blue magenta red
#' @importFrom lme4 glmer lmer
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom stats as.formula glm fitted predict
statsModel <- function(modelFn, uniqueEcoregionGroups, sumResponse, .specialData) {
  ## convert model call to vector of arguments
  modelArgs <- as.character(modelFn)
  names(modelArgs) <- names(modelFn)

  ## function name
  fun <- modelArgs[[1]]

  ## get formula and check
  form <- tryCatch(as.formula(modelArgs[2]),
                   error = function(e)
                     stop(paste("Could not convert '", modelArgs[2], "'to formula.",
                                "Check if formula is of type 'Y ~ X'")))

  ## check the no of grouping levels
  if (grepl("\\|", modelArgs[2])) {
    groupVar <- sub("\\).*", "", sub(".*\\| ", "",  modelArgs[2]))
    keepGrouping <- NROW(unique(.specialData[, ..groupVar])) >= 2

    if (!keepGrouping) {
      form <- sub("\\+ \\(.*\\|.*\\)", "", modelArgs[2])
      form <- as.formula(form)

      ## change to glm if dropping random effects
      fun <- "stats::glm"

      whereFam <- grep("family", names(modelArgs))
      modelFn2 <- if (length(whereFam))
        call(fun, form, family = modelArgs[[whereFam]]) else
          call(fun, form)

      message(blue("Grouping variable "), red("only has one level. "),
              blue("Formula changed to\n",
                   magenta(paste0(format(modelFn2, appendLF = FALSE), collapse = ""))))
    }
  }

  ## get function and check
  fun <- reproducible:::.extractFunction(fun)
  if (!is.function(fun))
    stop(paste0("Can't find the function '", modelArgs[1], "'.",
                " Is the function name correct and the package installed?"))

  ## prepare arguments, and strip function from list
  modelArgs <- as.list(modelArgs)
  modelArgs[[2]] <- form
  names(modelArgs)[2] <- "formula"
  modelArgs[[1]] <- NULL
  modelArgs$data <- quote(.specialData)

  mod <- do.call(fun, modelArgs)

  list(mod = mod, pred = fitted(mod), rsq = MuMIn::r.squaredGLMM(mod))
}

#' Default columns that define pixel groups
#'
#' @export
columnsForPixelGroups <- c("ecoregionGroup", "speciesCode", "age", "B")

#' Generate \code{cohortData} table per pixel:
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template doAssertion
#'
#' @return
#' An expanded \code{cohortData} \code{data.table} with a new \code{pixelIndex} column.
#'
#' @export
#' @importFrom raster getValues ncell
addPixels2CohortData <- function(cohortData, pixelGroupMap,
                                 doAssertion = getOption("LandR.assertions", TRUE)) {
  assertCohortData(cohortData, pixelGroupMap, doAssertion = doAssertion)

  pixelGroupTable <- na.omit(data.table(pixelGroup = getValues(pixelGroupMap),
                                        pixelIndex = 1:ncell(pixelGroupMap)))
  pixelCohortData <- cohortData[pixelGroupTable, on = "pixelGroup",
                                nomatch = 0, allow.cartesian = TRUE]

  assertPixelCohortData(pixelCohortData, pixelGroupMap, doAssertion = doAssertion)

  return(pixelCohortData)
}

#' Add number of pixels per \code{pixelGroup} and add it has a new column to \code{cohortData}
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template doAssertion
#'
#' @return
#' An \code{cohortData} \code{dat.table} with a new \code{noPixels}
#' column
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom raster maxValue
addNoPixel2CohortData <- function(cohortData, pixelGroupMap,
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  assertCohortData(cohortData, pixelGroupMap, doAssertion = doAssertion)

  noPixelsXGroup <- data.table(noPixels = tabulate(pixelGroupMap[]),
                               pixelGroup = c(1:maxValue(pixelGroupMap)))

  pixelCohortData <- cohortData[noPixelsXGroup, on = "pixelGroup", nomatch = 0]

  if (doAssertion) {
    test1 <- length(setdiff(pixelCohortData$pixelGroup, cohortData$pixelGroup)) > 0
    test2 <- sum(unique(pixelCohortData[, .(pixelGroup, noPixels)])$noPixels) !=
      sum(!is.na(pixelGroupMap[]) & pixelGroupMap[] != 0)  ## 0's have no cohorts.

    if (test1 | test2)
      stop("pixelGroups differ between pixelCohortData/pixelGroupMap and cohortData")
  }

  return(pixelCohortData)
}


#' Make the \code{cohortData} table, while modifying the temporary
#' \code{pixelCohortData} that will be used to prepare other files.
#'
#' Takes a \code{pixelCohortData} table (see \code{makeAndCleanInitialCohortData}),
#'   the \code{speciesEcoregion} list and returns a modified \code{pixelCohortData} and
#'   the \code{cohortData} tables to be used in the simulation.
#'   This function mainly removes unnecessary columns from \code{pixelCohortData},
#'   subsets pixels with \code{biomass > 0}, generates \code{pixelGroups},
#'   and adds \code{ecoregionGroup} and \code{totalBiomass} columns to \code{pixelCohortData}.
#'   \code{cohortData} is then created by subsetting unique combinations of \code{pixelGroup} and
#'   whatever columns are listed in \code{columnsForPixelGroups}.
#'   The resulting \code{cohortData} table has the following columns:
#' \itemize{
#'   \item \code{speciesCode} (factor)
#'   \item \code{ecoregionGroup} (factor)
#'   \item \code{pixelGroup} (integer)
#'   \item \code{age} (integer)
#'   \item \code{B} (integer)
#' }
#'
#' @template pixelCohortData
#' @param columnsForPixelGroups Default columns that define pixel groups
#' @template speciesEcoregion
#'
#' @return
#' A list with a modified \code{pixelCohortData} and \code{cohortData} \code{data.table}s.
#'
#' @export
#' @importFrom data.table melt setnames set
#' @importFrom reproducible Cache
makeCohortDataFiles <- function(pixelCohortData, columnsForPixelGroups, speciesEcoregion) {
  ## make ecoregioGroup a factor (again) and remove unnecessary cols.
  # refactor because the "_34" and "_35" ones are still levels
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))]
  cols <- intersect(c("logAge", "coverOrig", "totalBiomass",
                      "initialEcoregionCode", "cover", "lcc"),
                    names(pixelCohortData))
  set(pixelCohortData, j = cols, value = NULL)

  ## select pixels with biomass and generate pixel groups
  pixelCohortData <- pixelCohortData[B > 0]
  cd <- pixelCohortData[, .SD, .SDcols = c("pixelIndex", columnsForPixelGroups)]
  pixelCohortData[, pixelGroup := Cache(generatePixelGroups, cd, maxPixelGroup = 0,
                                        columns = columnsForPixelGroups)]

  ########################################################################
  ## rebuild ecoregion, ecoregionMap objects -- some initial ecoregions disappeared (e.g., 34, 35, 36)
  ## rebuild biomassMap object -- biomasses have been adjusted
  ecoregionsWeHaveParametersFor <- levels(speciesEcoregion$ecoregionGroup)

  pixelCohortData <- pixelCohortData[ecoregionGroup %in% ecoregionsWeHaveParametersFor] # keep only ones we have params for
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))]
  pixelCohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]

  cohortData <- unique(pixelCohortData, by = c("pixelGroup", columnsForPixelGroups))
  cohortData[, `:=`(pixelIndex = NULL)]

  assertUniqueCohortData(cohortData, c("pixelGroup", "ecoregionGroup", "speciesCode"))
  return(list(cohortData = cohortData, pixelCohortData = pixelCohortData))
}

#' Create new cohorts based on provenance table with unique \code{pixelGroup} and add to \code{cohortData}
#'
#' @param newPixelCohortData the cohorts that were harvested
#' @template cohortData
#' @template pixelGroupMap
#' @template currentTime
#' @param successionTimestep succession timestep used in the simulation
#'
#' @return A \code{data.table} with a new \code{cohortData}
#'
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
#' @export
plantNewCohorts <- function(newPixelCohortData, cohortData, pixelGroupMap,
                            currentTime, successionTimestep) {

  browser()
   ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum B per ecoregion, join to new cohort data
  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, pixelGroup := getValues(pixelGroupMap)[pixelIndex]]
    } else {
      stop("newPixelCohortData must have either pixelIndex or pixelGroup")
    }
  }

  if (!is.null(newPixelCohortData[["pixelIndex"]]))
    set(newPixelCohortData, NULL, "pixelIndex", NULL)

  #remove duplicate pixel groups that were harvested
  duplicates <- duplicated(newPixelCohortData[, .(speciesCode, ecoregionGroup, pixelGroup, age)])
  newCohortData <- newPixelCohortData[!duplicates]

  #Plant trees
  newCohortData[B == 0, age := 2]
  #Give the planted trees 2 * maxANPP - newly regenerating cohorts receive 1x maxANPP. 2x seems overkill
  newCohortData[B == 0, B := asInteger(2 * maxANPP)]

  #Here we subset cohortData instead of setting added columns to NULL. However, as these are 'new' cohorts, this is okay
  newCohortData <- newCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, age, B,
                                               mortality = 0L, aNPPAct = 0L)]

  if (getOption("LandR.assertions")) {
    if (isTRUE(NROW(unique(newCohortData, by = uniqueCohortDefinition)) != NROW(newCohortData))){
      browser()
      stop("Duplicated new cohorts in a pixelGroup. Please debug LandR:::.plantNewCohorts")
    #in this situation, it may be caused by not replanting all species. Not sure if this will come up.
    }
  }

  cohortData <- rbindlist(list(cohortData, newCohortData), fill = TRUE, use.names = TRUE)

  return(cohortData)
}

#' Add cohorts to \code{cohortData} and \code{pixelGroupMap}
#'
#' This is a wrapper for  \code{generatePixelGroups}, \code{initiateNewCohort} and updates to
#' \code{pixelGroupMap} via assignment to new \code{pixelIndex} values in \code{newPixelCohortData}.
#' By running these all together, there is less chance that they will diverge.
#' There are some checks internally for consistency.
#'
#' Does the following:
#' \enumerate{
#'   \item add new cohort data into \code{cohortData};
#'   \item assign initial \code{B} and \code{age} for new cohort;
#'   \item assign the new \code{pixelGroup} to the pixels that have new cohort;
#'   \item update the \code{pixelGroup} map.
#' }
#'
#' @template newPixelCohortData
#' @template cohortData
#' @template pixelGroupMap
#' @template currentTime
#' @template speciesEcoregion
#'
#' @param treedHarvestPixelTableSinceLastDisp A data.table with at least 2 columns, \code{pixelIndex} and \code{pixelGroup}.
#'   This will be used in conjunction with \code{cohortData} and \code{pixelGroupMap}
#'   to ensure that everything matches correctly.
#' @param successionTimestep The time between successive seed dispersal events.
#'   In LANDIS-II, this is called "Succession Timestep". This is used here
#' @param provenanceTable A \code{data.table} with three columns:
#' New cohorts are initiated at the \code{ecoregionGroup} \code{speciesEcoregion} from the
#' corresponding \code{speciesEcoregion} listed in the \code{Provenance} column
#'
#' @param verbose Integer, where increasing number is increasing verbosity. Currently,
#'    only level 1 exists; but this may change.
#'
#'
#' @template doAssertion
#'
#' @return
#' A list of length 2, \code{cohortData} and \code{pixelGroupMap}, with
#' \code{newPixelCohortData} inserted.
#'
#' @export
#' @rdname updateCohortDataPostHarvest
#' @importFrom crayon green magenta
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom SpaDES.core paddedFloatToChar
#' @importFrom stats na.omit
updateCohortDataPostHarvest <- function(newPixelCohortData, cohortData, pixelGroupMap, currentTime,
                             speciesEcoregion, treedHarvestPixelTableSinceLastDisp = NULL,
                             successionTimestep, provenanceTable,
                             verbose = getOption("LandR.verbose", TRUE),
                             doAssertion = getOption("LandR.assertions", TRUE)) {

  maxPixelGroup <- as.integer(maxValue(pixelGroupMap))

  if (!is.null(treedHarvestPixelTableSinceLastDisp)) {
    pixelGroupMap[treedHarvestPixelTableSinceLastDisp$pixelIndex] <- 0L
  }
  relevantPixels <- pixelGroupMap[][newPixelCohortData$pixelIndex]
  zeroOnPixelGroupMap <- relevantPixels == 0

  newPixelCohortData[B == 0, age := 1] #this is necessary for now.
  #Can't lose the biomass of non-zero (non harvested) trees

  if (all(zeroOnPixelGroupMap)) {
    #non-zeroes indicate partial harvest - this is not yet operational

    setkey(newPixelCohortData, ecoregionGroup, speciesCode)
    setkey(provenanceTable, ecoregionGroup, speciesCode)
    browser()
    newPixelCohortData <- newPixelCohortData[provenanceTable]

    #ecoregionGroup should remain the same, Provenance will be a newly added column

    #NA in provenance means this area would not be planted with this species

    #Remove all redundant age 0 species. These were created by multiple same-species cohorts
    duplicates <- duplicated(newPixelCohortData[, .(speciesCode, ecoregionGroup, pixelIndex, age)])
    newPixelCohortData <- newPixelCohortData[!duplicates]

    #this block must be run here, in case some pixelGroups change due to no existing maxB
    specieseco_current <- speciesEcoregionLatestYear(speciesEcoregion, currentTime)
    specieseco_current <- setkey(specieseco_current[, .(speciesCode, maxANPP, maxB, ecoregionGroup)],
                                 speciesCode, ecoregionGroup)

    specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]

    setkey(newPixelCohortData, speciesCode, ecoregionGroup)
    browser() #check that newPixelCohortData looks correct
    newPixelCohortData <- specieseco_current[newPixelCohortData]

    #Perhaps this needs to be parameterized, if Provenance is not a desired column?
    columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B", 'maxB', 'maxANPP', 'Provenance')

    cd <- newPixelCohortData[,c("pixelIndex", columnsForPG), with = FALSE]

    newPixelCohortData[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = maxPixelGroup,
                                                           columns = columnsForPG)]

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(newPixelCohortData[, c("pixelIndex", "pixelGroup")],
                             by = c("pixelIndex"))
  } else {
    # This is for situations where there are some empty pixels being filled,
    #   and some occupied pixels getting infilling. This requires a wholesale
    #   re-pixelGroup
    message("This harvest situation is untested and probably incorrect")
    pixelIndex <- which(pixelGroupMap[] %in% cohortData$pixelGroup)
    cohortDataPixelIndex <- data.table(pixelIndex = pixelIndex,
                                       pixelGroup = pixelGroupMap[][pixelIndex])
    cdLong <- cohortDataPixelIndex[cohortData, on = "pixelGroup", allow.cartesian = TRUE]
    cohorts <- rbindlist(list(cdLong, newPixelCohortData), use.names = TRUE, fill = TRUE)

    columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B")
    cd <- cohorts[, c("pixelIndex", columnsForPG), with = FALSE]

    cohorts[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L, columns = columnsForPG)]

    # Bring to pixelGroup level -- this will squash the data.table
    if (is.null(cohorts[["sumB"]])) {
      cohorts[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]
    }
    allCohortData <- cohorts[ , .(ecoregionGroup = ecoregionGroup[1],
                                  mortality = mortality[1],
                                  aNPPAct = aNPPAct[1],
                                  sumB = sumB[1]),
                              by = uniqueCohortDefinition]

    theNewOnes <- is.na(allCohortData$B)
    cohortData <- allCohortData[!theNewOnes]
    newPixelCohortData <- allCohortData[theNewOnes]

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(cohorts[, c("pixelIndex", "pixelGroup")], by = c("pixelIndex"))
  }

  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (doAssertion) {
    if (!isTRUE(all(pixelsToChange$pixelGroup == pixelGroupMap[][pixelsToChange$pixelIndex])))
      stop("pixelGroupMap and newPixelCohortData$pixelGroupMap don't match in updateCohortData fn")
  }

  ## give B in pixels that have serotiny/resprouting
  # newPixelCohortData[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]

  ##########################################################
  # Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone)
  ##########################################################

  cohortData <- plantNewCohorts(newPixelCohortData, cohortData,
                                pixelGroupMap, currentTime = currentTime,
                                successionTimestep = successionTimestep)

  outs <- rmMissingCohorts(cohortData, pixelGroupMap)

  assertCohortData(outs$cohortData, outs$pixelGroupMap,
                   doAssertion = doAssertion, verbose = verbose)

  if (doAssertion) {
    maxPixelGroupFromCohortData <- max(outs$cohortData$pixelGroup)
    maxPixelGroup <- as.integer(maxValue(outs$pixelGroupMap))
    test1 <- (!identical(maxPixelGroup, maxPixelGroupFromCohortData))
    if (test1) {
      stop("The sim$pixelGroupMap and cohortData have unmatching pixelGroup.",
           " They must be matching.",
           " If this occurs, please contact the module developers")
    }
  }

  if (verbose > 0) {
    nPixForest <- sum(!is.na(outs$pixelGroupMap[]))
    nPixGrps <- length(unique(outs$cohortData$pixelGroup))
    nPixNoPixGrp <- sum(outs$pixelGroupMap[] == 0, na.rm = TRUE)
    nPixTreed <- sum(outs$pixelGroupMap[] != 0, na.rm = TRUE)

    nDigits <- max(nchar(c(nPixForest, nPixGrps, nPixNoPixGrp))) + 3
    message(crayon::magenta("NUMBER OF FORESTED PIXELS          :",
                            paddedFloatToChar(nPixForest, padL = nDigits, pad = " ")))
    message(crayon::magenta("NUMBER OF PIXELS WITH TREES        :",
                            paddedFloatToChar(nPixTreed, padL = nDigits, pad = " ")))
    message(crayon::magenta("NUMBER OF UNIQUE PIXELGROUPS       :",
                            paddedFloatToChar(nPixGrps, padL = nDigits, pad = " ")))
    message(crayon::magenta("NUMBER OF PIXELS WITH NO PIXELGROUP:",
                            paddedFloatToChar(nPixNoPixGrp, padL = nDigits, pad = " ")))
  }

  return(list(cohortData = outs$cohortData,
              pixelGroupMap = outs$pixelGroupMap))
}
