utils::globalVariables(c(
  ".", "..cols", "..colsToSubset", ".I", ":=", "..groupVar", "age", "aNPPAct", "cover", "coverOrig",
  "ecoregion", "ecoregionGroup", "hasBadAge",
  "imputedAge", "initialEcoregion", "initialEcoregionCode", "initialPixels",
  "lcc", "maxANPP", "maxB", "maxB_eco", "mortality",
  "newPossLCC", "noPixels", "oldSumB", "ord", "outBiomass", "oldEcoregionGroup",
  "pixelGroup2", "pixelIndex", "pixels", "planted", "Provenance", "possERC",
  "speciesposition", "speciesGroup", "speciesInt", "state", "sumB",
  "temppixelGroup", "toDelete", "totalBiomass", "totalCover",
  "uniqueCombo", "uniqueComboByRow", "uniqueComboByPixelIndex", "V1", "year"
))

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
#' @template cohortDefinitionCols
#'
#' @param treedFirePixelTableSinceLastDisp A data.table with at least 2 columns, \code{pixelIndex}
#'   and \code{pixelGroup}.
#'   This will be used in conjunction with \code{cohortData} and \code{pixelGroupMap}
#'   to ensure that everything matches correctly.
#' @param successionTimestep The time between successive seed dispersal events.
#'   In LANDIS-II, this is called "Succession Timestep". This is used here
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
#' @importFrom reproducible paddedFloatToChar
#' @importFrom stats na.omit
updateCohortData <- function(newPixelCohortData, cohortData, pixelGroupMap, currentTime,
                             speciesEcoregion, treedFirePixelTableSinceLastDisp = NULL,
                             successionTimestep,
                             cohortDefinitionCols = c("pixelGroup", "age", "speciesCode"),
                             verbose = getOption("LandR.verbose", TRUE),
                             doAssertion = getOption("LandR.assertions", TRUE)) {
  maxPixelGroup <- as.integer(maxValue(pixelGroupMap))

  if (!is.null(treedFirePixelTableSinceLastDisp)) {
    ## only set to 0 the pixels that became empty (have no survivors)
    emptyPix <- setdiff(
      treedFirePixelTableSinceLastDisp$pixelIndex,
      newPixelCohortData[type == "survivor", pixelIndex]
    )
    pixelGroupMap[emptyPix] <- 0L
  }
  relevantPixels <- pixelGroupMap[][newPixelCohortData$pixelIndex]
  zeroOnPixelGroupMap <- relevantPixels == 0

  if (!"age" %in% colnames(newPixelCohortData)) {
    newPixelCohortData[, age := 1L]
  }

  if (all(zeroOnPixelGroupMap)) {
    # Deal with pixels on the map that have no pixelGroup -- these are burned
    # pixels --> the entirely newly regenerated pixels do not require a
    # re-pixelGroupMaping  -- can just add to existing pixelGroup values
    if (verbose > 0) {
      message(crayon::green("  Regenerating only burnt pixels with no survivors (i.e. resprouting & serotiny)"))
    }
    columnsForPG <- c("ecoregionGroup", "speciesCode", "age") # no Biomass because they all have zero
    cd <- newPixelCohortData[, c("pixelIndex", columnsForPG), with = FALSE]
    newPixelCohortData[, pixelGroup := generatePixelGroups(cd,
                                                           maxPixelGroup = maxPixelGroup,
                                                           columns = columnsForPG
    )] # ,
    # successionTimestep = successionTimestep)

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(newPixelCohortData[, c("pixelIndex", "pixelGroup")],
                             by = c("pixelIndex")
    )
  } else {
    # This is for situations where there are some empty pixels being filled,
    #   and some occupied pixels getting infilling. This requires a wholesale
    #   re-pixelGroup
    if (verbose > 0) {
      message(crayon::green("  Regenerating open and pixels with B (likely after seed dispersal, or partial mortality following disturbance)"))
    }

    pixelIndex <- which(pixelGroupMap[] %in% cohortData$pixelGroup)

    # remove unnecessary columns before making cohortDataLong
    if ("prevMortality" %in% names(cohortData)) {
      cohortData[, prevMortality := NULL]
    }

    if ("year" %in% names(newPixelCohortData)) {
      newPixelCohortData[, year := NULL]
    }

    cohortDataPixelIndex <- data.table(
      pixelIndex = pixelIndex,
      pixelGroup = pixelGroupMap[][pixelIndex]
    )
    cdLong <- cohortDataPixelIndex[cohortData, on = "pixelGroup", allow.cartesian = TRUE]
    cohorts <- rbindlist(list(cdLong, newPixelCohortData), use.names = TRUE, fill = TRUE)

    columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B")
    cd <- cohorts[, c("pixelIndex", columnsForPG), with = FALSE]
    cohorts[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L, columns = columnsForPG)]

    # Bring to pixelGroup level -- this will squash the data.table
    if (is.null(cohorts[["sumB"]])) {
      cohorts[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]
    }
    # Old way that does not preserve additional 'unknown' columns in cohortData
    # allCohortData <- cohorts[ , .(ecoregionGroup = ecoregionGroup[1],
    #                               mortality = mortality[1],
    #                               aNPPAct = aNPPAct[1],
    #                               sumB = sumB[1]),
    #                           by = uniqueCohortDefinition]

    # newer way that will potentially conflict with LandR.CS due to differing aNPP
    colsToSubset <- setdiff(colnames(cohortData), c("pixelIndex"))
    allCohortData <- cohorts[!duplicated(cohorts[, .(pixelGroup, speciesCode, ecoregionGroup, age)]), ..colsToSubset]

    theNewOnes <- is.na(allCohortData$B)
    cohortData <- allCohortData[!theNewOnes]
    newPixelCohortData <- allCohortData[theNewOnes]

    # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(cohorts[, c("pixelIndex", "pixelGroup")], by = c("pixelIndex"))
  }

  # update pixelGroupMap
  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (doAssertion) {
    if (!isTRUE(all(pixelsToChange$pixelGroup == pixelGroupMap[][pixelsToChange$pixelIndex]))) {
      stop("pixelGroupMap and newPixelCohortData$pixelGroupMap don't match in updateCohortData fn")
    }
  }

  ## give B in pixels that have serotiny/resprouting
  # newPixelCohortData[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]

  ##########################################################
  # Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone)
  ##########################################################
  cohortData <- .initiateNewCohorts(newPixelCohortData, cohortData,
                                    pixelGroupMap,
                                    currentTime = currentTime,
                                    speciesEcoregion = speciesEcoregion,
                                    successionTimestep = successionTimestep
  )

  outs <- rmMissingCohorts(cohortData, pixelGroupMap, cohortDefinitionCols = cohortDefinitionCols)

  if (!is.null(outs$cohortData$sumB)) {
    outs$cohortData[, sumB := NULL]
  }

  assertCohortData(outs$cohortData, outs$pixelGroupMap,
                   cohortDefinitionCols = cohortDefinitionCols,
                   doAssertion = doAssertion, verbose = verbose
  )

  if (doAssertion) {
    maxPixelGroupFromCohortData <- max(outs$cohortData$pixelGroup)
    maxPixelGroup <- as.integer(maxValue(outs$pixelGroupMap))
    test1 <- (!identical(maxPixelGroup, maxPixelGroupFromCohortData))
    if (test1) {
      stop(
        "The sim$pixelGroupMap and cohortData have unmatching pixelGroup.",
        " They must be matching.",
        " If this occurs, please contact the module developers"
      )
    }
  }

  if (verbose > 0) {
    nPixForest <- sum(!is.na(outs$pixelGroupMap[]))
    nPixGrps <- length(unique(outs$cohortData$pixelGroup))
    nPixNoPixGrp <- sum(outs$pixelGroupMap[] == 0, na.rm = TRUE)
    nPixTreed <- sum(outs$pixelGroupMap[] != 0, na.rm = TRUE)

    nDigits <- max(nchar(c(nPixForest, nPixGrps, nPixNoPixGrp))) + 3
    message(crayon::magenta(
      "NUMBER OF FORESTED PIXELS          :",
      paddedFloatToChar(nPixForest, padL = nDigits, pad = " ")
    ))
    message(crayon::magenta(
      "NUMBER OF PIXELS WITH TREES        :",
      paddedFloatToChar(nPixTreed, padL = nDigits, pad = " ")
    ))
    message(crayon::magenta(
      "NUMBER OF UNIQUE PIXELGROUPS       :",
      paddedFloatToChar(nPixGrps, padL = nDigits, pad = " ")
    ))
    message(crayon::magenta(
      "NUMBER OF PIXELS WITH NO PIXELGROUP:",
      paddedFloatToChar(nPixNoPixGrp, padL = nDigits, pad = " ")
    ))
  }

  return(list(
    cohortData = outs$cohortData,
    pixelGroupMap = outs$pixelGroupMap
  ))
}

#' Initiate new cohorts
#'
#' Calculate new values for \code{B}, add \code{age}, then \code{rbindlist} this
#' with \code{cohortData}.
#'
#' @template newPixelCohortData
#' @template cohortData
#' @template cohortDefinitionCols
#'
#' @return A \code{data.table} with a new \code{rbindlist}ed \code{cohortData}
#'
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
#' @rdname updateCohortData
.initiateNewCohorts <- function(newPixelCohortData, cohortData, pixelGroupMap, currentTime,
                                cohortDefinitionCols = c("pixelGroup", "speciesCode", "age"),
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
  if (!is.null(newPixelCohortData[["pixelIndex"]])) {
    set(newPixelCohortData, NULL, "pixelIndex", NULL)
  }
  newPixelCohortData <- newPixelCohortData[!duplicated(newPixelCohortData), ] ## faster than unique

  specieseco_current <- speciesEcoregionLatestYear(speciesEcoregion, currentTime)
  specieseco_current <- setkey(
    specieseco_current[, .(speciesCode, maxANPP, maxB, ecoregionGroup)],
    speciesCode, ecoregionGroup
  )

  ## Note that after the following join, some cohorts will be lost due to lack of
  ##  parameters in speciesEcoregion. These need to be modified in pixelGroupMap.
  # missingNewPixelCohortData <- newPixelCohortData[!specieseco_current, on = uniqueSpeciesEcoregionDefinition]
  specieseco_current <- specieseco_current[!is.na(maxB)]
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]
  newPixelCohortData <- specieseco_current[newPixelCohortData, on = uniqueSpeciesEcoregionDefinition]
  newPixelCohortData <- newPixelCohortData[!is.na(maxB)]

  if (any(newPixelCohortData$age > 1)) {
    stop("newPixelCohortData should only have new cohorts aged 1")
  }
  set(newPixelCohortData, NULL, "age", 1L) ## set age to 1

  ## Ceres: this was causing new cohorts to be initialized with maxANPP.
  ## instead, calculate total biomass of older cohorts
  # set(newPixelCohortData, NULL, "sumB", 0L)
  if (!is.null(newPixelCohortData[["sumB"]])) {
    set(newPixelCohortData, NULL, "sumB", NULL)
  }
  cohortData[age >= successionTimestep, oldSumB := sum(B, na.rm = TRUE), by = "pixelGroup"]

  newPixelCohortData <- unique(cohortData[, .(pixelGroup, oldSumB)],
                               by = "pixelGroup"
  )[newPixelCohortData, on = "pixelGroup"]
  set(newPixelCohortData, which(is.na(newPixelCohortData$oldSumB)), "oldSumB", 0) ## faster than [:=]
  setnames(newPixelCohortData, "oldSumB", "sumB")
  set(cohortData, NULL, "oldSumB", NULL)

  if ("B" %in% names(newPixelCohortData)) {
    newPixelCohortData[, B := NULL]
  }
  set(
    newPixelCohortData, NULL, "B",
    asInteger(pmax(1, newPixelCohortData$maxANPP *
                     exp(-1.6 * newPixelCohortData$sumB / newPixelCohortData$maxB_eco)))
  )
  set(newPixelCohortData, NULL, "B", asInteger(pmin(newPixelCohortData$maxANPP, newPixelCohortData$B)))

  newPixelCohortData <- newPixelCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, age, B,
                                               mortality = 0L, aNPPAct = 0L
  )]

  if (getOption("LandR.assertions")) {
    if (isTRUE(NROW(unique(newPixelCohortData, by = cohortDefinitionCols)) != NROW(newPixelCohortData))) {
      stop("Duplicated new cohorts in a pixelGroup. Please debug LandR:::.initiateNewCohorts")
    }
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
#' @template cohortDefinitionCols
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
                             cohortDefinitionCols = c("pixelGroup", "age", "speciesCode"),
                             doAssertion = getOption("LandR.assertions", TRUE)) {
  pgmValues <- data.table(
    pixelGroup = getValues(pixelGroupMap),
    pixelIndex = seq(ncell(pixelGroupMap))
  )

  pgmVals <- na.omit(pgmValues)
  pgmVals <- pgmVals[pixelGroup > 0]
  whPgsStillInCDGoneFromPGM <- !cohortData$pixelGroup %in% pgmVals$pixelGroup
  pgsStillInCDGoneFromPGM <- cohortData[whPgsStillInCDGoneFromPGM, ] # setdiff(unique(cohortData$pixelGroup), unique(pgmVals$pixelGroup))
  # pgsStillInPGMGoneFromCD <- setdiff(unique(pgmVals$pixelGroup), unique(cohortData$pixelGroup))
  whPgsStillInPGMGoneFromCD <- !pgmVals$pixelGroup %in% cohortData$pixelGroup
  pgsStillInPGMGoneFromCD <- pgmVals[whPgsStillInPGMGoneFromCD, ]

  # REMOVE lines in cohortData that are no longer in the pixelGroupMap
  cohortData <- cohortData[!pixelGroup %in% pgsStillInCDGoneFromPGM$pixelGroup]
  # REMOVE pixels in pixelGroupMap that are no longer in the cohortData
  pixelGroupMap[pgsStillInPGMGoneFromCD$pixelIndex] <- NA

  assertCohortData(cohortData, pixelGroupMap,
                   message = "rmMissingCohorts",
                   cohortDefinitionCols = cohortDefinitionCols,
                   doAssertion = doAssertion
  )

  if (NROW(unique(cohortData[pixelGroup == 67724]$ecoregionGroup)) > 1) stop()

  return(list(
    cohortData = cohortData,
    pixelGroupMap = pixelGroupMap
  ))
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
generatePixelGroups <- function(pixelDataTable, maxPixelGroup,
                                columns = c("ecoregionGroup", "speciesCode", "age", "B")) {
  columnsOrig <- columns
  columns <- columns[columns %in% names(pixelDataTable)]
  columns2 <- paste0(columns, "2")
  if (!all(columns == columnsOrig)) {
    message(
      "Creating pixelGroup values, but not using all columns requested. Only using, ",
      paste(columns, collapse = ", "), " instead of ", paste(columnsOrig, collapse = ", ")
    )
  }

  pcd <- pixelDataTable # no copy -- just for simpler name

  if (getOption("LandR.assertions")) {
    pcdOrig <- data.table::copy(pcd)
  }
  # concatenate within rows -- e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
  pcd[, uniqueComboByRow := do.call(paste, as.list(.SD)), .SDcols = columns]

  # concatenate within pixelIndex
  pcd[, c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"), by = "pixelIndex"]
  pcd[, c("pixelGroup") := as.integer(maxPixelGroup) + as.integer(factor(uniqueComboByPixelIndex))]

  if (getOption("LandR.assertions")) { # old algorithm
    # prepare object 1 (pcd) for checking below
    pcd[, ord := 1:.N]
    setorderv(pcd, c("pixelIndex"))
    uniqPG <- unique(pcd$pixelGroup)
    pcd[, pixelGroup2 := mapvalues2(pixelGroup, from = uniqPG, to = as.character(seq_along(uniqPG)))]
    # pcd[, pixelGroup2 := mapvalues(pixelGroup, from = unique(pixelGroup), to = as.character(seq_along(unique(pixelGroup))))]
    setorderv(pcd, "ord")

    pcdOld <- data.table::copy(pcdOrig)

    # Convert to unique numeric
    pcdOld[, c(columns2) := lapply(.SD, function(x) {
      a <- as.integer(factor(x))
    }), .SDcols = columns]

    # concatenate within rows -- e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
    pcdOld[, uniqueComboByRow := as.integer(factor(do.call(paste, as.list(.SD)))),
           .SDcols = columns2
    ]

    # concatenate within pixelIndex
    pcdOld[, c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"),
           by = "pixelIndex"
    ]
    pcdOld[, c("pixelGroup") := as.integer(maxPixelGroup) +
             as.integer(factor(uniqueComboByPixelIndex))]
    # prepare object 2 (pcdOld) for checking below
    pcdOld[, ord := 1:.N]
    setorderv(pcdOld, c("pixelIndex"))

    uniqPG <- unique(pcdOld$pixelGroup)
    pcdOld[, pixelGroup2 := mapvalues2(pixelGroup, from = uniqPG, to = as.character(seq_along(uniqPG)))]

    setorderv(pcdOld, "ord")

    # The check
    if (!identical(pcdOld$pixelGroup2, pcd$pixelGroup2)) {
      stop("new generatePixelGroups algorithm failing")
    }
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
  out <- lapply(vals, function(val) {
    .cohortMessages(cohortData, val)
  })
  message(magenta("Pixels with non-NA cover:, ", cohortData[!is.na(cover), length(unique(pixelIndex))]))
}

#' @keywords internal
.cohortMessages <- function(cohortData, val) {
  out <- list()
  if (val %in% colnames(cohortData)) {
    pixelsNA <- NROW(cohortData[is.na(get(val)), unique("pixelIndex"), with = FALSE])
    message(magenta("Pixels with missing", val, ":", format(pixelsNA, big.mark = ",")))
    pixelsZero <- NROW(cohortData[, all(get(val) == 0), by = "pixelIndex"][get("V1") == TRUE])
    message(magenta("Pixels with all(", val, " == 0): ", format(pixelsZero, big.mark = ",")))
    pixelsBiomassNonZero <- NROW(cohortData[, any(get(val) > 0), by = "pixelIndex"][get("V1") == TRUE])
    message(magenta("Pixels with all(", val, " > 0): ", format(pixelsBiomassNonZero, big.mark = ",")))
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
#' @importFrom reproducible paddedFloatToChar
#' @importFrom SpaDES.tools spread2
convertUnwantedLCC <- function(classesToReplace = 34:36, rstLCC,
                               availableERC_by_Sp, theUnwantedPixels,
                               ecoregionGroupVec, speciesEcoregion, pixelClassesToReplace,
                               doAssertion = getOption("LandR.assertions", TRUE)) {
  if (!missing(pixelClassesToReplace)) {
    stop("pixelClassesToReplace is deprecated. Please use classesToReplace")
  }
  if (!missing(ecoregionGroupVec)) {
    stop("ecoregionGroupVec is deprecated. Please use availableERC_by_Sp")
  }
  if (!missing(speciesEcoregion)) {
    stop("speciesEcoregion is deprecated. Please use availableERC_by_Sp")
  }

  if (!is.data.table(availableERC_by_Sp)) {
    if (is.data.frame(availableERC_by_Sp)) {
      availableERC_by_Sp <- as.data.table(availableERC_by_Sp)
    } else {
      stop("availableERC_by_Sp must be a data.table or data.frame")
    }
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
    theUnwantedRows <- gsub(".*_", "", availableERC_by_Sp$initialEcoregionCode) %in%
      as.character(classesToReplace)
    theUnwantedPixels <- sort(unique(availableERC_by_Sp[theUnwantedRows, "pixelIndex"])[[1]])
  }

  if (doAssertion) {
    #  stop("values of 34 and 35 on pixelCohortData and sim$LCC2005 don't match")
  }
  iterations <- 1
  # remove the lines that have the code "classesToReplace"
  availableERG2 <- if (hasPreDash) {
    availableERC_by_Sp[-which(gsub(".*_", "", initialEcoregionCode) %in% classesToReplace)]
  } else {
    availableERC_by_Sp[-which(initialEcoregionCode %in% classesToReplace)]
  }

  availableERG2 <- unique(availableERG2, by = c("speciesCode", "initialEcoregionCode"))
  if (doAssertion) {
    if (any(gsub(".*_", "", availableERG2$initialEcoregionCode) %in% classesToReplace)) {
      stop("classesToReplace are still considered 'available' forest classes")
    }
  }

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
    out <- spread2(rstLCC,
                   start = theUnwantedPixels, asRaster = FALSE,
                   iterations = iterations, allowOverlap = TRUE, spreadProb = 1
    )
    out <- out[initialPixels != pixels] # rm pixels which are same as initialPixels --> these are known wrong
    iterations <- iterations + 1
    out[, lcc := rstLCC[][pixels]]
    out[lcc %in% c(classesToReplace), lcc := NA]
    out <- na.omit(out)
    out5 <- availableERC_by_Sp[out[, state := NULL],
                               allow.cartesian = TRUE,
                               on = c("pixelIndex" = "initialPixels"), nomatch = NA
    ] # join the availableERC_by_Sp which has initialEcoregionCode

    if (hasPreDash) {
      out5[, possERC := paste0(
        ecoregion, "_",
        paddedFloatToChar(as.integer(lcc), padL = numCharLCCCodes, padR = 0)
      )]
    } else {
      out5[, possERC := lcc]
    }
    out7 <- out5[availableERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode"), nomatch = NA]
    out6 <- na.omit(out7)

    # These ones are missing at least something in the new possERC
    possERCToRm <- out5[!availableERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode")]
    out6 <- out6[!possERC %in% unique(possERCToRm$possERC)]

    # sanity check -- don't let an infinite loop
    if (currentLenUnwantedPixels == length(theUnwantedPixels)) {
      repeatsOnSameUnwanted <- repeatsOnSameUnwanted + 1
    } else {
      currentLenUnwantedPixels <- length(theUnwantedPixels)
      repeatsOnSameUnwanted <- 0
    }

    if (repeatsOnSameUnwanted > 5) {
      out2 <- data.table(newPossLCC = NA, pixelIndex = theUnwantedPixels, ecoregionGroup = NA)
      message(
        "  removing ", NROW(theUnwantedPixels), " pixel of class ",
        paste(rstLCC[theUnwantedPixels], collapse = ", "), " because couldn't",
        " find a suitable replacement"
      )
      pixelsToNA <- theUnwantedPixels
      theUnwantedPixels <- integer()
    }

    if (NROW(out6) > 0) {
      ## take random sample of available, weighted by abundance
      rowsToKeep <- out6[, list(keep = .resample(.I, 1)), by = c("pixelIndex")]
      out8 <- out6[rowsToKeep$keep]
      out2 <- out8[, list(newPossLCC = lcc, pixelIndex)]
      if (hasPreDash) {
        out2[, initialEcoregion := substr(out8[, initialEcoregionCode], 1, numCharEcoregion)]
        out2[, ecoregionGroup := paste0(
          initialEcoregion, "_",
          paddedFloatToChar(as.integer(newPossLCC), padL = 2, padR = 0)
        )] # nolint
        out2[, initialEcoregion := NULL]
      } else {
        out2[, ecoregionGroup := as.integer(newPossLCC)] # nolint
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

  if (exists("pixelsToNA")) {
    ## make sure these pixels get an NA ecoregion by rm them in case they are present
    if (any(out3$pixelIndex %in% pixelsToNA)) {
      out3 <- out3[!pixelIndex %in% pixelsToNA]
    }
    out3 <- rbind(out3, data.table(pixelIndex = pixelsToNA, ecoregionGroup = NA))
  }

  if (doAssertion) {
    if (any(gsub(".*_", "", out3$ecoregionGroup) %in% classesToReplace)) {
      stop("classesToReplace we're not fully removed")
    }
  }

  out3
}

#' Generate template \code{cohortData} table
#'
#' Internal function used by \code{\link{makeAndCleanInitialCohortData}}.
#'
#' @param inputDataTable A \code{data.table} with columns described above.
#'
#' @template doAssertion
#'
#' @param rescale Logical. If \code{TRUE}, the default, cover for each species will be rescaled
#'   so all cover in \code{pixelGroup} or pixel sums to 100.
#'
#' @importFrom crayon blue
#' @importFrom data.table melt setnames
#' @keywords internal
.createCohortData <- function(inputDataTable, # pixelGroupBiomassClass,
                              doAssertion = getOption("LandR.assertions", TRUE), rescale = TRUE,
                              minCoverThreshold = 5) {
  coverColNames <- grep(colnames(inputDataTable), pattern = "cover", value = TRUE)
  newCoverColNames <- gsub("cover\\.", "", coverColNames)
  setnames(inputDataTable, old = coverColNames, new = newCoverColNames)
  message(blue("Create initial cohortData object, with no pixelGroups yet"))
  message(green("-- Begin reconciling data inconsistencies"))

  inputDataTable[, totalCover := rowSums(.SD), .SDcols = newCoverColNames]
  whEnoughCover <- inputDataTable$totalCover > minCoverThreshold
  message(green(
    "  -- Removing all pixels with totalCover <= minCoverThreshold (affects",
    sum(!whEnoughCover),
    "of", NROW(inputDataTable), "pixels)"
  ))
  message(green("     --> resulting in", sum(whEnoughCover), "pixels)"))
  inputDataTable <- inputDataTable[whEnoughCover]

  whAgeEqZero <- which(inputDataTable$age == 0)
  message(green(
    "  -- Setting TotalBiomass in pixel to 0 where age == 0 (affects", length(whAgeEqZero),
    "of", NROW(inputDataTable), "pixels)"
  ))
  message(green("     --> keeping ", NROW(inputDataTable), "pixels)"))
  inputDataTable[whAgeEqZero, `:=`(totalBiomass = 0)]

  whTotalBEqZero <- which(inputDataTable$totalBiomass == 0)
  message(green(
    "  -- Setting age in pixel to 0 where totalBiomass == 0 (affects", length(whTotalBEqZero),
    "of", NROW(inputDataTable), "pixels)"
  ))
  message(green("     --> keeping ", NROW(inputDataTable), "pixels)"))
  inputDataTable[whTotalBEqZero, `:=`(age = 0)]

  cohortData <- data.table::melt(inputDataTable,
                                 value.name = "cover",
                                 measure.vars = newCoverColNames,
                                 variable.name = "speciesCode"
  )

  # Remove all cover <= minCoverThreshold
  whCoverGTMinCover <- which(cohortData$cover > minCoverThreshold)
  message(green(
    "  -- Removing all cohorts with cover <= minCoverThreshold (affects", NROW(cohortData) - length(whCoverGTMinCover),
    "of", NROW(cohortData), "cohorts"
  ))
  message(green("     --> resulting in", length(whCoverGTMinCover), "cohorts)"))
  cohortData <- cohortData[whCoverGTMinCover]
  message(green("     --> resulting in", length(unique(cohortData$pixelIndex)), "pixels)"))

  cohortData[, coverOrig := cover]
  if (isTRUE(doAssertion)) {
    if (any(duplicated(cohortData))) {
      warning(".createCohortData: cohortData contains duplicate rows.")
    }
  }

  # if (doAssertion)
  # describeCohortData(cohortData)
  # message(green("  -- Assign B = 0 and age = 0 for pixels where cover = 0,\n",
  #             "because cover is most reliable dataset"))

  # hasCover0 <- which(cohortData[["cover"]] == 0)
  cncd <- colnames(cohortData)
  # if (any(c("age", "logAge") %in% cncd)) {
  #   set(cohortData, hasCover0, "age", 0L)
  #   set(cohortData, hasCover0, "logAge", .logFloor(0))
  # }
  if (any(c("B", "totalBiomass") %in% cncd)) {
    # set(cohortData, hasCover0, "B", 0)
    # message(green("  -- Assign totalBiomass = 0 if sum(cover) = 0 in a pixel, ",
    #             "  because cover is most reliable dataset"))
    # cohortData <- cohortData[, sum(cover) == 0, by = "pixelIndex"][V1 == TRUE][
    #  cohortData, on = "pixelIndex"][V1 == TRUE, totalBiomass := 0L]
    # cohortData[, V1 := NULL]
  }

  ## CRAZY TODO: DIVIDE THE COVER BY 2 for DECIDUOUS -- will only affect mixed stands
  # message(crayon::green(paste("POSSIBLE ALERT:",
  #                            "assume deciduous cover is 1/2 the conversion to B as conifer")))
  # cohortData[speciesCode == "Popu_sp", cover := asInteger(cover / 2)]

  # Change temporarily to numeric for following calculation
  set(cohortData, NULL, "cover", as.numeric(cohortData[["cover"]]))
  if (isTRUE(rescale)) {
    cohortData[, totalCover := sum(cover), by = "pixelIndex"]
    cohortData[, cover := cover / totalCover * 100]
  }

  # cohortData[ , cover := {
  #   sumCover <- sum(cover)
  #   if (sumCover > 100) {
  #     cover <- cover / (sumCover + 0.0001) * 100L
  #   }
  #   cover
  # }, by = "pixelIndex"]

  # set(cohortData, NULL, "cover", asInteger(cohortData[["cover"]]))

  if (FALSE) {
    if (any(c("B", "totalBiomass") %in% cncd)) {
      # Biomass -- by cohort (NOTE: divide by 100 because cover is percent)
      # set(cohortData, NULL, "B", as.numeric(cohortData[["B"]]))
      set(cohortData, NULL, "B", cohortData[["totalBiomass"]] * cohortData[["cover"]] / 100)
      message(green("  -- Divide total B of each pixel by the relative cover of the cohorts"))

      # cohortData[ , B := mean(totalBiomass) * cover / 100, by = "pixelIndex"]
      # message(blue("Round B to nearest P(sim)$pixelGroupBiomassClass"))
      # cohortData[ , B := ceiling(B / pixelGroupBiomassClass) * pixelGroupBiomassClass]

      message(green("Set B to 0 where cover > 0 and age = 0, because B is least quality dataset"))
      cohortData[cover > 0 & age == 0, B := 0L]
      cohortData[, totalBiomass := asInteger(totalBiomass)]
      set(cohortData, NULL, "B", asInteger(cohortData[["B"]]))
    }
  }

  # clean up
  set(cohortData, NULL, c("totalCover", "coverOrig"), NULL)
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
#' @param imputeBadAgeModel DESCRIPTION NEEDED
#'
#' @param minCoverThreshold DESCRIPTION NEEDED
#'
#' @template doAssertion
#'
#' @param doSubset Turns on/off subsetting. Defaults to \code{TRUE}.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom crayon blue green
#' @importFrom data.table melt setnames
#' @importFrom reproducible Cache .sortDotsUnderscoreFirst messageDF
#' @importFrom pemisc termsInData
#' @rdname makeAndCleanInitialCohortData
makeAndCleanInitialCohortData <- function(inputDataTable, sppColumns,
                                          # pixelGroupBiomassClass,
                                          # pixelGroupAgeClass = 1,
                                          imputeBadAgeModel = quote(lme4::lmer(age ~ B * speciesCode + cover * speciesCode + (1 | initialEcoregionCode))),
                                          minCoverThreshold,
                                          doAssertion = getOption("LandR.assertions", TRUE),
                                          doSubset = TRUE) {
  ### Create groupings
  if (doAssertion) {
    expectedColNames <- c(
      "age", "logAge", "initialEcoregionCode", "totalBiomass",
      "lcc", "pixelIndex"
    )
    if (!all(expectedColNames %in% colnames(inputDataTable))) {
      stop(
        "Column names for inputDataTable must include ",
        paste(expectedColNames, collapse = " ")
      )
    }
    if (!all(sppColumns %in% colnames(inputDataTable))) {
      stop("Species names are incorrect")
    }
    if (!all(unlist(lapply(
      inputDataTable[, sppColumns, with = FALSE],
      function(x) all(x >= 0 & x <= 100)
    )))) {
      stop(
        "Species columns are not percent cover between 0 and 100. This may",
        " be because they more NA values than the Land Cover raster"
      )
    }
  }

  cohortData <- Cache(.createCohortData,
                      inputDataTable = inputDataTable,
                      # pixelGroupBiomassClass = pixelGroupBiomassClass,
                      minCoverThreshold = minCoverThreshold,
                      doAssertion = doAssertion
  )

  ######################################################
  # Impute missing ages on poor age dataset
  ######################################################
  # Cases:
  #  All species cover = 0 yet totalB > 0

  cohortDataMissingAge <- cohortData[
    , hasBadAge :=
      # (age == 0 & cover > 0)#| # ok because cover can be >0 with biomass = 0
      (age > 0 & cover == 0) |
      is.na(age) #|
    # (B > 0 & age == 0) |
    # (B == 0 & age > 0)
  ][hasBadAge == TRUE] # , by = "pixelIndex"]

  if (NROW(cohortDataMissingAge) > 0) {
    cohortDataMissingAgeUnique <- unique(cohortDataMissingAge,
                                         by = c("initialEcoregionCode", "speciesCode")
    )[
      , .(initialEcoregionCode, speciesCode)
    ]
    cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[
      cohortData,
      on = c("initialEcoregionCode", "speciesCode"), nomatch = 0
    ]
    cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[!is.na(cohortDataMissingAgeUnique$age)]
    cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[, .(
      totalBiomass, age, speciesCode,
      initialEcoregionCode, cover
    )]
    zeros <- sapply(cohortDataMissingAgeUnique, function(x) sum(x == 0))
    if (sum(zeros, na.rm = TRUE)) {
      hasZeros <- zeros[zeros > 0]
      message(
        " ", paste(names(hasZeros), collapse = ", "), " had ",
        paste(hasZeros, collapse = ", "), " zeros, respectively"
      )
      warning(" These are being removed from the dataset. If this is not desired; please fix.")
      # terms <- strsplit(gsub(" ", "", as.character(imputeBadAgeModel)), split = "[[:punct:]]+")[[2]][-1] # remove response
      # terms <- unique(terms)
      # terms <- terms[terms %in% colnames(cohortDataMissingAgeUnique)]
      terms <- termsInData(imputeBadAgeModel, cohortDataMissingAgeUnique)
      lapply(terms, function(x) {
        cohortDataMissingAgeUnique <<- cohortDataMissingAgeUnique[get(x) != 0]
      })
    }
    cohortDataMissingAgeUnique <- subsetDT(cohortDataMissingAgeUnique,
                                           by = c("initialEcoregionCode", "speciesCode"),
                                           doSubset = doSubset
    )
    message(blue("Impute missing age values: started", Sys.time()))

    outAge <- Cache(statsModel,
                    modelFn = imputeBadAgeModel,
                    uniqueEcoregionGroups = .sortDotsUnderscoreFirst(
                      as.character(unique(cohortDataMissingAgeUnique$initialEcoregionCode))
                    ),
                    .specialData = cohortDataMissingAgeUnique,
                    omitArgs = ".specialData"
    )
    message(blue("                           completed", Sys.time()))

    # paste with capture.output keeps table structure intact
    messageDF(outAge$rsq, 3, "blue")

    ## allow.new.levels = TRUE because some groups will have only NA for age for all species
    cohortDataMissingAge[
      , imputedAge := pmax(0L, asInteger(predict(outAge$mod,
                                                 newdata = cohortDataMissingAge,
                                                 allow.new.levels = TRUE
      )))
    ]

    cohortData <- cohortDataMissingAge[, .(pixelIndex, imputedAge, speciesCode)][
      cohortData,
      on = c("pixelIndex", "speciesCode")
    ]
    cohortData[!is.na(imputedAge), `:=`(age = imputedAge, logAge = .logFloor(imputedAge))]
    cohortData[, `:=`(imputedAge = NULL)]
  }
  # # Round ages to nearest pixelGroupAgeClass
  # set(cohortData, NULL, "age", asInteger(cohortData$age / pixelGroupAgeClass) *
  #       as.integer(pixelGroupAgeClass))

  cohortData[, `:=`(hasBadAge = NULL)]

  # #######################################################
  # # set B to zero if age is zero because B is lowest quality dataset
  # #######################################################
  # message(blue("Set recalculate totalBiomass as sum(B);",
  #              "many biomasses will have been set to 0 in previous steps"))
  # cohortData[cover > 0 & age == 0, B := 0L]
  # cohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]

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
#' @param indices Logical. If \code{TRUE}, this will return vector of row indices only. Defaults
#'   to \code{FALSE}, i.e., return the subsampled \code{data.table}
#'
#' @export
#' @examples
#' library(data.table)
#' dt <- data.table(Lett = sample(LETTERS, replace = TRUE, size = 1000), Nums = 1:100)
#' dt1 <- subsetDT(dt, by = "Lett", doSubset = 3)
subsetDT <- function(DT, by, doSubset = TRUE, indices = FALSE) {
  if (!is.null(doSubset)) {
    if (!isFALSE(doSubset)) {
      sam <- if (is.numeric(doSubset)) doSubset else 50
      message(
        "subsampling initial dataset for faster model estimation: ",
        "using maximum of ", sam, " samples per combination of ecoregionGroup and speciesCode. ",
        "Change 'doSubset' to a different number if this is not enough"
      )
      # subset -- add line numbers of those that were sampled
      a <- DT[, list(lineNum = .I[sample(.N, size = min(.N, sam))]), by = by]
      # Select only those row numbers from whole dataset
      if (isFALSE(indices)) {
        DT <- DT[a$lineNum]
      } else {
        DT <- a$lineNum
      }
    } else {
      if (isTRUE(indices)) {
        DT <- 1:nrow(DT)
      }
    }
  } else {
    if (isTRUE(indices)) {
      DT <- 1:nrow(DT)
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
                   error = function(e) {
                     stop(paste(
                       "Could not convert '", modelArgs[2], "'to formula.",
                       "Check if formula is of type 'Y ~ X'"
                     ))
                   }
  )

  ## check the no of grouping levels
  if (grepl("\\|", modelArgs[2])) {
    groupVar <- sub("\\).*", "", sub(".*\\| ", "", modelArgs[2]))
    keepGrouping <- NROW(unique(.specialData[, ..groupVar])) >= 2

    if (!keepGrouping) {
      form <- sub("\\+ \\(.*\\|.*\\)", "", modelArgs[2])
      form <- as.formula(form)

      ## change to glm if dropping random effects
      fun <- "stats::glm"

      whereFam <- grep("family", names(modelArgs))
      modelFn2 <- if (length(whereFam)) {
        call(fun, form, family = modelArgs[[whereFam]])
      } else {
        call(fun, form)
      }

      message(
        blue("Grouping variable "), red("only has one level. "),
        blue(
          "Formula changed to\n",
          magenta(paste0(format(modelFn2, appendLF = FALSE), collapse = ""))
        )
      )
    }
  }

  ## get function and check
  fun <- reproducible:::.extractFunction(fun)
  if (!is.function(fun)) {
    stop(paste0(
      "Can't find the function '", modelArgs[1], "'.",
      " Is the function name correct and the package installed?"
    ))
  }

  ## prepare arguments, and strip function from list
  modelArgs <- as.list(modelArgs)
  modelArgs[[2]] <- form
  names(modelArgs)[2] <- "formula"
  modelArgs[[1]] <- NULL
  modelArgs$data <- quote(.specialData)

  isChar <- unlist(lapply(modelArgs, is.character))
  if (any(isChar)) {
    modelArgs[isChar] <- lapply(modelArgs[isChar], function(yyy) {
      eval(parse(text = yyy))
    })
  }

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
#' @template cohortDefinitionCols
#'
#' @template doAssertion
#'
#' @return
#' An expanded \code{cohortData} \code{data.table} with a new \code{pixelIndex} column.
#'
#' @export
#' @importFrom raster getValues ncell
addPixels2CohortData <- function(cohortData, pixelGroupMap,
                                 cohortDefinitionCols = c("pixelGroup", "age", "speciesCode"),
                                 doAssertion = getOption("LandR.assertions", TRUE)) {
  assertCohortData(cohortData, pixelGroupMap,
                   cohortDefinitionCols = cohortDefinitionCols,
                   doAssertion = doAssertion
  )

  pixelGroupTable <- na.omit(data.table(
    pixelGroup = getValues(pixelGroupMap),
    pixelIndex = 1:ncell(pixelGroupMap)
  ))
  pixelCohortData <- cohortData[pixelGroupTable,
                                on = "pixelGroup",
                                nomatch = 0, allow.cartesian = TRUE
  ]

  assertPixelCohortData(pixelCohortData, pixelGroupMap, doAssertion = doAssertion)

  return(pixelCohortData)
}

#' Add number of pixels per \code{pixelGroup} and add it has a new column to \code{cohortData}
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#' @template cohortDefinitionCols
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
                                  cohortDefinitionCols = c("pixelGroup", "age", "speciesCode"),
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  assertCohortData(cohortData, pixelGroupMap,
                   cohortDefinitionCols = cohortDefinitionCols, doAssertion = doAssertion
  )

  noPixelsXGroup <- data.table(
    noPixels = tabulate(pixelGroupMap[]),
    pixelGroup = c(1:maxValue(pixelGroupMap))
  )

  pixelCohortData <- cohortData[noPixelsXGroup, on = "pixelGroup", nomatch = 0]

  if (doAssertion) {
    test1 <- length(setdiff(pixelCohortData$pixelGroup, cohortData$pixelGroup)) > 0
    test2 <- sum(unique(pixelCohortData[, .(pixelGroup, noPixels)])$noPixels) !=
      sum(!is.na(pixelGroupMap[]) & pixelGroupMap[] != 0) ## 0's have no cohorts.

    if (test1 | test2) {
      stop("pixelGroups differ between pixelCohortData/pixelGroupMap and cohortData")
    }
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
#' @param pixelGroupAgeClass Integer. When assigning \code{pixelGroup} membership, this defines the
#'   resolution of ages that will be considered 'the same \code{pixelGroup}', e.g., if it is 10,
#'   then 6 and 14 will be the same.
#' @param pixelGroupBiomassClass Integer. When assigning \code{pixelGroup} membership, this defines
#'   the resolution of biomass that will be considered 'the same \code{pixelGroup}', e.g., if it is
#'   100, then 5160 and 5240 will be the same
#' @param minAgeForGrouping Minimum age for regrouping. This may be because there is a source of
#'   ages for young stands/trees that is very reliable, such as a fire database.
#'   Ages below this will not be grouped together. Defaults to -1, meaning treat all ages equally.
#'   If this is related to known ages from a high quality database, then use age of the oldest
#'   trees in that database.
#' @param pixelFateDT A \code{data.table} of \code{pixelFateDT}; if none provided, will make an empty one.
#'
#' @template speciesEcoregion
#'
#' @return A list with a modified \code{pixelCohortData}, \code{cohortData}, and \code{pixelFateDT}
#' \code{data.table}s.
#'
#' @export
#' @importFrom data.table melt setnames set
#' @importFrom reproducible Cache
#' @importFrom utils tail
makeCohortDataFiles <- function(pixelCohortData, columnsForPixelGroups, speciesEcoregion,
                                pixelGroupBiomassClass, pixelGroupAgeClass, minAgeForGrouping = 0,
                                pixelFateDT) {
  ## make ecoregioGroup a factor (again) and remove unnecessary cols.
  # refactor because the "_34" and "_35" ones are still levels
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))]
  cols <- intersect(
    c(
      "logAge", "coverOrig", "totalBiomass",
      "initialEcoregionCode", "cover", "lcc"
    ),
    names(pixelCohortData)
  )
  set(pixelCohortData, j = cols, value = NULL)


  # Round ages to nearest pixelGroupAgeClass
  pixelCohortData[
    age > minAgeForGrouping,
    age := asInteger(age / pixelGroupAgeClass) *
      as.integer(pixelGroupAgeClass)
  ]

  # Round Biomass to nearest pixelGroupBiomassClass
  message(blue("Round B to nearest P(sim)$pixelGroupBiomassClass"))
  pixelCohortData[ # age > minAgeForGrouping,
    , B := asInteger(B / pixelGroupBiomassClass) * as.integer(pixelGroupBiomassClass)
  ]

  # Remove B == 0 cohorts after young removals
  message(green("  -- Removing cohorts with B = 0 and age > 0 -- these were likely poor predictions from updateYoungBiomasses"))
  whBEqZeroAgeGT0 <- which(pixelCohortData$B == 0 & pixelCohortData$age > 0)

  if (length(whBEqZeroAgeGT0) > 0) {
    pixelCohortData2 <- pixelCohortData[-whBEqZeroAgeGT0]
  } else {
    pixelCohortData2 <- pixelCohortData
  }

  lostPixels <- setdiff(pixelCohortData$pixelIndex, pixelCohortData2$pixelIndex)
  message(green("     affected", length(whBEqZeroAgeGT0), "cohorts, in", length(lostPixels), "pixels;"))
  lenUniquePix <- length(unique(pixelCohortData2$pixelIndex))
  message(green("     leaving", lenUniquePix, "pixels"))
  pixelFateDT <- pixelFate(pixelFateDT,
                           fate = "rm pixels with Biomass == 0, after updating young cohort B",
                           length(lostPixels), runningPixelTotal = lenUniquePix
  )

  pixelCohortData <- pixelCohortData2
  # # Set B to 0 if age is 0
  # whAgeZero <- which(pixelCohortData$age == 0)
  # if (length(whAgeZero)) {
  #   message(green("    -- There were", length(whAgeZero), "pixels with age = 0; forcing B to zero"))
  #   pixelCohortData[whAgeZero, B := 0L]
  # }

  ## select pixels with biomass and generate pixel groups
  # pixelCohortData <- pixelCohortData[B >= 0]

  ########################################################################
  ## rebuild ecoregion, ecoregionMap objects -- some initial ecoregions disappeared (e.g., 34, 35, 36)
  ## rebuild biomassMap object -- biomasses have been adjusted
  ## There will be some speciesCode * ecoregionGroups combinations for which we only had age = 0 and B = 0
  if (is.factor(speciesEcoregion$ecoregionGroup)) {
    ## we need to check available params for speciesXecoregionGroup combos, this will be achieved with a
    # ecoregionsWeHaveParametersFor <- levels(speciesEcoregion$ecoregionGroup)
  } else {
    stop(
      "speciesEcoregion$ecoregionGroup is supposed to be a factor; please go back in this",
      "module and correct this."
    )
  }

  message(blue("Removing some pixels because their species * ecoregionGroup combination has no age or B data to estimate ecoregion traits:"))
  # message(blue(paste(sort(unique(pixelCohortData[!ecoregionGroup %in% ecoregionsWeHaveParametersFor]$ecoregionGroup)), collapse = ", ")))
  cols <- c("speciesCode", "ecoregionGroup")
  messageDF(
    colour = "blue",
    pixelCohortData[!speciesEcoregion, on = cols][, ..cols][, list(numPixelsRemoved = .N), by = cols], # anti-join
  )

  # pixelCohortData <- pixelCohortData[ecoregionGroup %in% ecoregionsWeHaveParametersFor] # keep only ones we have params for
  pixelCohortData <- speciesEcoregion[, ..cols][pixelCohortData, on = cols, nomatch = 0]

  # Lost some ecoregionGroups -- refactor
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))]

  cd <- pixelCohortData[, .SD, .SDcols = c("pixelIndex", columnsForPixelGroups)]
  pixelCohortData[, pixelGroup := Cache(generatePixelGroups, cd,
                                        maxPixelGroup = 0,
                                        columns = columnsForPixelGroups
  )]

  pixelCohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]

  cohortData <- unique(pixelCohortData, by = c("pixelGroup", columnsForPixelGroups))
  cohortData[, `:=`(pixelIndex = NULL)]

  pixelFateDT <- pixelFate(
    pixelFateDT, "removing ecoregionGroups without enough data to est. maxBiomass",
    tail(pixelFateDT$runningPixelTotal, 1) - NROW(unique(pixelCohortData$pixelIndex))
  )

  assertUniqueCohortData(cohortData, c("pixelGroup", "ecoregionGroup", "speciesCode"))
  return(list(cohortData = cohortData, pixelCohortData = pixelCohortData, pixelFateDT = pixelFateDT))
}

#' Create new cohorts based on provenance table with unique \code{pixelGroup} and add to \code{cohortData}
#'
#' @param newPixelCohortData the cohorts that were harvested
#' @template cohortData
#' @template pixelGroupMap
#' @template currentTime
#' @param successionTimestep succession timestep used in the simulation
#' @param trackPlanting adds column that tracks planted cohorts if \code{TRUE}
#'
#' @return A \code{data.table} with a new \code{cohortData}
#'
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom stats na.omit
#' @export
plantNewCohorts <- function(newPixelCohortData, cohortData, pixelGroupMap,
                            currentTime, successionTimestep, trackPlanting = FALSE) {
  ## get spp "productivity traits" per ecoregion/present year

  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, pixelGroup := getValues(pixelGroupMap)[pixelIndex]]
    } else {
      stop("newPixelCohortData must have either pixelIndex or pixelGroup")
    }
  }

  if (!is.null(newPixelCohortData[["pixelIndex"]])) {
    set(newPixelCohortData, NULL, "pixelIndex", NULL)
  }

  # remove duplicate pixel groups that were harvested
  duplicates <- duplicated(newPixelCohortData[, .(speciesCode, ecoregionGroup, pixelGroup, age)])
  newCohortData <- newPixelCohortData[!duplicates]

  # Plant trees
  newCohortData[, age := 2]
  # Give the planted trees 2 * maxANPP - newly regenerating cohorts receive 1x maxANPP. 2x seems overkill
  newCohortData[, B := asInteger(2 * maxANPP)]

  # Here we subset cohortData instead of setting added columns to NULL. However, as these are 'new' cohorts, this is okay
  newCohortData <- newCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, age, B, Provenance,
                                     mortality = 0L, aNPPAct = 0L
  )]

  if (getOption("LandR.assertions")) {
    if (isTRUE(NROW(unique(newCohortData, by = c("pixelGroup", "age", "speciesCode", "Provenance")))
               != NROW(newCohortData))) {
      stop("Duplicated new cohorts in a pixelGroup. Please debug LandR:::.plantNewCohorts")
      # in this situation, it may be caused by not replanting all species. Not sure if this will come up.
    }
  }

  if (trackPlanting) {
    newCohortData[, planted := TRUE]
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
#' @template cohortDefinitionCols
#' @template currentTime
#' @template speciesEcoregion
#'
#' @param treedHarvestPixelTable A data.table with at least 2 columns, \code{pixelIndex} and \code{pixelGroup}.
#'   This will be used in conjunction with \code{cohortData} and \code{pixelGroupMap}
#'   to ensure that everything matches correctly.
#' @param successionTimestep The time between successive seed dispersal events.
#'   In LANDIS-II, this is called "Succession Timestep". This is used here
#' @param provenanceTable A \code{data.table} with three columns:
#' New cohorts are initiated at the \code{ecoregionGroup} \code{speciesEcoregion} from the
#' corresponding \code{speciesEcoregion} listed in the \code{Provenance} column
#' @param trackPlanting if true, planted cohorts in \code{cohortData} are tracked with \code{TRUE}
#' in column 'planted'
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
#' @importFrom crayon green magenta
#' @importFrom data.table copy rbindlist set setkey
#' @importFrom raster getValues
#' @importFrom reproducible paddedFloatToChar
#' @importFrom stats na.omit
#' @rdname updateCohortDataPostHarvest
updateCohortDataPostHarvest <- function(newPixelCohortData, cohortData, pixelGroupMap, currentTime,
                                        speciesEcoregion, treedHarvestPixelTable = NULL,
                                        successionTimestep, provenanceTable, trackPlanting = FALSE,
                                        cohortDefinitionCols = c("pixelGroup", "age", "speciesCode"),
                                        verbose = getOption("LandR.verbose", TRUE),
                                        doAssertion = getOption("LandR.assertions", TRUE)) {
  cohortData <- copy(cohortData)
  provenanceTable <- copy(provenanceTable)

  maxPixelGroup <- as.integer(maxValue(pixelGroupMap))

  if (!is.null(treedHarvestPixelTable)) {
    pixelGroupMap[treedHarvestPixelTable$pixelIndex] <- 0L
  }
  relevantPixels <- pixelGroupMap[][newPixelCohortData$pixelIndex]
  zeroOnPixelGroupMap <- relevantPixels == 0

  newPixelCohortData[B == 0, age := 1] # this is necessary for now.
  # Can't lose the biomass of non-zero (non harvested) trees

  setkey(newPixelCohortData, ecoregionGroup, speciesCode)
  setkey(provenanceTable, ecoregionGroup, speciesCode)

  newPixelCohortData <- provenanceTable[newPixelCohortData]

  # ecoregionGroup should remain the same, Provenance will be a newly added column

  # NA in provenance means this area would not be planted with this species

  # Remove duplicate species. These were created by cohorts of different age but same species
  duplicates <- duplicated(newPixelCohortData[, .(speciesCode, ecoregionGroup, pixelIndex)])
  newPixelCohortData <- newPixelCohortData[!duplicates]

  # this block must be run here, in case some pixelGroups change due to no existing maxB
  specieseco_current <- speciesEcoregionLatestYear(speciesEcoregion, currentTime)
  specieseco_current <- setkey(
    specieseco_current[, .(speciesCode, maxANPP, maxB, ecoregionGroup)],
    speciesCode, ecoregionGroup
  )

  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]

  setkey(newPixelCohortData, speciesCode, ecoregionGroup)

  newPixelCohortData <- specieseco_current[newPixelCohortData]

  columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "B", "maxB", "maxANPP", "Provenance")

  cd <- newPixelCohortData[, c("pixelIndex", columnsForPG), with = FALSE]

  newPixelCohortData[, pixelGroup := generatePixelGroups(cd,
                                                         maxPixelGroup = maxPixelGroup,
                                                         columns = columnsForPG
  )]

  # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
  pixelsToChange <- unique(newPixelCohortData[, c("pixelIndex", "pixelGroup")],
                           by = c("pixelIndex")
  )

  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (doAssertion) {
    if (!isTRUE(all(pixelsToChange$pixelGroup == pixelGroupMap[][pixelsToChange$pixelIndex]))) {
      stop("pixelGroupMap and newPixelCohortData$pixelGroupMap don't match in updateCohortData fn")
    }
  }


  ##########################################################
  # Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone)
  ##########################################################

  cohortData <- plantNewCohorts(newPixelCohortData, cohortData,
                                pixelGroupMap,
                                currentTime = currentTime,
                                successionTimestep = successionTimestep,
                                trackPlanting = trackPlanting
  )

  outs <- rmMissingCohorts(cohortData, pixelGroupMap, cohortDefinitionCols = cohortDefinitionCols)

  assertCohortData(outs$cohortData, outs$pixelGroupMap,
                   cohortDefinitionCols = cohortDefinitionCols,
                   doAssertion = doAssertion, verbose = verbose
  )

  if (doAssertion) {
    maxPixelGroupFromCohortData <- max(outs$cohortData$pixelGroup)
    maxPixelGroup <- as.integer(maxValue(outs$pixelGroupMap))
    test1 <- (!identical(maxPixelGroup, maxPixelGroupFromCohortData))
    if (test1) {
      stop(
        "The sim$pixelGroupMap and cohortData have unmatching pixelGroup.",
        " They must be matching.",
        " If this occurs, please contact the module developers"
      )
    }
  }

  if (verbose > 0) {
    nPixForest <- sum(!is.na(outs$pixelGroupMap[]))
    nPixGrps <- length(unique(outs$cohortData$pixelGroup))
    nPixNoPixGrp <- sum(outs$pixelGroupMap[] == 0, na.rm = TRUE)
    nPixTreed <- sum(outs$pixelGroupMap[] != 0, na.rm = TRUE)

    nDigits <- max(nchar(c(nPixForest, nPixGrps, nPixNoPixGrp))) + 3
    message(crayon::magenta(
      "NUMBER OF FORESTED PIXELS          :",
      paddedFloatToChar(nPixForest, padL = nDigits, pad = " ")
    ))
    message(crayon::magenta(
      "NUMBER OF PIXELS WITH TREES        :",
      paddedFloatToChar(nPixTreed, padL = nDigits, pad = " ")
    ))
    message(crayon::magenta(
      "NUMBER OF UNIQUE PIXELGROUPS       :",
      paddedFloatToChar(nPixGrps, padL = nDigits, pad = " ")
    ))
    message(crayon::magenta(
      "NUMBER OF PIXELS WITH NO PIXELGROUP:",
      paddedFloatToChar(nPixNoPixGrp, padL = nDigits, pad = " ")
    ))
  }

  return(list(
    cohortData = outs$cohortData,
    pixelGroupMap = outs$pixelGroupMap
  ))
}

#' Create or amend data to a \code{pixelFateDT} object
#'
#' @param pixelFateDT A \code{pixelFateDT} \code{data.table} with 3 columns: \code{fate},
#'   \code{pixelsRemoted}, and \code{runningPixelTotal}.
#' @param fate A character string (length 1) describing in words the change
#' @param pixelsRemoved A numeric indicating how many pixels were removed due to the \code{fate}.
#' @param runningPixelTotal an optional numeric with new, running total. If not supplied,
#'   it will be calculated from the last row of \code{pixelFateDT} \code{runningTotal} minus the
#'   \code{pixelsRemoved}
#'
#' @return A \code{pixelFateDT} object, updated with one extra row.
#'
#' @export
pixelFate <- function(pixelFateDT, fate = NA_character_, pixelsRemoved = 0,
                      runningPixelTotal = NA_integer_) {
  if (missing(pixelFateDT)) {
    pixelFateDT <- data.table(fate = character(), pixelsRemoved = integer(), runningPixelTotal = integer())
  }
  if (is.na(runningPixelTotal)) {
    runningPixelTotal <- tail(pixelFateDT$runningPixelTotal, 1) - pixelsRemoved
  }
  pixelFateDT <- rbindlist(list(pixelFateDT, data.table(
    fate = fate, pixelsRemoved = pixelsRemoved,
    runningPixelTotal = runningPixelTotal
  )))
  pixelFateDT
}

#' Generate and add vegetation type column to \code{cohortData}
#'
#' This function is a simplification of \code{vegTypeMapGenerator}
#' and instead of generating a map, it adds the vegetation type column
#' to the \code{cohortData} table.
#'
#' @param x A \code{cohortData} object
#'
#' @param vegLeadingProportion Numeric between 0-1, determining the relative biomass
#'                             threshold a species needs to pass to be considered "leading".
#'
#' @param mixedType An integer defining whether mixed stands are of any kind of species
#'                  admixture (1), or only when deciduous mixed with conifer (2).
#'                  Defaults to 2.
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param pixelGroupColName Name of the column in \code{pixelGroup} to use.
#'
#' @template doAssertion
#'
#' @param ... Additional arguments.
#'
#' @return \code{x} with a new column, 'leading', coding the vegetation type
#'    of each group defined by \code{pixelGroupColName}
#'
#' @author Eliot McIntire, Ceres Barros, Alex Chubaty
#' @export
#' @importFrom data.table copy data.table setkey setorderv
#' @importFrom utils data
#'
#' @rdname vegTypeGenerator
#' @examples
#' library(data.table)
#' x <- data.table(
#'   pixelGroup = rep(1:2, each = 2), B = c(100, 200, 20, 400),
#'   speciesCode = rep(c("Pice_Gla", "Popu_Tre"), 2)
#' )
#' vegTypeGenerator(x)
vegTypeGenerator <- function(x, vegLeadingProportion = 0.8,
                             mixedType = 2, sppEquiv = NULL, sppEquivCol,
                             pixelGroupColName = "pixelGroup",
                             doAssertion = getOption("LandR.assertions", TRUE), ...) {
  nrowCohortData <- NROW(x)

  leadingBasedOn <- preambleVTG(x, vegLeadingProportion, doAssertion, nrowCohortData)

  if (mixedType == 2) {
    if (is.null(sppEquiv)) {
      sppEquiv <- get(data("sppEquivalencies_CA", package = "LandR", envir = environment()),
                      inherits = FALSE
      )

      # Find the sppEquivCol that best matches what you have in x
      sppEquivCol <- names(sort(sapply(sppEquiv, function(xx) sum(xx %in% unique(x$species))),
                                decreasing = TRUE
      )[1])
      message(paste0(
        "Using mixedType == 2, but no sppEquiv provided. ",
        "Attempting to use data('sppEquivalencies_CA', 'LandR') ",
        "and sppEquivCol == '", sppEquivCol, "'"
      ))
    }
  }

  ## use new vs old algorithm based on size of x. new one (2) is faster in most cases.
  ## enable assertions to view timings for each algorithm before deciding which to use.
  algo <- ifelse(nrowCohortData > 3.5e6, 1, 2)

  pgdAndSc <- c(pixelGroupColName, "speciesCode")
  pgdAndScAndLeading <- c(pgdAndSc, leadingBasedOn)
  totalOfLeadingBasedOn <- paste0("total", leadingBasedOn)
  speciesOfLeadingBasedOn <- paste0("speciesGroup", leadingBasedOn)
  if (algo == 1 || isTRUE(doAssertion)) {
    # slower -- older, but simpler Eliot June 5, 2019
    # 1. Find length of each pixelGroup -- don't include pixelGroups in "by" that have only 1 cohort: N = 1
    cohortData1 <- copy(x)
    systimePre1 <- Sys.time()
    pixelGroupData1 <- cohortData1[, list(N = .N), by = pixelGroupColName]

    # Calculate speciesProportion from cover or B
    pixelGroupData1 <- cohortData1[, ..pgdAndScAndLeading][pixelGroupData1, on = pixelGroupColName]
    set(pixelGroupData1, NULL, totalOfLeadingBasedOn, pixelGroupData1[[leadingBasedOn]])

    if (identical(leadingBasedOn, "cover")) {
      pixelGroupData1[N != 1, (totalOfLeadingBasedOn) := sum(cover, na.rm = TRUE), by = pixelGroupColName]
      pixelGroupData1 <- pixelGroupData1[, list(sum(cover, na.rm = TRUE), totalcover[1]), by = pgdAndSc]
    } else {
      pixelGroupData1[N != 1, (totalOfLeadingBasedOn) := sum(B, na.rm = TRUE), by = pixelGroupColName]
      pixelGroupData1 <- pixelGroupData1[, list(sum(B, na.rm = TRUE), totalB[1]), by = pgdAndSc]
    }
    setnames(pixelGroupData1, old = c("V1", "V2"), new = c(speciesOfLeadingBasedOn, totalOfLeadingBasedOn))

    set(pixelGroupData1, NULL, "speciesProportion", pixelGroupData1[[speciesOfLeadingBasedOn]] /
          pixelGroupData1[[totalOfLeadingBasedOn]])
    systimePost1 <- Sys.time()

    setorderv(pixelGroupData1, pixelGroupColName)
  }

  if (algo == 2 || isTRUE(doAssertion)) {
    # Replacement algorithm to calculate speciesProportion
    #  Logic is similar to above --
    #  1. sort by pixelGroup
    #  2. calculate N, use this to repeat itself (instead of a join above)
    #  3. calculate speciesProportion, noting to calculate with by only if N > 1, otherwise
    #     it is a simpler non-by calculation
    cohortData2 <- copy(x)
    systimePre2 <- Sys.time()
    setkeyv(cohortData2, pgdAndSc)
    # setorderv(x, pixelGroupColName)
    pixelGroupData2 <- cohortData2[, list(N = .N), by = pixelGroupColName]
    cohortData2 <- cohortData2[, ..pgdAndScAndLeading]

    N <- rep.int(pixelGroupData2$N, pixelGroupData2$N)
    wh1 <- N == 1
    set(cohortData2, which(wh1), totalOfLeadingBasedOn, cohortData2[[leadingBasedOn]][wh1])
    if (identical(leadingBasedOn, "cover")) {
      totalBNot1 <- cohortData2[!wh1, list(N = .N, totalcover = sum(cover, na.rm = TRUE)), by = pixelGroupColName]
    } else {
      totalBNot1 <- cohortData2[!wh1, list(N = .N, totalB = sum(B, na.rm = TRUE)), by = pixelGroupColName]
    }
    totalBNot1 <- rep.int(totalBNot1[[totalOfLeadingBasedOn]], totalBNot1[["N"]])
    set(cohortData2, which(!wh1), totalOfLeadingBasedOn, totalBNot1)

    b <- cohortData2[, list(N = .N), by = pgdAndSc]
    b <- rep.int(b[["N"]], b[["N"]])
    GT1 <- (b > 1)
    if (any(GT1)) {
      pixelGroupData2List <- list()
      cohortData2[GT1, speciesProportion := sum(B, na.rm = TRUE) / totalB[1], by = pgdAndSc]
      cohortData2[!GT1, speciesProportion := B / totalB]
      # pixelGroupData2List[[2]] <- cohortData2[!GT1]
      # pixelGroupData2 <- rbindlist(pixelGroupData2List)
    } else {
      # cols <- c(pixelGroupColName, "speciesCode", "speciesProportion")
      set(cohortData2, NULL, "speciesProportion", cohortData2[[leadingBasedOn]] /
            cohortData2[[totalOfLeadingBasedOn]])
      # pixelGroupData2[[NROW(pixelGroupData2) + 1]] <- cohortData2[!GT1, ..cols]
    }
    pixelGroupData2 <- cohortData2
    systimePost2 <- Sys.time()
  }

  if (isTRUE(doAssertion)) {
    ## slower -- older, but simpler Eliot June 5, 2019
    ## TODO: these algorithm tests should be deleted after a while. See date on prev line.
    if (!exists("oldAlgoVTM", envir = .pkgEnv)) .pkgEnv$oldAlgoVTM <- 0
    if (!exists("newAlgoVTM", envir = .pkgEnv)) .pkgEnv$newAlgoVTM <- 0
    .pkgEnv$oldAlgoVTM <- .pkgEnv$oldAlgoVTM + (systimePost1 - systimePre1)
    .pkgEnv$newAlgoVTM <- .pkgEnv$newAlgoVTM + (systimePost2 - systimePre2)
    message("LandR::vegTypeMapGenerator: new algo ", .pkgEnv$newAlgoVTM)
    message("LandR::vegTypeMapGenerator: old algo ", .pkgEnv$oldAlgoVTM)
    setorderv(pixelGroupData2, pgdAndSc)
    whNA <- unique(unlist(sapply(pixelGroupData2, function(xx) which(is.na(xx)))))
    pixelGroupData1 <- pixelGroupData1[!pixelGroupData2[whNA], on = pgdAndSc]
    setkeyv(pixelGroupData1, pgdAndSc)
    setkeyv(pixelGroupData2, pgdAndSc)
    aa <- pixelGroupData1[pixelGroupData2, on = pgdAndSc]
    if (!isTRUE(all.equal(aa[["speciesProportion"]], aa[["i.speciesProportion"]]))) {
      stop("Old algorithm in vegMapGenerator is different than new map")
    }
  }

  if (algo == 1) {
    pixelGroupData <- pixelGroupData1
    rm(pixelGroupData1)
  } else if (algo == 2) {
    pixelGroupData <- pixelGroupData2
    rm(pixelGroupData2)
  }

  ########################################################
  #### Determine "mixed"
  ########################################################
  if (mixedType == 1) {
    ## create "mixed" class #    -- Eliot May 28, 2019 -- faster than previous below
    ## 1. anything with >= vegLeadingProportion is "pure"
    ## 2. sort on pixelGroup and speciesProportion, reverse so that 1st row of each pixelGroup is the largest
    ## 3. Keep only first row in each pixelGroup
    ## 4. change column names and convert pure to mixed ==> mixed <- !pure
    pixelGroupData3 <- pixelGroupData[, list(
      pure = speciesProportion >= vegLeadingProportion,
      speciesCode, pixelGroup, speciesProportion
    )]
    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName]
    pixelGroupData3[pure == FALSE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    pixelGroupData3[, pure := !pure]
    setnames(pixelGroupData3, "pure", "mixed")

    ## Old algorithm for above, this is ~43 times slower
    # a2 <- Sys.time()
    # pixelGroupData2 <- pixelGroupData[, list(mixed = all(speciesProportion < vegLeadingProportion),
    #                                          leading = speciesCode[which.max(speciesProportion)]),
    #                                   by = pixelGroupColName]
    # pixelGroupData2[mixed == TRUE, leading := "Mixed"]
    # b2 <- Sys.time()
  } else if (mixedType == 2) {
    if (!sppEquivCol %in% colnames(sppEquiv)) {
      stop(sppEquivCol, " is not in sppEquiv. Please pass an existing sppEquivCol")
    }

    sppEq <- data.table(sppEquiv[[sppEquivCol]], sppEquiv[["Type"]])

    names(sppEq) <- c("speciesCode", "Type")
    setkey(pixelGroupData, speciesCode)

    # don't need all columns now
    colsToDelete <- c("rasterToMatch", leadingBasedOn, totalOfLeadingBasedOn)
    colsToDelete <- colsToDelete[colsToDelete %in% colnames(pixelGroupData)]
    set(pixelGroupData, NULL, colsToDelete, NULL)
    pixelGroupData3 <- merge(pixelGroupData, sppEq[!duplicated(sppEq)], all.x = TRUE)
    setkeyv(pixelGroupData, pgdAndSc)

    setkeyv(pixelGroupData3, pgdAndSc)
    mixedType2Condition <- quote(Type == "Deciduous" &
                                   speciesProportion < vegLeadingProportion &
                                   speciesProportion > 1 - vegLeadingProportion)
    pixelGroupData3[, mixed := FALSE]

    if (algo == 2 || isTRUE(doAssertion)) {
      b <- pixelGroupData3[, list(N = .N), by = pixelGroupColName]
      b <- rep.int(b[["N"]], b[["N"]])
      GT1 <- b > 1


      pgd3GT1 <- pixelGroupData3[GT1]
      pgd3NGT1 <- pixelGroupData3[!GT1]

      pgd3GT1[eval(mixedType2Condition), mixed := TRUE, by = pixelGroupColName]
      pgd3GT1[, mixed := any(mixed), by = pixelGroupColName]
      pixelGroupData3 <- rbindlist(list(pgd3NGT1, pgd3GT1))
    } else {
      pixelGroupData3[eval(mixedType2Condition), mixed := TRUE, by = pixelGroupColName]
      pixelGroupData3[, mixed := any(mixed), by = pixelGroupColName]
    }

    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    set(pixelGroupData3, NULL, "Type", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName] ## sp. w/ highest prop. per pixelGroup
    pixelGroupData3[mixed == TRUE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    set(pixelGroupData3, NULL, "leading", factor(pixelGroupData3[["leading"]]))
  } else {
    stop("invalid mixedType! Must be one of '1' or '2'.")
  }

  cols <- c(pixelGroupColName, "leading")
  xx <- pixelGroupData3[, ..cols][x, on = pixelGroupColName]

  return(xx)
}



#' @importFrom SpaDES.tools inRange
preambleVTG <- function(x, vegLeadingProportion, doAssertion, nrowCohortData) {
  if (!inRange(vegLeadingProportion, 0, 1)) {
    stop("vegLeadingProportion must be a proportion")
  }


  leadingBasedOn <- if ("B" %in% colnames(x)) {
    message("Using B to derive leading type")
    "B"
  } else if ("cover" %in% colnames(x)) {
    message("Using cover to derive leading type, as there is no B column")
    "cover"
  } else {
    stop("x must have either B or cover to determine leading species/type")
  }
  if (!nrowCohortData > 0) stop("cohortData is empty")

  if (isTRUE(doAssertion)) {
    message("LandR::vegTypeMapGenerator: NROW(x) == ", nrowCohortData)
  }

  leadingBasedOn
}

# derived from plyr::mapvalues
mapvalues2 <- function(x, from, to) { #
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}
