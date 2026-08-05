utils::globalVariables(c(
  ".", "..cols", "..colsToSubset", ".I", ":=", "..groupVar",
  "age", "age2", "allValid", "aNPPAct", "cover", "coverOrig", "cumCells", "drawTarget",
  "dist", "distCells", "ecoregion", "ecoregionGroup",
  "hasBadAge", "imputedAge", "initialEcoregion", "initialEcoregionCode", "initialPixels",
  "lcc", "maxAge", "maxANPP", "maxB", "maxB_eco", "mortality", "nCells", "new", "newPossLCC",
  "noPixels",
  "oldSumB", "ord", "outBiomass", "oldEcoregionGroup", "selfLcc",
  "pixelGroup2", "pixelIndex", "pixels", "planted", "Provenance", "possERC",
  "speciesposition", "speciesGroup", "speciesInt", "state", "sumB",
  "temppixelGroup", "toDelete", "totalBiomass", "totalBiomass2", "totalCover",
  "uniqueCombo", "uniqueComboByRow", "uniqueComboByPixelIndex", "V1", "valid", "year"
))

#' Add cohorts to `cohortData` and `pixelGroupMap`
#'
#' This is a wrapper for `generatePixelGroups`, `initiateNewCohort` and updates to
#' `pixelGroupMap` via assignment to new `pixelIndex` values in `newPixelCohortData`.
#' By running these all together, there is less chance that they will diverge.
#' There are some checks internally for consistency.
#'
#' Does the following:
#' 1. add *new cohort* (not survivor) data into `cohortData`;
#' 2. assign initial `B` and `age` for new cohort;
#' 3. assign the new `pixelGroup` to the pixels that have new cohort;
#' 4. update the `pixelGroup` map.
#'
#' Note that if `newPixelCohortData` is generated after a disturbance
#'   it must contain a `type` column indicating the origin of the cohorts
#'   (e.g. "survivor", "serotiny", "resprouting"). "Survivor" cohorts will
#'   not be added to the output objects, as they are assumed to be accounted
#'   for in the input `cohortData` and, correspondingly, `pixelGroupMap`.
#'
#' @template newPixelCohortData
#' @template cohortData
#' @template pixelGroupMap
#' @template currentTime
#' @template speciesEcoregion
#' @template cohortDefinitionCols
#'
#' @param treedFirePixelTableSinceLastDisp A data.table with at least 2 columns, `pixelIndex`
#'   and `pixelGroup`.
#'   This will be used in conjunction with `cohortData` and `pixelGroupMap`
#'   to ensure that everything matches correctly.
#' @template successionTimestep
#' @template verbose
#' @template initialB
#' @template doAssertion
#'
#' @return
#' A list of length 2, `cohortData` and `pixelGroupMap`, with
#' `newPixelCohortData` inserted.
#'
#' @export
#' @rdname updateCohortData
updateCohortData <- function(
  newPixelCohortData,
  cohortData,
  pixelGroupMap,
  currentTime,
  speciesEcoregion,
  treedFirePixelTableSinceLastDisp = NULL,
  successionTimestep,
  cohortDefinitionCols = LandR::cohortDefinitionCols(),
  initialB = 10,
  verbose = getOption("LandR.verbose", TRUE),
  doAssertion = getOption("LandR.assertions", TRUE)
) {
  maxPixelGroup <- as.integer(maxFn(pixelGroupMap))

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
    ## Deal with pixels on the map that have no pixelGroup -- these are burned
    ## pixels --> the entirely newly regenerated pixels do not require a
    ## re-pixelGroupMaping  -- can just add to existing pixelGroup values
    if (verbose > 0) {
      message(cli::col_green(
        "  Regenerating only burnt pixels with no survivors (i.e. resprouting & serotiny)"
      ))
    }

    ## NOTE: no B in this columnsForPG b/c they all have zero
    columnsForPG <- LandR::columnsForPixelGroups()[-which(LandR::columnsForPixelGroups() == "B")]
    cd <- newPixelCohortData[, c("pixelIndex", columnsForPG), with = FALSE]
    newPixelCohortData[,
      pixelGroup := generatePixelGroups(cd, maxPixelGroup = maxPixelGroup, columns = columnsForPG)
    ] # ,
    # successionTimestep = successionTimestep)

    ## Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(
      newPixelCohortData[, c("pixelIndex", "pixelGroup")],
      by = c("pixelIndex")
    )
  } else {
    ## This is for situations where there are some empty pixels being filled,
    ## and some occupied pixels getting infilling. This requires a wholesale re-pixelGroup
    if (verbose > 0) {
      message(cli::col_green(
        "  Regenerating open and pixels with B (likely after seed dispersal, or partial mortality following disturbance)"
      ))
    }

    pixelIndex <- which(pixelGroupMap[] %in% cohortData$pixelGroup)

    ## remove unnecessary columns before making cohortDataLong
    suppressWarnings(set(cohortData, j = "prevMortality", value = NULL))
    suppressWarnings(set(newPixelCohortData, j = "year", value = NULL))

    cohortDataPixelIndex <- data.table(
      pixelIndex = pixelIndex,
      pixelGroup = pixelGroupMap[][pixelIndex]
    )
    cdLong <- cohortDataPixelIndex[cohortData, on = "pixelGroup", allow.cartesian = TRUE]
    if (!is.null(newPixelCohortData$type)) {
      ## survivor cohorts are assumed to be in cohort data already.
      newPixelCohortData <- newPixelCohortData[type != "survivor"]
    }
    cohorts <- rbindlist(list(cdLong, newPixelCohortData), use.names = TRUE, fill = TRUE)

    columnsForPG <- LandR::columnsForPixelGroups()
    cd <- cohorts[, c("pixelIndex", columnsForPG), with = FALSE]
    cohorts[, pixelGroup := generatePixelGroups(cd, maxPixelGroup = 0L, columns = columnsForPG)]

    ## Bring to pixelGroup level -- this will squash the data.table
    if (is.null(cohorts[["sumB"]])) {
      cohorts[, sumB := sum(B, na.rm = TRUE), by = pixelGroup]
    }
    ## Old way that does not preserve additional 'unknown' columns in cohortData
    # allCohortData <- cohorts[ , .(ecoregionGroup = ecoregionGroup[1],
    #                               mortality = mortality[1],
    #                               aNPPAct = aNPPAct[1],
    #                               sumB = sumB[1]),
    #                           by = uniqueCohortDefinition]

    ## newer way that will potentially conflict with LandR.CS due to differing aNPP
    colsToSubset <- setdiff(colnames(cohortData), c("pixelIndex"))
    allCohortData <- cohorts[
      !duplicated(cohorts[, .(pixelGroup, speciesCode, ecoregionGroup, age)]),
      ..colsToSubset
    ]

    theNewOnes <- is.na(allCohortData$B)
    cohortData <- allCohortData[!theNewOnes]
    newPixelCohortData <- allCohortData[theNewOnes]

    ## Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
    pixelsToChange <- unique(cohorts[, c("pixelIndex", "pixelGroup")], by = c("pixelIndex"))
  }

  ## update pixelGroupMap
  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (doAssertion) {
    if (!isTRUE(all(pixelsToChange$pixelGroup == pixelGroupMap[][pixelsToChange$pixelIndex]))) {
      stop("pixelGroupMap and newPixelCohortData$pixelGroupMap don't match in updateCohortData fn")
    }
  }

  ## give B in pixels that have serotiny/resprouting
  # newPixelCohortData[, sumB := asInteger(sum(B, na.rm = TRUE)), by = pixelGroup]

  ## Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone) -----------
  cohortData <- .initiateNewCohorts(
    newPixelCohortData,
    cohortData,
    pixelGroupMap,
    currentTime = currentTime,
    speciesEcoregion = speciesEcoregion,
    successionTimestep = successionTimestep,
    initialB = initialB
  )

  outs <- rmMissingCohorts(
    cohortData,
    pixelGroupMap,
    cohortDefinitionCols = LandR::cohortDefinitionCols()
  )

  if (!is.null(outs$cohortData$sumB)) {
    outs$cohortData[, sumB := NULL]
  }

  assertCohortData(
    outs$cohortData,
    outs$pixelGroupMap,
    cohortDefinitionCols = LandR::cohortDefinitionCols(),
    doAssertion = doAssertion,
    verbose = verbose
  )

  if (verbose > 0) {
    nPixForest <- sum(!is.na(outs$pixelGroupMap[]))
    nPixGrps <- length(unique(outs$cohortData$pixelGroup))
    nPixNoPixGrp <- sum(outs$pixelGroupMap[] == 0, na.rm = TRUE)
    nPixTreed <- sum(outs$pixelGroupMap[] != 0, na.rm = TRUE)

    nDigits <- max(nchar(c(nPixForest, nPixGrps, nPixNoPixGrp))) + 3
    message(cli::col_magenta(
      "NUMBER OF FORESTED PIXELS          :",
      paddedFloatToChar(nPixForest, padL = nDigits, pad = " ")
    ))
    message(cli::col_magenta(
      "NUMBER OF PIXELS WITH TREES        :",
      paddedFloatToChar(nPixTreed, padL = nDigits, pad = " ")
    ))
    message(cli::col_magenta(
      "NUMBER OF UNIQUE PIXELGROUPS       :",
      paddedFloatToChar(nPixGrps, padL = nDigits, pad = " ")
    ))
    message(cli::col_magenta(
      "NUMBER OF PIXELS WITH NO PIXELGROUP:",
      paddedFloatToChar(nPixNoPixGrp, padL = nDigits, pad = " ")
    ))
  }

  return(list(cohortData = outs$cohortData, pixelGroupMap = outs$pixelGroupMap))
}

#' Initiate new cohorts
#'
#' Calculate new values for `B`, add `age`, then `rbindlist` this
#' with `cohortData`.
#'
#' @template newPixelCohortData
#' @template cohortData
#' @template cohortDefinitionCols
#' @template initialB
#'
#' @return A `data.table` with a new `rbindlist`ed `cohortData`
#'
#' @rdname updateCohortData
.initiateNewCohorts <- function(
  newPixelCohortData,
  cohortData,
  pixelGroupMap,
  currentTime,
  cohortDefinitionCols = LandR::cohortDefinitionCols(),
  speciesEcoregion,
  successionTimestep,
  initialB = 10
) {
  ## get spp "productivity traits" per ecoregion/present year
  ## calculate maximum B per ecoregion, join to new cohort data
  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, pixelGroup := as.vector(pixelGroupMap[])[pixelIndex]]
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
  specieseco_current <- specieseco_current[, .(speciesCode, maxANPP, maxB, ecoregionGroup)] |>
    setkey(speciesCode, ecoregionGroup)

  ## Note that after the following join, some cohorts will be lost due to lack of
  ##  parameters in speciesEcoregion. These need to be modified in pixelGroupMap.
  # missingNewPixelCohortData <- newPixelCohortData[!specieseco_current,
  #                                                 on = uniqueSpeciesEcoregionDefinition]
  specieseco_current <- specieseco_current[!is.na(maxB)]
  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]
  newPixelCohortData <- specieseco_current[
    newPixelCohortData,
    on = uniqueSpeciesEcoregionDefinition
  ]
  newPixelCohortData <- newPixelCohortData[!is.na(maxB)]

  if (any(newPixelCohortData$age > 1)) {
    stop("newPixelCohortData should only have new cohorts aged 1")
  }
  set(newPixelCohortData, NULL, "age", 1L) ## set age to 1

  ## Ceres: this was causing new cohorts to be initialized with maxANPP.
  ## instead, remove column and calculate total biomass of older cohorts in cohortData
  # set(newPixelCohortData, NULL, "sumB", 0L)
  suppressWarnings(set(newPixelCohortData, j = "sumB", value = NULL))

  cohortData[age >= successionTimestep, oldSumB := sum(B, na.rm = TRUE), by = "pixelGroup"]

  newPixelCohortData <- unique(cohortData[, .(pixelGroup, oldSumB)], by = "pixelGroup")[
    newPixelCohortData,
    on = "pixelGroup"
  ]
  ## using set() is faster than [:=]
  set(newPixelCohortData, which(is.na(newPixelCohortData$oldSumB)), "oldSumB", 0)
  setnames(newPixelCohortData, "oldSumB", "sumB")
  set(cohortData, NULL, "oldSumB", NULL)

  if ("B" %in% names(newPixelCohortData)) {
    newPixelCohortData[, B := NULL]
  }

  if (isTRUE(is.na(initialB)) || is.null(initialB)) {
    set(
      newPixelCohortData,
      NULL,
      "B",
      asInteger(pmax(
        1,
        newPixelCohortData$maxANPP *
          exp(-1.6 * newPixelCohortData$sumB / newPixelCohortData$maxB_eco)
      ))
    )
    set(
      newPixelCohortData,
      NULL,
      "B",
      asInteger(pmin(newPixelCohortData$maxANPP, newPixelCohortData$B))
    )
  } else {
    ## 2022-02: change to 10 - maxANPP is unrealistic, particularly with
    ##          high maxANPP needed to produce realistic growth curves
    set(newPixelCohortData, NULL, "B", asInteger(initialB))
  }

  newPixelCohortData <- newPixelCohortData[, .(
    pixelGroup,
    ecoregionGroup,
    speciesCode,
    age,
    B,
    mortality = 0L,
    aNPPAct = 0L
  )]

  if (getOption("LandR.assertions")) {
    if (NROW(unique(newPixelCohortData, by = cohortDefinitionCols)) != NROW(newPixelCohortData)) {
      stop("Duplicated new cohorts in a pixelGroup. Please debug LandR:::.initiateNewCohorts")
    }
  }

  ## keep track of new cohorts to ensure they get the correct ecoregionGroup (from existing ones)
  set(cohortData, NULL, "new", FALSE)
  set(newPixelCohortData, NULL, "new", TRUE)
  cohortData <- rbindlist(list(cohortData, newPixelCohortData), fill = TRUE, use.names = TRUE)
  if (
    NROW(unique(cohortData, by = c("pixelGroup"))) !=
      NROW(unique(cohortData, by = c("pixelGroup", "ecoregionGroup")))
  ) {
    message("Found pixelGroup with multiple ecoregionGroups when initiating cohorts.")
    message("Adjusting new ecoregionGroups to match those of existing pixelGroups.")
    cohortData[,
      ecoregionGroup := unique(.SD[new == FALSE, ecoregionGroup]),
      by = "pixelGroup",
      .SDcols = c("new", "ecoregionGroup")
    ] ## TODO: very slow!!
  }
  set(cohortData, NULL, "new", NULL)
  set(newPixelCohortData, NULL, "new", NULL)

  ## recalculate sumB
  cohortData[, sumB := asInteger(sum(B, na.rm = TRUE)), by = "pixelGroup"]

  assertCohortDataERG(cohortData)

  return(cohortData)
}

#' Remove missing cohorts from `cohortData` based on `pixelGroupMap`
#'
#' @template cohortData
#'
#' @template cohortDefinitionCols
#'
#' @template pixelGroupMap
#'
#' @template doAssertion
#'
#' @return
#' A `list` with 2 `data.table` objects, `cohortData` and `pixelGroupMap`,
#' each updated based on missing `pixelGroups` in the other.
#'
#' @export
rmMissingCohorts <- function(
  cohortData,
  pixelGroupMap,
  cohortDefinitionCols = LandR::cohortDefinitionCols(),
  doAssertion = getOption("LandR.assertions", TRUE)
) {
  pgmValues <- data.table(
    pixelGroup = as.vector(pixelGroupMap[]),
    pixelIndex = seq(ncell(pixelGroupMap))
  )

  pgmVals <- pgmValues[!is.na(pixelGroup), ] ## na.omit() doesn't omit NaN, which terra returns
  pgmVals <- pgmVals[pixelGroup > 0]
  whPgsStillInCDGoneFromPGM <- !cohortData$pixelGroup %in% pgmVals$pixelGroup
  pgsStillInCDGoneFromPGM <- cohortData[whPgsStillInCDGoneFromPGM, ]
  whPgsStillInPGMGoneFromCD <- !pgmVals$pixelGroup %in% cohortData$pixelGroup
  pgsStillInPGMGoneFromCD <- pgmVals[whPgsStillInPGMGoneFromCD, ]

  ## REMOVE lines in cohortData that are no longer in the pixelGroupMap
  cohortData <- cohortData[!pixelGroup %in% pgsStillInCDGoneFromPGM$pixelGroup]
  ## REMOVE pixels in pixelGroupMap that are no longer in the cohortData
  pixelGroupMap[pgsStillInPGMGoneFromCD$pixelIndex] <- NA

  assertCohortData(
    cohortData,
    pixelGroupMap,
    message = "rmMissingCohorts",
    cohortDefinitionCols = LandR::cohortDefinitionCols(),
    doAssertion = doAssertion
  )

  assertCohortDataERG(cohortData, doAssertion = doAssertion)

  return(list(cohortData = cohortData, pixelGroupMap = pixelGroupMap))
}

#' Add the correct `pixelGroups` to a `pixelDataTable` object
#'
#' Generates unique groupings of a `data.table` object where one or more rows can
#' all belong to the same `pixelIndex`. Pixel groups will be identical pixels based
#' on unique combinations of `columns`.
#'
#' @param pixelDataTable  A `data.table` with column-based descriptions.
#'   Must have a column called `pixelIndex`, which allows for multiple rows to be associated
#'   with a single pixel.
#' @param maxPixelGroup A length 1 numeric indicating the current maximum `pixelGroup` value;
#'    the `pixelGroup` numbers returned will start at `maxPixelGroup + 1`.
#' @param columns A character vector of column names to use as part of the generation of unique
#'   combinations of features. Default is `c("ecoregionGroup", "speciesCode", "age", "B")`
#'
#' @return
#' Returns a vector of `pixelGroup` in the original order of the input `pixelDataTable`.
#' This should likely be added to the `pixelDataTable` object immediately.
#'
#' @export
generatePixelGroups <- function(
  pixelDataTable,
  maxPixelGroup,
  columns = c("ecoregionGroup", "speciesCode", "age", "B")
) {
  columnsOrig <- columns
  columns <- columns[columns %in% names(pixelDataTable)]
  columns2 <- paste0(columns, "2")
  if (!all(columns == columnsOrig)) {
    message(
      "Creating pixelGroup values, but not using all columns requested. Only using, ",
      paste(columns, collapse = ", "),
      " instead of ",
      paste(columnsOrig, collapse = ", ")
    )
  }

  pcd <- pixelDataTable # no copy -- just for simpler name

  if (getOption("LandR.assertions")) {
    pcdOrig <- data.table::copy(pcd)
  }
  ## concatenate within rows:
  ## e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
  pcd[, uniqueComboByRow := do.call(paste, as.list(.SD)), .SDcols = columns]

  ## concatenate within pixelIndex
  pcd[, c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"), by = "pixelIndex"]
  pcd[, c("pixelGroup") := as.integer(maxPixelGroup) + as.integer(factor(uniqueComboByPixelIndex))]

  ## old algorithm:
  if (isTRUE(getOption("LandR.assertions"))) {
    ## prepare object 1 (pcd) for checking below
    pcd[, ord := seq_len(.N)]
    setorderv(pcd, c("pixelIndex"))
    uniqPG <- unique(pcd$pixelGroup)
    pcd[,
      pixelGroup2 := mapvalues2(pixelGroup, from = uniqPG, to = as.character(seq_along(uniqPG)))
    ]
    # pcd[, pixelGroup2 := mapvalues(pixelGroup, from = unique(pixelGroup),
    #                                to = as.character(seq_along(unique(pixelGroup))))]
    setorderv(pcd, "ord")

    pcdOld <- data.table::copy(pcdOrig)

    # Convert to unique numeric
    pcdOld[,
      c(columns2) := lapply(.SD, function(x) {
        a <- as.integer(factor(x))
      }),
      .SDcols = columns
    ]

    ## concatenate within rows:
    ## e.g., ecoregionCode_speciesCode_age_biomass or 647_11_Abie_sp_100_2000
    pcdOld[,
      uniqueComboByRow := as.integer(factor(do.call(paste, as.list(.SD)))),
      .SDcols = columns2
    ]

    ## concatenate within pixelIndex
    pcdOld[,
      c("uniqueComboByPixelIndex") := paste(uniqueComboByRow, collapse = "__"),
      by = "pixelIndex"
    ]
    pcdOld[,
      c("pixelGroup") := as.integer(maxPixelGroup) + as.integer(factor(uniqueComboByPixelIndex))
    ]
    ## prepare object 2 (pcdOld) for checking below
    pcdOld[, ord := seq_len(.N)]
    setorderv(pcdOld, c("pixelIndex"))

    uniqPG <- unique(pcdOld$pixelGroup)
    pcdOld[,
      pixelGroup2 := mapvalues2(pixelGroup, from = uniqPG, to = as.character(seq_along(uniqPG)))
    ]

    setorderv(pcdOld, "ord")

    ## the check
    if (!identical(pcdOld$pixelGroup2, pcd$pixelGroup2)) {
      stop("new generatePixelGroups algorithm failing")
    }
  }

  return(pcd$pixelGroup)
}

#' Pull out the values from `speciesEcoregion` table for current time
#'
#' @template speciesEcoregion
#' @template currentTime
#'
#' @return
#' The `speciesEcoregion` input object, but with data from only one year, the year
#' that is less than or equal to the `currentTime`
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

#' The columns in a `cohortData` that define "unique"
#'
#' If two pixels have identical values in all of these columns, they are the same `pixelGroup`.
#'
#' @export
#' @rdname uniqueDefinitions
uniqueCohortDefinition <- c("pixelGroup", "speciesCode", "age", "B")

#' @export
#' @rdname uniqueDefinitions
uniqueSpeciesEcoregionDefinition <- c("speciesCode", "ecoregionGroup")

#' Summary for `cohortData`
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
  message(cli::col_magenta(
    "Pixels with non-NA cover:, ",
    cohortData[!is.na(cover), length(unique(pixelIndex))]
  ))
}

#' @keywords internal
.cohortMessages <- function(cohortData, val) {
  out <- list()
  if (val %in% colnames(cohortData)) {
    pixelsNA <- NROW(cohortData[is.na(get(val)), unique("pixelIndex"), with = FALSE])
    message(cli::col_magenta("Pixels with missing", val, ":", format(pixelsNA, big.mark = ",")))
    pixelsZero <- NROW(cohortData[, all(get(val) == 0), by = "pixelIndex"][get("V1") == TRUE])
    message(cli::col_magenta("Pixels with all(", val, " == 0): ", format(pixelsZero, big.mark = ",")))
    pixelsBiomassNonZero <- NROW(cohortData[, any(get(val) > 0), by = "pixelIndex"][
      get("V1") == TRUE
    ])
    message(cli::col_magenta(
      "Pixels with all(",
      val,
      " > 0): ",
      format(pixelsBiomassNonZero, big.mark = ",")
    ))
    out <- list(
      pixelsNA = pixelsNA,
      pixelsZero = pixelsZero,
      pixelsBiomassNonZero = pixelsBiomassNonZero
    )
  }
  return(invisible(out))
}

#' Convert Land Cover Classes (LCC) to another value in its neighbourhood
#'
#' This searches around the pixels on `rstLCC` that have `classesToReplace` for the
#' nearest neighbouring Land Cover Class other than those indicated in
#' `classesToReplace`, and assigns the cohorts in those pixels the corresponding new
#' `ecoregionCode`. Replacement codes must already exist in `availableERC_by_Sp`,
#' i.e., no new `ecoregionCode` values are created. See Details.
#'
#' @details
#' This function is designed to be used in highly constrained situations, where it is not
#' just replacing a Land Cover Class by a neighbouring Land Cover Class. But it can
#' be used for the simpler cases of simply replacing a Land Cover Class.
#'
#' Each pixel whose value is in `classesToReplace` is assigned an available
#' `initialEcoregionCode` taken from its neighbourhood. For each candidate land-cover
#' class present in `rstLCC`, the distance from every pixel to the nearest cell of that
#' class is computed with a single vectorized [terra::distance()] call, and each unwanted
#' pixel keeps only those classes whose implied `initialEcoregionCode` is available for
#' that pixel's ecoregion (and all of its `speciesCode`s). One of the surviving classes is
#' then drawn with probability proportional to how many cells of that class fall inside
#' the pixel's neighbourhood -- the smallest window reaching that pixel's nearest
#' available class. `method` decides only where the draw comes from:
#'
#' * `"nearestWeighted"` (default) derives it from the pixel's own ground position, so the
#'   result is deterministic: it needs no `set.seed()`, is stable under [reproducible::Cache()]
#'   (which does not key on RNG state), and, because the key is the cell *centre* rather than
#'   the cell *index*, a grid-aligned crop reproduces its parent raster cell for cell --
#'   so a small development subset agrees with the full run. Reprojecting or changing
#'   resolution changes the pixels, and so changes the result.
#' * `"nearestRandom"` draws from the RNG instead, so it varies with `set.seed()`. Use it
#'   when replicates should differ; note that `Cache()` will replay a single draw unless
#'   the seed is part of the cache key.
#'
#' Both are `O(nClasses)` distance transforms, independent of the geometry of the
#' `classesToReplace` blobs. They replace an earlier iterative `spread2()`-based search
#' whose run time grew with the square of the radius of the largest contiguous block of
#' `classesToReplace`, which could take hours -- or never finish -- on large study areas
#' with big lakes or burns, or an irregular masked boundary.
#'
#' The abundance weighting restores the behaviour of that former search, which sampled
#' among all valid cells within the radius at which it first found one and so favoured
#' classes in proportion to their local abundance. The window is rectangular, matching the
#' Chebyshev neighbourhood of `spread2(directions = 8)`. What it does *not* reproduce is
#' that implementation's exact output: the radius here comes from a Euclidean distance
#' transform rather than from counting spread iterations, and the draws differ. Runs prior
#' to LandR 1.2.0.9004 cannot be reproduced bit-for-bit by any current `method`.
#'
#' LandR 1.2.0.9004's deterministic nearest-class rule, which broke ties to the lowest
#' land-cover class, has been removed. Ties are common -- 35-41% of unwanted pixels on real
#' landscapes -- so always taking the lowest code biased the result toward low-numbered
#' classes, which in the Canada LCC coding means toward sparse, non-forest cover: it
#' assigned shrubs 1.69x as often as the previous implementation, and broadleaf 0.58x and
#' mixedwood 0.28x as often.
#'
#' @param pixelClassesToReplace Deprecated. Use `classesToReplace`
#'
#' @param classesToReplace Integer vector of classes that are are to be replaced,
#'     e.g., 34, 35, 36 on LCC2005, which are burned young, burned 10 year, and cities.
#'
#' @param rstLCC LCC raster, e.g., LCC2005
#'
#' @param theUnwantedPixels An optional vector of pixel IDs that need to be changed.
#'   If not provided, then pixels to change will be taken from the match between
#'   `availableERC_by_Sp` and `classesToReplace`. Supplying this allows
#'   the user to only replace some of the pixels with a given class.
#'
#' @param ecoregionGroupVec Deprecated. Use `availableERC_by_Sp`
#'
#' @param speciesEcoregion Deprecated. Use `availableERC_by_Sp`
#'
#' @param availableERC_by_Sp A `data.table` or `data.frame` with 3 columns:
#'   `speciesCode`, `initialEcoregionCode` and `pixelIndex`.
#'   `pixelIndex` is the pixel id for each line in the `data.table`;
#'   `speciesCode` is the species name in the pixel (can have more than one
#'     species per pixel, so multiple rows per pixel); and,
#'   `initialEcoregionCode` is the unique codes that are "available" to be
#'     used as a replacement for `classesToReplace`. `initialEcoregionCode`
#'     must be a character vector, with one or no "_" used as a separator, with the last
#'     component being the Land Cover Class that matches `classesToReplace`, e.g.,
#'     `"242_18"`. If there is no "_" in this code, then the codes must match the
#'     `classesToReplace` exactly, e.g., `"11"`.
#'     If `pixelIndex` is missing, the function will fill it
#'     with `seq(ncell(rstLCC))`. If `speciesCode` is missing, the function
#'     will replace it with a dummy value (`"allSpecies"`).
#'
#' @template doAssertion
#'
#' @param method Character; where the draw among an unwanted pixel's available land-cover
#'   classes comes from. Both options weight the classes by their abundance in the pixel's
#'   neighbourhood and differ only in reproducibility: `"nearestWeighted"` (default) keys
#'   the draw on the pixel's ground position, so it is deterministic and seed-free, while
#'   `"nearestRandom"` draws from the RNG and so varies with `set.seed()`. See Details.
#'
#' @return
#' A `data.table` with three columns, `newPossLCC`, `pixelIndex` and `ecoregionGroup`.
#' This represents the new codes to used in the `pixelIndex` locations.
#' These should have no values overlapping with `classesToReplace`.
#' `newPossLCC` is the assigned land-cover class itself, i.e. `ecoregionGroup` without
#' its ecoregion prefix; it is `NA` for pixels that could not be assigned.
#'
#' @author Eliot McIntire
#' @export
convertUnwantedLCC <- function(
  classesToReplace = 34:36,
  rstLCC,
  availableERC_by_Sp,
  theUnwantedPixels,
  ecoregionGroupVec,
  speciesEcoregion,
  pixelClassesToReplace,
  doAssertion = getOption("LandR.assertions", TRUE),
  method = c("nearestWeighted", "nearestRandom")
) {
  method <- match.arg(method)

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
    numCharEcoregion <- if (nrow(availableERG2) > 0) {
      max(nchar(availableERG2$ecoregion), na.rm = TRUE)
    } else {
      ## above gives -Inf when no classes to replace, so use zero instead
      0
    }
    numCharLCCCodes <- numCharIEC - numCharEcoregion - 1 # minus 1 for the dash
    availableERG2[, `:=`(ecoregion = NULL)]
  } else {
    numCharLCCCodes <- numCharIEC
  }

  ## Assign each unwanted pixel an available `initialEcoregionCode` from its neighbourhood.
  ## For each candidate land-cover class present in `rstLCC` (i.e. not a class to replace),
  ## the distance from every pixel to the nearest cell of that class is computed with a
  ## single vectorized `terra::distance()`; each unwanted pixel then keeps those classes
  ## whose implied `possERC` is available for its ecoregion (and all of its `speciesCode`s),
  ## and one is drawn weighted by local abundance -- `method` only decides whether the draw
  ## is keyed on the pixel's position or on the RNG. Both cost `O(nClasses)` distance transforms,
  ## independent of blob geometry -- unlike the former iterative `spread2()` search, whose
  ## cost grew with the square of the radius of the largest contiguous block of
  ## `classesToReplace`.
  if (length(theUnwantedPixels) == 0L) {
    out3 <- data.table(newPossLCC = NA, pixelIndex = NA, ecoregionGroup = NA)[!is.na(pixelIndex)]
  } else {
    message(
      "Converting ",
      length(theUnwantedPixels),
      " unwanted LCC pixels to a nearby available class, weighted by local abundance (",
      method,
      ")."
    )
    lccVals <- as.vector(rstLCC[])
    candClasses <- sort(unique(lccVals[!is.na(lccVals) & !lccVals %in% classesToReplace]))

    if (length(candClasses) == 0L) {
      out3 <- data.table(newPossLCC = NA, pixelIndex = theUnwantedPixels, ecoregionGroup = NA)
    } else {
      ## distance from each unwanted pixel to the nearest cell of each candidate class
      distByClass <- vapply(
        candClasses,
        function(cc) {
          m <- terra::setValues(terra::rast(rstLCC), ifelse(lccVals == cc, 1L, NA_integer_))
          terra::values(terra::distance(m))[theUnwantedPixels, 1]
        },
        numeric(length(theUnwantedPixels))
      )
      cand <- data.table(
        pixelIndex = rep(theUnwantedPixels, times = length(candClasses)),
        lcc = rep(candClasses, each = length(theUnwantedPixels)),
        dist = as.vector(distByClass)
      )[is.finite(dist)]

      ## attach each unwanted pixel's ecoregion + speciesCode(s)
      uwCols <- c("pixelIndex", "speciesCode", if (hasPreDash) "ecoregion")
      uwInfo <- unique(availableERC_by_Sp[pixelIndex %in% theUnwantedPixels, uwCols, with = FALSE])
      cand <- cand[uwInfo, on = "pixelIndex", allow.cartesian = TRUE]

      ## candidate ecoregion-group code, formed exactly as the former loop did
      if (hasPreDash) {
        cand[,
          possERC := paste0(
            ecoregion,
            "_",
            paddedFloatToChar(as.integer(lcc), padL = numCharLCCCodes, padR = 0)
          )
        ]
      } else {
        cand[, possERC := as.character(lcc)]
      }

      ## keep only available (speciesCode, possERC) combinations; a (pixel, class) is
      ## assignable only if valid for *all* of the pixel's speciesCodes (cf. possERCToRm)
      validKey <- paste(availableERG2$speciesCode, availableERG2$initialEcoregionCode)
      cand[, valid := paste(speciesCode, possERC) %in% validKey]
      ok <- cand[,
        list(allValid = all(valid), possERC = possERC[1L]),
        by = c("pixelIndex", "lcc", "dist")
      ][allValid == TRUE]

      chosen <- .chooseByLocalAbundance(
        ok,
        rstLCC = rstLCC,
        lccVals = lccVals,
        candClasses = candClasses,
        deterministic = identical(method, "nearestWeighted")
      )
      out3 <- chosen[, list(
        newPossLCC = as.integer(lcc),
        pixelIndex,
        ecoregionGroup = if (hasPreDash) possERC else as.integer(lcc)
      )]

      ## unwanted pixels with no available replacement anywhere -> NA
      missed <- setdiff(theUnwantedPixels, out3$pixelIndex)
      if (length(missed) > 0L) {
        out3 <- rbind(
          out3,
          data.table(newPossLCC = NA, pixelIndex = missed, ecoregionGroup = NA),
          fill = TRUE
        )
      }
      out3 <- unique(out3)
    }
  }

  if (doAssertion) {
    if (any(gsub(".*_", "", out3$ecoregionGroup) %in% classesToReplace)) {
      stop("classesToReplace were not fully removed")
    }
  }

  out3
}

## Choose, for each unwanted pixel, one of the available land-cover classes in its
## neighbourhood, with probability proportional to how many cells of that class the
## neighbourhood contains. The neighbourhood is the smallest window that reaches the
## pixel's nearest available class -- i.e. the window at which the former `spread2()`
## search would have stopped -- so this reproduces that search's abundance weighting
## without its O(radius^2) cost. Counts come from a summed-area table, so a window 1500
## cells across costs the same as one 3 cells across.
##
## `ok` is the table of available (pixelIndex, lcc, dist, possERC) combinations built by
## `convertUnwantedLCC()`; every candidate class is offered to every pixel it is available
## for, not only the nearest one, because a class whose nearest cell lies beyond the window
## radius can still have cells in the window's corners -- as it could under `spread2()`,
## whose `directions = 8` neighbourhood was itself square.
.chooseByLocalAbundance <- function(ok, rstLCC, lccVals, candClasses, deterministic) {
  nRow <- terra::nrow(rstLCC)
  nCol <- terra::ncol(rstLCC)
  nCand <- length(candClasses)

  ## The window has to be sized in *cells*, but `ok$dist` is in map units, and dividing it
  ## by `res(rstLCC)` is only exact for square cells and outright wrong for a lon/lat raster
  ## (terra::distance() returns great-circle metres there, while res() is in degrees). So
  ## redo the distance transforms on a 1 m-per-cell projected grid of the same dimensions:
  ## the distances then *are* cell counts, whatever `rstLCC`'s CRS and resolution.
  uwPixels <- sort(unique(ok$pixelIndex))
  unitRas <- terra::rast(
    nrows = nRow,
    ncols = nCol,
    xmin = 0,
    xmax = nCol,
    ymin = 0,
    ymax = nRow,
    crs = "EPSG:3978"
  )
  distCells <- vapply(
    candClasses,
    function(cc) {
      m <- terra::setValues(unitRas, ifelse(lccVals == cc, 1L, NA_integer_))
      terra::values(terra::distance(m))[uwPixels, 1]
    },
    numeric(length(uwPixels))
  )
  inCells <- data.table(
    pixelIndex = rep(uwPixels, times = nCand),
    lcc = rep(candClasses, each = length(uwPixels)),
    distCells = as.vector(distCells)
  )

  ## window half-widths: far enough to reach the pixel's nearest *available* class, but
  ## always at least the 8 adjacent cells (one `spread2()` iteration)
  rad <- ok[inCells, on = c("pixelIndex", "lcc"), nomatch = NULL][,
    list(rmin = min(distCells)),
    by = "pixelIndex"
  ]
  set(rad, NULL, "kx", pmin(nCol, pmax(1L, as.integer(ceiling(rad$rmin)))))
  set(rad, NULL, "ky", pmin(nRow, pmax(1L, as.integer(ceiling(rad$rmin)))))

  counts <- windowCountsByClassCpp(
    lccVals = as.integer(lccVals),
    candClasses = as.integer(candClasses),
    nrow = nRow,
    ncol = nCol,
    cells0 = as.integer(rad$pixelIndex - 1L),
    kx = rad$kx,
    ky = rad$ky
  )
  wts <- data.table(
    pixelIndex = rep(rad$pixelIndex, times = nCand),
    lcc = rep(candClasses, each = nrow(rad)),
    nCells = as.vector(counts),
    selfLcc = rep(as.integer(lccVals[rad$pixelIndex]), times = nCand)
  )
  ## a pixel is not its own neighbour (cf. the former `out[initialPixels != pixels]`)
  wts[!is.na(selfLcc) & lcc == selfLcc, nCells := nCells - 1L]
  set(wts, NULL, "selfLcc", NULL)

  ## the nearest available class always has >= 1 cell in the window (its nearest cell is
  ## within `rmin`, hence within `ceiling(rmin / res)` cells), so no pixel loses all weight
  ok <- ok[wts, on = c("pixelIndex", "lcc"), nomatch = NULL][nCells > 0L]
  setorderv(ok, c("pixelIndex", "lcc"))
  set(ok, NULL, "cumCells", ok[, cumsum(nCells), by = "pixelIndex"]$V1)

  ## One uniform draw per pixel, scaled to that pixel's total weight. `"nearestWeighted"`
  ## takes it from the pixel's ground position, so the result is reproducible without a seed
  ## and a grid-aligned crop reproduces its parent raster cell for cell; `"nearestRandom"`
  ## takes it from the RNG, so it varies with `set.seed()` across replicates.
  draws <- ok[, list(drawTarget = sum(nCells)), by = "pixelIndex"]
  u <- if (isTRUE(deterministic)) {
    xy <- terra::xyFromCell(rstLCC, draws$pixelIndex)
    resXY <- terra::res(rstLCC)
    pixelUnifCpp(x = xy[, 1], y = xy[, 2], resx = resXY[1], resy = resXY[2])
  } else {
    runif(nrow(draws))
  }
  set(draws, NULL, "drawTarget", u * draws$drawTarget)

  ok <- ok[draws, on = "pixelIndex"]
  setorderv(ok, c("pixelIndex", "lcc"))
  ok[cumCells >= drawTarget, .SD[1L], by = "pixelIndex"]
}


#' Assess non-forested pixels based on species cover data and land-cover
#'
#' @template speciesLayers
#'
#' @param omitNonTreedPixels logical. Should pixels with classes in `forestedLCCClasses` be
#'                           included as non-forested?
#'
#' @param forestedLCCClasses vector of forested land-cover classes in `rstLCC`
#'
#' @template rstLCC
#'
#' @return a logical vector of length `ncell(rstLCC)`
#'   where `TRUE` indicates non-forested pixels where there is no
#'   species cover data, or a non-forested land-cover class
#'
#' @export
nonForestedPixels <- function(speciesLayers, omitNonTreedPixels, forestedLCCClasses, rstLCC) {
  # pixelsToRm <- rowSums(!is.na(sim$speciesLayers[])) == 0 # keep
  pixelsToRm <- is.na(as.vector(speciesLayers[[1]][])) # seems to be OK because seem to be NA on each layer for a given pixel

  ## remove non-forested if asked by user
  if (omitNonTreedPixels) {
    if (is.null(forestedLCCClasses)) {
      stop(
        "No `forestedLCCClasses` provided, but `omitNonTreedPixels` is TRUE.\n",
        "Please provide a vector of forested classes in `forestedLCCClasses`"
      )
    }
    lccPixelsRemoveTF <- !(as.vector(rstLCC[]) %in% forestedLCCClasses)
    pixelsToRm <- lccPixelsRemoveTF | pixelsToRm
  }
  pixelsToRm
}

#' Generate template `cohortData` table
#'
#' Internal function used by [makeAndCleanInitialCohortData()].
#'
#' @param inputDataTable A `data.table` with columns described above.
#'
#' @param sppColumns A vector of the names of the columns in `inputDataTable` that
#'   represent percent cover by species, rescaled to sum up to 100%%.
#'
#' @param minCoverThreshold minimum total cover percentage necessary to consider the pixel
#'  vegetated, or a cohort present in a pixel.
#'
#' @template doAssertion
#'
#' @param rescale Logical. If `TRUE`, the default, cover for each species will be rescaled
#'   so all cover in `pixelGroup` or pixel sums to 100.
#'
#' @return `cohortData` (`data.table`) with attribute `"imputedPixID"`
#'
#' @keywords internal
.createCohortData <- function(
  inputDataTable,
  sppColumns,
  minCoverThreshold = 5,
  doAssertion = getOption("LandR.assertions", TRUE),
  rescale = TRUE
) {
  newCoverColNames <- gsub("cover\\.", "", sppColumns)
  setnames(inputDataTable, old = sppColumns, new = newCoverColNames)
  message(cli::col_blue("Create initial cohortData object, with no pixelGroups yet"))
  message(cli::col_green("-- Begin reconciling data inconsistencies"))

  imputedPixID <- integer(0)

  inputDataTable[, totalCover := rowSums(.SD), .SDcols = newCoverColNames]
  whEnoughCover <- inputDataTable$totalCover > minCoverThreshold
  message(cli::col_green(
    "  -- Removing all pixels with totalCover <= minCoverThreshold (affects",
    sum(!whEnoughCover),
    "of",
    NROW(inputDataTable),
    "pixels)"
  ))
  message(cli::col_green("     --> resulting in", sum(whEnoughCover), "pixels)"))
  inputDataTable <- inputDataTable[whEnoughCover]

  whAgeEqZero <- which(inputDataTable$age == 0)

  if (!is.null(inputDataTable[["totalBiomass"]])) {
    message(cli::col_green(
      "  -- Setting TotalBiomass in pixel to 0 where age == 0 (affects",
      length(whAgeEqZero),
      "of",
      NROW(inputDataTable),
      "pixels)"
    ))
    message(cli::col_green("     --> keeping ", NROW(inputDataTable), "pixels)"))
    ## correct B in a separate column to keep track of imputed pixels, then replace column
    inputDataTable[, `:=`(totalBiomass2 = totalBiomass)]
    inputDataTable[whAgeEqZero, `:=`(totalBiomass2 = 0)]
    imputedPixID <- inputDataTable[totalBiomass != totalBiomass2, pixelIndex]
    inputDataTable[, totalBiomass := totalBiomass2]
    inputDataTable[, totalBiomass2 := NULL]

    whTotalBEqZero <- which(inputDataTable$totalBiomass == 0)
    message(cli::col_green(
      "  -- Setting age in pixel to 0 where totalBiomass == 0 (affects",
      length(whTotalBEqZero),
      "of",
      NROW(inputDataTable),
      "pixels)"
    ))
    message(cli::col_green("     --> keeping ", NROW(inputDataTable), "pixels)"))
    ## correct age in a separate column to keep track of imputed pixels, then replace column
    inputDataTable[, `:=`(age2 = age)]
    inputDataTable[whTotalBEqZero, `:=`(age2 = 0)]
    imputedPixID <- c(imputedPixID, inputDataTable[age != age2, pixelIndex])
    inputDataTable[, age := age2]
    inputDataTable[, age2 := NULL]
  }

  cohortData <- data.table::melt(
    inputDataTable,
    value.name = "cover",
    measure.vars = newCoverColNames,
    variable.name = "speciesCode"
  )

  ## Remove all cover <= minCoverThreshold
  whCoverGTMinCover <- which(cohortData$cover > minCoverThreshold)
  message(cli::col_green(
    "  -- Removing all cohorts with cover <= minCoverThreshold (affects",
    NROW(cohortData) - length(whCoverGTMinCover),
    "of",
    NROW(cohortData),
    "cohorts"
  ))
  message(cli::col_green("     --> resulting in", length(whCoverGTMinCover), "cohorts)"))
  cohortData <- cohortData[whCoverGTMinCover]
  message(cli::col_green("     --> resulting in", length(unique(cohortData$pixelIndex)), "pixels)"))

  cohortData[, coverOrig := cover]
  if (isTRUE(doAssertion)) {
    if (any(duplicated(cohortData))) {
      warning(".createCohortData: cohortData contains duplicate rows.")
    }
  }

  # if (doAssertion)
  # describeCohortData(cohortData)
  # message(cli::col_green("  -- Assign B = 0 and age = 0 for pixels where cover = 0,\n",
  #             "because cover is most reliable dataset"))

  # hasCover0 <- which(cohortData[["cover"]] == 0)
  cncd <- colnames(cohortData)
  # if (any(c("age", "logAge") %in% cncd)) {
  #   set(cohortData, hasCover0, "age", 0L)
  #   set(cohortData, hasCover0, "logAge", .logFloor(0))
  # }
  if (any(c("B", "totalBiomass") %in% cncd)) {
    # set(cohortData, hasCover0, "B", 0)
    # message(cli::col_green("  -- Assign totalBiomass = 0 if sum(cover) = 0 in a pixel, ",
    #             "  because cover is most reliable dataset"))
    # cohortData <- cohortData[, sum(cover) == 0, by = "pixelIndex"][V1 == TRUE][
    #  cohortData, on = "pixelIndex"][V1 == TRUE, totalBiomass := 0L]
    # cohortData[, V1 := NULL]
  }

  ## CRAZY TODO: DIVIDE THE COVER BY 2 for DECIDUOUS -- will only affect mixed stands
  # message(cli::col_cli::col_green(paste("POSSIBLE ALERT:",
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
      ## Biomass -- by cohort (NOTE: divide by 100 because cover is percent)
      # set(cohortData, NULL, "B", as.numeric(cohortData[["B"]]))
      set(cohortData, NULL, "B", cohortData[["totalBiomass"]] * cohortData[["cover"]] / 100)
      message(cli::col_green("  -- Divide total B of each pixel by the relative cover of the cohorts"))

      # cohortData[ , B := mean(totalBiomass) * cover / 100, by = "pixelIndex"]
      # message(cli::col_blue("Round B to nearest P(sim)$pixelGroupBiomassClass"))
      # cohortData[ , B := ceiling(B / pixelGroupBiomassClass) * pixelGroupBiomassClass]

      message(cli::col_green("Set B to 0 where cover > 0 and age = 0, because B is least quality dataset"))
      cohortData[cover > 0 & age == 0, B := 0L]
      cohortData[, totalBiomass := asInteger(totalBiomass)]
      set(cohortData, NULL, "B", asInteger(cohortData[["B"]]))
    }
  }

  ## clean up
  set(cohortData, NULL, c("totalCover", "coverOrig"), NULL)
  setattr(cohortData, "imputedPixID", imputedPixID)
  return(cohortData)
}

#' Generate initial `cohortData` table
#'
#' Takes a single `data.table` input, which has the following columns in addition to
#' others that will be labelled with species name, and contain percent cover of each:
#'
#' \itemize{
#'   \item `pixelIndex` (integer)
#'   \item `age` (integer)
#'   \item `logAge` (numeric)
#'   \item `initialEcoregionCode` (factor)
#'   \item `totalBiomass` (integer)
#'   \item `lcc` (integer)
#'   \item `rasterToMatch` (integer)
#'   \item `speciesCode` (factor)
#'   \item `cover` (integer)
#'   \item `coverOrig` (integer)
#'   \item `B` (integer)
#' }
#'
#' Several data correction/imputation operations are also performed. Namely, age is imputed
#' in pixels where age data is missing (but not cover) and where `cover == 0` but `age > 0`,
#' total biomass is zeroed if `age == 0`, and age is zeroed if `biomass == 0`.
#'
#' @param inputDataTable A `data.table` with columns described above.
#'
#' @param sppColumns A vector of the names of the columns in `inputDataTable` that
#'   represent percent cover by species, rescaled to sum up to 100%%.
#'
#' @param imputeBadAgeModel statistical model used to impute ages in pixels with missing
#'  data or with cover == 0. If set to NULL no imputation will be attempted, and pixels with
#'  missing age are excluded.
#'
#' @param minCoverThreshold minimum total cover percentage necessary to consider the pixel
#'  vegetated, or a cohort present in a pixel.
#'
#' @template doAssertion
#'
#' @param doSubset Turns on/off subsetting. Defaults to `TRUE`.
#'
#' @return a `cohortData` `data.table` with attribute `"imputedPixID"`
#'     (a vector of pixel IDs that suffered imputation).
#'
#' @author Eliot McIntire
#' @export
makeAndCleanInitialCohortData <- function(
  inputDataTable,
  sppColumns,
  imputeBadAgeModel = quote(lme4::lmer(
    age ~ B * speciesCode + cover * speciesCode + (1 | initialEcoregionCode)
  )),
  minCoverThreshold,
  doAssertion = getOption("LandR.assertions", TRUE),
  doSubset = TRUE
) {
  ## Create groupings
  if (doAssertion) {
    expectedColNames <- c(
      "age",
      "logAge",
      "initialEcoregionCode",
      "totalBiomass",
      "lcc",
      "pixelIndex"
    )
    if (!all(expectedColNames %in% colnames(inputDataTable))) {
      stop("Column names for inputDataTable must include ", paste(expectedColNames, collapse = " "))
    }
    if (!all(sppColumns %in% colnames(inputDataTable))) {
      stop("Species names are incorrect")
    }
    if (
      !all(unlist(lapply(inputDataTable[, sppColumns, with = FALSE], function(x) {
        all(x >= 0 & x <= 100)
      })))
    ) {
      stop(
        "Species columns are not percent cover between 0 and 100. ",
        "This may be because there are more NA values than the Land Cover raster."
      )
    }
  }

  ## .createCohortData deals with the following inconsistencies:
  ## 1. Pixels with total cover (across all species) < 5% are removed;
  ## 2. Pixels with 0 stand age, are assigned 0 stand biomass;
  ## 3. Pixels with 0 stand biomass, are assigned 0 stand age.
  ## 4. (after calcualting cohort B and age): if `cover > 0` and `age == 0`, `B` is set to 0
  cohortData <- Cache(
    .createCohortData,
    inputDataTable = inputDataTable,
    # pixelGroupBiomassClass = pixelGroupBiomassClass,
    minCoverThreshold = minCoverThreshold,
    sppColumns = sppColumns,
    doAssertion = doAssertion
  )
  assertCohortDataAttr(cohortData, doAssertion = doAssertion)
  imputedPixID <- attr(cohortData, "imputedPixID")

  ## Impute missing ages on poor age dataset ----------------------------------------
  ## Cases:
  ##  All species cover = 0 yet totalB > 0
  ## (see other age inconsistencies solved above)

  cohortDataMissingAge <- cohortData[,
    hasBadAge := (age > 0 & cover == 0) | is.na(age) # (age == 0 & cover > 0)#| # ok because cover can be >0 with biomass = 0 #|
    # (B > 0 & age == 0) |
    # (B == 0 & age > 0)
  ][hasBadAge == TRUE] # , by = "pixelIndex"]

  if (NROW(cohortDataMissingAge) > 0) {
    if (!is.null(imputeBadAgeModel)) {
      cohortDataMissingAgeUnique <- unique(
        cohortDataMissingAge,
        by = c("initialEcoregionCode", "speciesCode")
      )[, .(initialEcoregionCode, speciesCode)]
      cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[
        cohortData,
        on = c("initialEcoregionCode", "speciesCode"),
        nomatch = 0
      ]
      cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[
        !is.na(cohortDataMissingAgeUnique$age)
      ]
      cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[, .(
        totalBiomass,
        age,
        speciesCode,
        initialEcoregionCode,
        cover
      )]
      zeros <- sapply(cohortDataMissingAgeUnique, function(x) sum(x == 0))
      if (sum(zeros, na.rm = TRUE)) {
        hasZeros <- zeros[zeros > 0]
        message(
          " ",
          paste(names(hasZeros), collapse = ", "),
          " had ",
          paste(hasZeros, collapse = ", "),
          " zeros, respectively"
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
      cohortDataMissingAgeUnique <- subsetDT(
        cohortDataMissingAgeUnique,
        by = c("initialEcoregionCode", "speciesCode"),
        doSubset = doSubset
      )
      message(cli::col_blue("Impute missing age values: started", Sys.time()))

      outAge <- Cache(
        statsModel,
        modelFn = imputeBadAgeModel,
        uniqueEcoregionGroups = .sortDotsUnderscoreFirst(as.character(unique(
          cohortDataMissingAgeUnique$initialEcoregionCode
        ))),
        .specialData = cohortDataMissingAgeUnique,
        omitArgs = ".specialData"
      )
      message(cli::col_blue("                           completed", Sys.time()))

      # paste with capture.output keeps table structure intact
      messageDF(outAge$rsq, 3, "blue")

      ## A species in cohortDataMissingAge (needs an age imputed) can be absent from the fit set
      ## cohortDataMissingAgeUnique when all its known-age cohorts were dropped just above for zero
      ## biomass and/or cover. The age model then has no coefficient for that fixed-effect speciesCode,
      ## so predict.merMod() below errors with the cryptic "non-conformable arguments" (allow.new.levels
      ## = TRUE only forgives new RANDOM-effect levels, not fixed-effect speciesCode). Name the species.
      ## TODO: handle gracefully rather than erroring -- e.g. enlarge studyArea_biomassParam (only helps
      ## if the species is under-sampled, not structurally cover==0 in the data); or restrict the predict
      ## newdata to species present in the fit set and fallback-impute the rest (e.g. ecoregion mean age);
      ## or make imputeBadAgeModel tolerant of unseen fixed-effect speciesCode levels.
      droppedSpecies <- setdiff(unique(as.character(cohortDataMissingAge$speciesCode)),
                                unique(as.character(cohortDataMissingAgeUnique$speciesCode)))
      if (length(droppedSpecies) > 0L) {
        stop("Cannot impute missing cohort ages: species ", paste(shQuote(droppedSpecies), collapse = ", "),
             " have no usable rows to fit the age model (all their known-age cohorts had zero biomass ",
             "and/or cover), so predict() below would fail with 'non-conformable arguments'. ",
             "See the TODO above for fixes.", call. = FALSE)
      }

      ## allow.new.levels = TRUE because some groups will have only NA for age for all species
      cohortDataMissingAge[,
        imputedAge := pmax(
          0L,
          asInteger(predict(outAge$mod, newdata = cohortDataMissingAge, allow.new.levels = TRUE))
        )
      ]

      cohortData <- cohortDataMissingAge[, .(pixelIndex, imputedAge, speciesCode)][
        cohortData,
        on = c("pixelIndex", "speciesCode")
      ]
      cohortData[!is.na(imputedAge), `:=`(age = imputedAge, logAge = .logFloor(imputedAge))]
      imputedPixID <- c(imputedPixID, unique(cohortData[!is.na(imputedAge), pixelIndex]))
      cohortData[, `:=`(imputedAge = NULL)]
    } else {
      ## if not imputing bad ages, then exclude bad data entries.
      cohortData <- cohortData[
        !cohortDataMissingAge[, .(pixelIndex, speciesCode)],
        on = c("pixelIndex", "speciesCode")
      ]
    }
  }
  ## Round ages to nearest pixelGroupAgeClass
  # set(cohortData, NULL, "age", asInteger(cohortData$age / pixelGroupAgeClass) *
  #       as.integer(pixelGroupAgeClass))

  cohortData[, `:=`(hasBadAge = NULL)]

  ## set B to zero if age is zero because B is lowest quality dataset ---------------
  # message(cli::col_blue("Set recalculate totalBiomass as sum(B);",
  #              "many biomasses will have been set to 0 in previous steps"))
  # cohortData[cover > 0 & age == 0, B := 0L]
  # cohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]

  # This was unused, but for beta regression, this can allow 0s and 1s without
  #   needing a separate model for zeros and ones
  # https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
  # cohortData[ , coverProp := (cover/100 * (NROW(cohortData) - 1) + 0.5) / NROW(cohortData)]

  setattr(cohortData, "imputedPixID", imputedPixID)
  return(cohortData)
}

#' Subset a `data.table` with random subsampling within `by` groups
#'
#' @param DT A `data.table`
#' @param by Character vector of column names to use for groups
#' @param doSubset Logical or numeric indicating the number of subsamples to use
#' @param indices Logical. If `TRUE`, this will return vector of row indices only. Defaults
#'   to `FALSE`, i.e., return the subsampled `data.table`
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
        "using maximum of ",
        sam,
        " samples per combination of ecoregionGroup and speciesCode. ",
        "Change 'doSubset' to a different number if this is not enough"
      )
      ## subset -- add line numbers of those that were sampled
      a <- DT[, list(lineNum = .I[sample(.N, size = min(.N, sam))]), by = by]
      ## select only those row numbers from whole dataset
      if (isFALSE(indices)) {
        DT <- DT[a$lineNum]
      } else {
        DT <- a$lineNum
      }
    } else {
      if (isTRUE(indices)) {
        DT <- seq_len(nrow(DT))
      }
    }
  } else {
    if (isTRUE(indices)) {
      DT <- seq_len(nrow(DT))
    }
  }
  return(DT)
}

#' Drop factor term including interactions from a model formula
#'
#' Based on <https://stackoverflow.com/a/23382097/1380598>.
#'
#' @param form A model formula.
#'
#' @param term Character vector giving the name of the term to drop.
#'
#' @param dropRanEff Logical. If `TRUE` (the default), then the `term` to drop
#'   will also be dropped from the random effects terms.
#'   If `FALSE`, it will only be dropped from the fixed terms.
#'
#' @return An updated model formula.
#'
#' @export
dropTerm <- function(form, term, dropRanEff = TRUE) {
  if (!is(form, "formula")) {
    form <- as.formula(form)
  }

  fterms <- terms(form)
  fac <- attr(fterms, "factors")

  new_form <- form
  termsInner <- rownames(fac)

  for (tt in term) {
    idr <- grepl(tt, termsInner)
    facPartial <- fac[idr, , drop = FALSE]
    toDrop <- list()
    ## Cycle through 1 row at a time of the matrix
    for (rn in seq_len(NROW(facPartial))) {
      ranEff <- grepl("\\|", termsInner[idr][rn])
      if (any(ranEff)) {
        if (isTRUE(dropRanEff)) {
          for (whRE in which(ranEff)) {
            old <- termsInner[idr][rn][whRE]
            if (any(grepl("\\*", old))) {
              stop(
                "This dropTerm function does not work for interaction terms inside the random effects; ",
                "Please rewrite formula or update this source code"
              )
            }
            newRe <- deparse(update(
              Formula(as.formula(paste0("~", old))),
              as.formula(paste0("~ . -", tt))
            ))
            newRe <- gsub("~", "", newRe) ## remove the ~ part to convert to string
            termsInner[idr][rn][whRE] <- newRe
            oldWithParenth <- paste0("(", old, ")") ## random effects must have ( )
            newWithParenth <- paste0("(", newRe, ")") ## random effects must have ( )
            new_form <- update(new_form, paste0(". ~ . - ", oldWithParenth)) ## remove old
            new_form <- update(new_form, paste0(". ~ . + ", newWithParenth)) ## add new
          }
        }
      } else {
        ## Fixed effect terms
        idc <- which(as.logical(facPartial[rn, ]))
        toDrop <- names(facPartial[rn, ][idc])
        needsParenth <- vapply(
          paste0("(", toDrop, ")"),
          FUN = grepl,
          FUN.VALUE = logical(1),
          x = as.character(new_form)[3],
          fixed = TRUE
        )
        if (any(needsParenth)) {
          toDrop[needsParenth] <- paste0("(", toDrop[needsParenth], ")")
        }

        new_form <- update(new_form, paste0(". ~ . -", paste(toDrop, collapse = " - ")))
      }
    }
  }

  return(new_form)
}

#' The generic statistical model to run (`lmer` or `glmer`)
#'
#' This does a few things including R squared, gets the fitted values.
#' It appears that running the models "as is" without this wrapper does not work with `Cache`.
#' The return of the model in a list solves this problem.
#' For Caching, the `.specialData` should be "omitted" via `omitArgs`, and
#' `uniqueEcoregionGroups` should not be omitted.
#'
#' @param modelFn A quoted expression of type `package::model(Y ~ X, ...)`, omitting
#'   the `data` argument. E.g., `lme4::glmer(Y ~ X + (X|G), family = poisson)`.
#'
#' @param uniqueEcoregionGroups Unique values of `ecoregionGroups`.
#'   This is the basis for the statistics, and can be used to optimize caching,
#'   e.g., ignore `.specialData` in `.omitArgs`.
#'
#' @param sumResponse a sum of all the response variable values.
#'   Also to be used to optimize caching, e.g. ignore `.specialData` in `.omitArgs`.
#' @param .specialData The custom dataset required for the model.
#'
#' @export
statsModel <- function(modelFn, uniqueEcoregionGroups, sumResponse, .specialData) {
  .requireNamespace("MuMIn", stopOnFALSE = TRUE)

  ## convert model call to vector of arguments
  modelArgs <- as.character(modelFn)
  names(modelArgs) <- names(modelFn)

  ## function name
  fun <- modelArgs[[1]]

  ## get formula and check
  form <- tryCatch(
    as.formula(modelArgs[2], env = .GlobalEnv), # .GlobalEnv keeps object small
    error = function(e) {
      stop(paste(
        "Could not convert '",
        modelArgs[2],
        "'to formula.",
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
      form <- as.formula(form, env = .GlobalEnv)

      ## change to glm if dropping random effects
      fun <- "stats::glm"

      whereFam <- grep("family", names(modelArgs))
      modelFn2 <- if (length(whereFam)) {
        call(fun, form, family = modelArgs[[whereFam]])
      } else {
        call(fun, form)
      }

      message(
        cli::col_blue("Grouping variable "),
        cli::col_red("only has one level"),
        cli::col_blue(". "),
        cli::col_blue(
          "Formula changed to\n",
          cli::col_magenta(paste0(format(modelFn2, appendLF = FALSE), collapse = ""))
        )
      )
    }
  }

  ## get function and check
  fun <- .extractFunction(fun)
  if (!is.function(fun)) {
    stop(paste0(
      "Can't find the function '",
      modelArgs[1],
      "'.",
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

  ## drop factor terms with a single level
  singles <- names(which(sapply(lapply(.specialData, unique), length) == 1))
  keep <- which(unname(vapply(
    singles,
    grepl,
    x = paste(modelArgs$formula, collapse = " "),
    logical(1)
  )))
  singles <- singles[keep]
  if (length(singles) > 0) {
    modelArgs$formula <- dropTerm(modelArgs$formula, singles)
  }

  mod <- do.call(fun, modelArgs)

  list(mod = mod, pred = fitted(mod), rsq = MuMIn::r.squaredGLMM(mod))
}

#' Default columns that define cohorts
#'
#' @note because the name `cohortDefinitionCols` is also used as a function argument,
#' be sure to use `LandR::cohortDefinitionCols()` in those functions or you'll get a
#' "promise already under evaluation" error.
#'
#' @export
cohortDefinitionCols <- function() {
  ## 2024-08-08: do not include ecoregionGroup and B when defining cohorts:
  ## - ecoregionGroup already taken into account with pixelGroup;
  ## - species with same age *should* have the same B already;
  c("pixelGroup", "speciesCode", "age")
}

#' Default columns that define pixel groups
#'
#' @note because the name `columnsForPixelGroups` is also used as a function argument,
#' be sure to use `LandR::columnsForPixelGroups()` in those functions or you'll get a
#' "promise already under evaluation" error.
#'
#' @export
columnsForPixelGroups <- function() {
  c("ecoregionGroup", "speciesCode", "age", "B")
}

#' Generate `cohortData` table per pixel:
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
#' An expanded `cohortData` `data.table` with a new `pixelIndex` column.
#'
#' @export
addPixels2CohortData <- function(
  cohortData,
  pixelGroupMap,
  cohortDefinitionCols = LandR::cohortDefinitionCols(),
  doAssertion = getOption("LandR.assertions", TRUE)
) {
  assertCohortData(
    cohortData,
    pixelGroupMap,
    cohortDefinitionCols = LandR::cohortDefinitionCols(),
    doAssertion = doAssertion
  )

  pixelGroupTable <- na.omit(data.table(
    pixelGroup = as.vector(pixelGroupMap[]),
    pixelIndex = 1:ncell(pixelGroupMap)
  ))
  pixelCohortData <- cohortData[
    pixelGroupTable,
    on = "pixelGroup",
    nomatch = 0,
    allow.cartesian = TRUE
  ]

  assertPixelCohortData(pixelCohortData, pixelGroupMap, doAssertion = doAssertion)

  return(pixelCohortData)
}

#' Add number of pixels per `pixelGroup` and add it has a new column to `cohortData`
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template cohortDefinitionCols
#'
#' @template doAssertion
#'
#'
#' @return
#' An `cohortData` `dat.table` with a new `noPixels`
#' column
#'
#' @export
addNoPixel2CohortData <- function(
  cohortData,
  pixelGroupMap,
  cohortDefinitionCols = LandR::cohortDefinitionCols(),
  doAssertion = getOption("LandR.assertions", TRUE)
) {
  assertCohortData(
    cohortData,
    pixelGroupMap,
    cohortDefinitionCols = LandR::cohortDefinitionCols(),
    doAssertion = doAssertion
  )

  noPixelsXGroup <- data.table(
    noPixels = tabulate(pixelGroupMap[]),
    pixelGroup = c(1:maxFn(pixelGroupMap))
  )

  pixelCohortData <- cohortData[noPixelsXGroup, on = "pixelGroup", nomatch = 0]

  if (doAssertion) {
    test1 <- length(setdiff(pixelCohortData$pixelGroup, cohortData$pixelGroup)) > 0
    test2 <- sum(unique(pixelCohortData[, .(pixelGroup, noPixels)])$noPixels) !=
      sum(!is.na(pixelGroupMap[]) & pixelGroupMap[] != 0) ## 0's have no cohorts.

    if (test1 || test2) {
      stop("pixelGroups differ between pixelCohortData/pixelGroupMap and cohortData")
    }
  }

  return(pixelCohortData)
}

#' Make the `cohortData` table, while modifying the temporary
#' `pixelCohortData` that will be used to prepare other files.
#'
#' Takes a `pixelCohortData` table (see `makeAndCleanInitialCohortData`),
#'   the `speciesEcoregion` list and returns a modified `pixelCohortData` and
#'   the `cohortData` tables to be used in the simulation.
#'   This function mainly removes unnecessary columns from `pixelCohortData`,
#'   subsets pixels with `biomass > 0`, generates `pixelGroups`,
#'   and adds `ecoregionGroup` and `totalBiomass` columns to `pixelCohortData`.
#'   `cohortData` is then created by subsetting unique combinations of `pixelGroup` and
#'   whatever columns are listed in `columnsForPixelGroups`.
#'   The resulting `cohortData` table has the following columns:
#' \itemize{
#'   \item `speciesCode` (factor)
#'   \item `ecoregionGroup` (factor)
#'   \item `pixelGroup` (integer)
#'   \item `age` (integer)
#'   \item `B` (integer)
#' }
#'
#' @template pixelCohortData
#' @param columnsForPixelGroups Default columns that define pixel groups
#' @param pixelGroupAgeClass Integer. When assigning `pixelGroup` membership, this defines the
#'   resolution of ages that will be considered 'the same `pixelGroup`', e.g., if it is 10,
#'   then 6 and 14 will be the same.
#' @param pixelGroupBiomassClass Integer. When assigning `pixelGroup` membership, this defines
#'   the resolution of biomass that will be considered 'the same `pixelGroup`', e.g., if it is
#'   100, then 5160 and 5240 will be the same
#' @param minAgeForGrouping Minimum age for regrouping. This may be because there is a source of
#'   ages for young stands/trees that is very reliable, such as a fire database.
#'   Ages below this will not be grouped together. Defaults to -1, meaning treat all ages equally.
#'   If this is related to known ages from a high quality database, then use age of the oldest
#'   trees in that database.
#' @param pixelFateDT A `data.table` of `pixelFateDT`; if none provided, will make an empty one.
#' @param rmImputedPix Should imputed pixels be removed
#' @param imputedPixID a vector of IDs of pixels that suffered data imputation.
#'
#' @template speciesEcoregion
#'
#' @return A list with a modified `pixelCohortData`, `cohortData`, and `pixelFateDT`
#' `data.table`s.
#'
#' @export
makeCohortDataFiles <- function(
  pixelCohortData,
  columnsForPixelGroups,
  speciesEcoregion,
  pixelGroupBiomassClass,
  pixelGroupAgeClass,
  minAgeForGrouping = 0,
  rmImputedPix = FALSE,
  imputedPixID,
  pixelFateDT
) {
  ## make ecoregioGroup a factor (again) and remove unnecessary cols.
  ## refactor because the "_34" and "_35" ones are still levels
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))]
  cols <- intersect(
    c("logAge", "coverOrig", "totalBiomass", "initialEcoregionCode", "cover", "lcc"),
    names(pixelCohortData)
  )
  set(pixelCohortData, j = cols, value = NULL)

  ## Round ages to nearest pixelGroupAgeClass
  pixelCohortData[
    age > minAgeForGrouping,
    age := asInteger(age / pixelGroupAgeClass) * as.integer(pixelGroupAgeClass)
  ]

  ## Round Biomass to nearest pixelGroupBiomassClass
  message(cli::col_blue("Round B to nearest P(sim)$pixelGroupBiomassClass"))
  pixelCohortData[
    # age > minAgeForGrouping,
    ,
    B := asInteger(B / pixelGroupBiomassClass) * as.integer(pixelGroupBiomassClass)
  ]

  ## Remove B == 0 cohorts after young removals
  message(cli::col_green(
    "  -- Removing cohorts with B = 0 and age > 0 -- these were likely poor predictions from updateYoungBiomasses"
  ))
  whBEqZeroAgeGT0 <- which(pixelCohortData$B == 0 & pixelCohortData$age > 0)

  if (length(whBEqZeroAgeGT0) > 0) {
    pixelCohortData2 <- pixelCohortData[-whBEqZeroAgeGT0]
  } else {
    pixelCohortData2 <- pixelCohortData
  }

  lostPixels <- setdiff(pixelCohortData$pixelIndex, pixelCohortData2$pixelIndex)
  message(cli::col_green(
    "     affected",
    length(whBEqZeroAgeGT0),
    "cohorts, in",
    length(lostPixels),
    "pixels;"
  ))
  lenUniquePix <- length(unique(pixelCohortData2$pixelIndex))
  message(cli::col_green("     leaving", lenUniquePix, "pixels"))
  pixelFateDT <- pixelFate(
    pixelFateDT,
    fate = "rm pixels with Biomass == 0, after updating young cohort B",
    length(lostPixels),
    runningPixelTotal = lenUniquePix
  )

  pixelCohortData <- pixelCohortData2
  ## Set B to 0 if age is 0
  # whAgeZero <- which(pixelCohortData$age == 0)
  # if (length(whAgeZero)) {
  #   message(cli::col_green("    -- There were", length(whAgeZero), "pixels with age = 0; forcing B to zero"))
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

  message(cli::col_blue(
    "Removing some pixels because their species * ecoregionGroup combination has no age or B data to estimate ecoregion traits:"
  ))
  # message(cli::col_blue(paste(sort(unique(pixelCohortData[!ecoregionGroup %in% ecoregionsWeHaveParametersFor][["ecoregionGroup]])), collapse = ", ")))
  cols <- c("speciesCode", "ecoregionGroup")
  messageDF(
    colour = "blue",
    pixelCohortData[!speciesEcoregion, on = cols][, ..cols][,
      list(numPixelsRemoved = .N),
      by = cols
    ], # anti-join
  )

  ## REMOVE PIXELS IN ECOREGION GROUPS THAT ENDED UP WITHOUT PARAMS
  # pixelCohortData <- pixelCohortData[ecoregionGroup %in% ecoregionsWeHaveParametersFor] # keep only ones we have params for
  pixelCohortData <- speciesEcoregion[, ..cols][pixelCohortData, on = cols, nomatch = 0]

  pixelFateDT <- pixelFate(
    pixelFateDT,
    "removing ecoregionGroups without enough data to est. maxBiomass",
    tail(pixelFateDT$runningPixelTotal, 1) - NROW(unique(pixelCohortData$pixelIndex))
  )

  ## REMOVE IMPUTED PIXELS FROM THE SIMULATION IF NEED BE
  if (rmImputedPix) {
    pixelCohortData <- pixelCohortData[!pixelIndex %in% imputedPixID]

    pixelFateDT <- pixelFate(
      pixelFateDT,
      "removing pixels that suffered data imputation",
      tail(pixelFateDT$runningPixelTotal, 1) - NROW(unique(pixelCohortData$pixelIndex))
    )
  }

  ## Lost some ecoregionGroups -- refactor
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))]

  cd <- pixelCohortData[, .SD, .SDcols = c("pixelIndex", columnsForPixelGroups())]
  pixelCohortData[,
    pixelGroup := Cache(
      generatePixelGroups,
      cd,
      maxPixelGroup = 0,
      columns = LandR::columnsForPixelGroups()
    )
  ]

  pixelCohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]

  cohortData <- unique(pixelCohortData, by = c("pixelGroup", columnsForPixelGroups()))
  cohortData[, `:=`(pixelIndex = NULL)]

  assertUniqueCohortData(cohortData, c("pixelGroup", "ecoregionGroup", "speciesCode"))
  return(list(
    cohortData = cohortData,
    pixelCohortData = pixelCohortData,
    pixelFateDT = pixelFateDT
  ))
}

#' Create new cohorts based on provenance table with unique `pixelGroup` and add to `cohortData`
#'
#' @param newPixelCohortData the cohorts that were harvested
#' @template cohortData
#' @template pixelGroupMap
#' @template currentTime
#' @template successionTimestep
#' @param trackPlanting adds column that tracks planted cohorts if `TRUE`
#' @param initialB the initial biomass of new cohorts. Defaults to ten.
#'
#' @return A `data.table` with a new `cohortData`
#'
#' @export
plantNewCohorts <- function(
  newPixelCohortData,
  cohortData,
  pixelGroupMap,
  initialB = 10,
  currentTime,
  successionTimestep,
  trackPlanting = FALSE
) {
  ## get spp "productivity traits" per ecoregion/present year

  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, pixelGroup := as.vector(pixelGroupMap[])[pixelIndex]]
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

  ## Ceres: temporary workaround to ensure the default value is used even
  ## if LANDIS behaviour is being used to calculate initial B in non-planted cohorts
  if (isTRUE(is.na(initialB)) || is.null(initialB)) {
    initialB <- formals(plantNewCohorts)$initialB
  }
  newCohortData[, B := initialB]

  # Here we subset cohortData instead of setting added columns to NULL. However, as these are 'new' cohorts, this is okay
  newCohortData <- newCohortData[, .(
    pixelGroup,
    ecoregionGroup,
    speciesCode,
    age,
    B,
    Provenance,
    mortality = 0L,
    aNPPAct = 0L
  )]

  if (getOption("LandR.assertions")) {
    if (
      isTRUE(
        NROW(unique(newCohortData, by = c("pixelGroup", "age", "speciesCode", "Provenance"))) !=
          NROW(newCohortData)
      )
    ) {
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

#' Add cohorts to `cohortData` and `pixelGroupMap`
#'
#' This is a wrapper for  `generatePixelGroups`, `initiateNewCohort` and updates to
#' `pixelGroupMap` via assignment to new `pixelIndex` values in `newPixelCohortData`.
#' By running these all together, there is less chance that they will diverge.
#' There are some checks internally for consistency.
#'
#' Does the following:
#' \enumerate{
#'   \item add new cohort data into `cohortData`;
#'   \item assign initial `B` and `age` for new cohort;
#'   \item assign the new `pixelGroup` to the pixels that have new cohort;
#'   \item update the `pixelGroup` map.
#' }
#'
#' @template newPixelCohortData
#' @template cohortData
#' @template pixelGroupMap
#' @template cohortDefinitionCols
#' @template currentTime
#' @template speciesEcoregion
#'
#' @param treedHarvestPixelTable A data.table with at least 2 columns, `pixelIndex` and `pixelGroup`.
#'   This will be used in conjunction with `cohortData` and `pixelGroupMap`
#'   to ensure that everything matches correctly.
#' @template successionTimestep
#' @param provenanceTable A `data.table` with three columns:
#' New cohorts are initiated at the `ecoregionGroup` `speciesEcoregion` from the
#' corresponding `speciesEcoregion` listed in the `Provenance` column
#' @param initialB the initial biomass of new cohorts. Defaults to ten, even if NA/NULL is passed.
#' @param trackPlanting if true, planted cohorts in `cohortData` are tracked with `TRUE`
#' in column 'planted'
#'
#' @template verbose
#'
#' @template doAssertion
#'
#' @return
#' A list of length 2, `cohortData` and `pixelGroupMap`, with
#' `newPixelCohortData` inserted.
#'
#' @export
#' @rdname updateCohortDataPostHarvest
updateCohortDataPostHarvest <- function(
  newPixelCohortData,
  cohortData,
  pixelGroupMap,
  currentTime,
  speciesEcoregion,
  treedHarvestPixelTable = NULL,
  successionTimestep,
  provenanceTable,
  trackPlanting = FALSE,
  initialB = 10,
  cohortDefinitionCols = LandR::cohortDefinitionCols(),
  verbose = getOption("LandR.verbose", TRUE),
  doAssertion = getOption("LandR.assertions", TRUE)
) {
  cohortData <- copy(cohortData)
  provenanceTable <- copy(provenanceTable)

  maxPixelGroup <- as.integer(maxFn(pixelGroupMap))

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
    speciesCode,
    ecoregionGroup
  )

  specieseco_current[, maxB_eco := max(maxB), by = ecoregionGroup]

  setkey(newPixelCohortData, speciesCode, ecoregionGroup)

  newPixelCohortData <- specieseco_current[newPixelCohortData]

  columnsForPG <- c(LandR::columnsForPixelGroups(), "maxB", "maxANPP", "Provenance")

  cd <- newPixelCohortData[, c("pixelIndex", columnsForPG), with = FALSE]

  newPixelCohortData[,
    pixelGroup := generatePixelGroups(cd, maxPixelGroup = maxPixelGroup, columns = columnsForPG)
  ]

  # Remove the duplicated pixels within pixelGroup (i.e., 2+ species in the same pixel)
  pixelsToChange <- unique(
    newPixelCohortData[, c("pixelIndex", "pixelGroup")],
    by = c("pixelIndex")
  )

  pixelGroupMap[pixelsToChange$pixelIndex] <- pixelsToChange$pixelGroup

  if (doAssertion) {
    if (!isTRUE(all(pixelsToChange$pixelGroup == pixelGroupMap[][pixelsToChange$pixelIndex]))) {
      stop("pixelGroupMap and newPixelCohortData$pixelGroupMap don't match in updateCohortData fn")
    }
  }

  ## Add new cohorts and rm missing cohorts (i.e., those pixelGroups that are gone) -----------
  cohortData <- plantNewCohorts(
    newPixelCohortData,
    cohortData,
    pixelGroupMap,
    currentTime = currentTime,
    successionTimestep = successionTimestep,
    initialB = initialB,
    trackPlanting = trackPlanting
  )

  outs <- rmMissingCohorts(
    cohortData,
    pixelGroupMap,
    cohortDefinitionCols = LandR::cohortDefinitionCols()
  )

  assertCohortData(
    outs$cohortData,
    outs$pixelGroupMap,
    cohortDefinitionCols = LandR::cohortDefinitionCols(),
    doAssertion = doAssertion,
    verbose = verbose
  )

  if (verbose > 0) {
    nPixForest <- sum(!is.na(outs$pixelGroupMap[]))
    nPixGrps <- length(unique(outs$cohortData$pixelGroup))
    nPixNoPixGrp <- sum(outs$pixelGroupMap[] == 0, na.rm = TRUE)
    nPixTreed <- sum(outs$pixelGroupMap[] != 0, na.rm = TRUE)

    nDigits <- max(nchar(c(nPixForest, nPixGrps, nPixNoPixGrp))) + 3
    message(cli::col_magenta(
      "NUMBER OF FORESTED PIXELS          :",
      paddedFloatToChar(nPixForest, padL = nDigits, pad = " ")
    ))
    message(cli::col_magenta(
      "NUMBER OF PIXELS WITH TREES        :",
      paddedFloatToChar(nPixTreed, padL = nDigits, pad = " ")
    ))
    message(cli::col_magenta(
      "NUMBER OF UNIQUE PIXELGROUPS       :",
      paddedFloatToChar(nPixGrps, padL = nDigits, pad = " ")
    ))
    message(cli::col_magenta(
      "NUMBER OF PIXELS WITH NO PIXELGROUP:",
      paddedFloatToChar(nPixNoPixGrp, padL = nDigits, pad = " ")
    ))
  }

  return(list(cohortData = outs$cohortData, pixelGroupMap = outs$pixelGroupMap))
}

#' Create or amend data to a `pixelFateDT` object
#'
#' @param pixelFateDT A `pixelFateDT` `data.table` with 3 columns: `fate`,
#'   `pixelsRemoted`, and `runningPixelTotal`.
#' @param fate A character string (length 1) describing in words the change
#' @param pixelsRemoved A numeric indicating how many pixels were removed due to the `fate`.
#' @param runningPixelTotal an optional numeric with new, running total. If not supplied,
#'   it will be calculated from the last row of `pixelFateDT` `runningTotal` minus the
#'   `pixelsRemoved`
#'
#' @return A `pixelFateDT` object, updated with one extra row.
#'
#' @export
pixelFate <- function(
  pixelFateDT,
  fate = NA_character_,
  pixelsRemoved = 0,
  runningPixelTotal = NA_integer_
) {
  if (missing(pixelFateDT)) {
    pixelFateDT <- data.table(
      fate = character(),
      pixelsRemoved = integer(),
      runningPixelTotal = integer()
    )
  }
  if (is.na(runningPixelTotal)) {
    runningPixelTotal <- tail(pixelFateDT$runningPixelTotal, 1) - pixelsRemoved
  }
  pixelFateDT <- rbindlist(list(
    pixelFateDT,
    data.table(fate = fate, pixelsRemoved = pixelsRemoved, runningPixelTotal = runningPixelTotal)
  ))
  pixelFateDT
}

#' Generate and add vegetation type column to `cohortData`
#'
#' This function is a simplification of `vegTypeMapGenerator` and instead of generating a map,
#' it adds the vegetation type column to the `cohortData` table.
#'
#' @param x A `cohortData` object
#'
#' @template vegLeadingProportion
#'
#' @param mixedType An integer defining whether mixed stands are:
#'                  not differentiated (0), any kind of species admixture (1),
#'                  or deciduous mixed with conifer (2; default).
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param pixelGroupColName Name of the column in `pixelGroup` to use.
#'
#' @template doAssertion
#'
#' @param ... Additional arguments.
#'
#' @return `x` with a new column, 'leading', coding the vegetation type
#'    of each group defined by `pixelGroupColName`
#'
#' @author Eliot McIntire, Ceres Barros, Alex Chubaty
#' @export
#' @rdname vegTypeGenerator
#'
#' @examples
#' library(data.table)
#' x <- data.table(
#'   pixelGroup = rep(1:2, each = 2), B = c(100, 200, 20, 400),
#'   speciesCode = rep(c("Pice_Gla", "Popu_Tre"), 2)
#' )
#' vegTypeGenerator(x)
vegTypeGenerator <- function(
  x,
  vegLeadingProportion = 0.8,
  mixedType = 2,
  sppEquiv = NULL,
  sppEquivCol,
  pixelGroupColName = "pixelGroup",
  doAssertion = getOption("LandR.assertions", TRUE),
  ...
) {
  stopifnot(mixedType %in% 0:2, length(mixedType) == 1)

  nrowCohortData <- NROW(x)

  leadingBasedOn <- preambleVTG(x, vegLeadingProportion, doAssertion, nrowCohortData)

  if (mixedType == 2) {
    if (is.null(sppEquiv)) {
      sppEquiv <- get(
        data("sppEquivalencies_CA", package = "LandR", envir = environment()),
        inherits = FALSE
      )

      # Find the sppEquivCol that best matches what you have in x
      sppEquivCol <- names(sort(
        sapply(sppEquiv, function(xx) sum(xx %in% unique(x$species))),
        decreasing = TRUE
      )[1])
      message(paste0(
        "Using mixedType == 2, but no sppEquiv provided. ",
        "Attempting to use data('sppEquivalencies_CA', 'LandR') ",
        "and sppEquivCol = '",
        sppEquivCol,
        "'"
      ))
    }
  }

  ## use new vs old algorithm based on size of x. new one (2) is faster in most cases.
  ## enable assertions to view timings for each algorithm before deciding which to use.
  ## 2025-01: old algo faster for larger tables and should be used by default (see #39)

  pgdAndSc <- c(pixelGroupColName, "speciesCode")
  pgdAndScAndLeading <- c(pgdAndSc, leadingBasedOn)
  totalOfLeadingBasedOn <- paste0("total", leadingBasedOn)
  speciesOfLeadingBasedOn <- paste0("speciesGroup", leadingBasedOn)

  ## 1. Find length of each pixelGroup -- don't include pixelGroups in "by" that have only 1 cohort: N = 1
  cohortData1 <- copy(x)
  systimePre1 <- Sys.time()
  pixelGroupData1 <- cohortData1[, list(N = .N), by = pixelGroupColName]

  ## Calculate speciesProportion from cover or B
  pixelGroupData1 <- cohortData1[, ..pgdAndScAndLeading][pixelGroupData1, on = pixelGroupColName]
  set(pixelGroupData1, NULL, totalOfLeadingBasedOn, pixelGroupData1[[leadingBasedOn]])

  if (identical(leadingBasedOn, "cover")) {
    pixelGroupData1[
      N != 1,
      (totalOfLeadingBasedOn) := sum(cover, na.rm = TRUE),
      by = pixelGroupColName
    ]
    pixelGroupData1 <- pixelGroupData1[,
      list(sum(cover, na.rm = TRUE), totalcover[1]),
      by = pgdAndSc
    ]
  } else {
    pixelGroupData1[N != 1, (totalOfLeadingBasedOn) := sum(B, na.rm = TRUE), by = pixelGroupColName]
    pixelGroupData1 <- pixelGroupData1[, list(sum(B, na.rm = TRUE), totalB[1]), by = pgdAndSc]
  }
  setnames(
    pixelGroupData1,
    old = c("V1", "V2"),
    new = c(speciesOfLeadingBasedOn, totalOfLeadingBasedOn)
  )

  set(
    pixelGroupData1,
    NULL,
    "speciesProportion",
    pixelGroupData1[[speciesOfLeadingBasedOn]] / pixelGroupData1[[totalOfLeadingBasedOn]]
  )
  systimePost1 <- Sys.time()

  setorderv(pixelGroupData1, pixelGroupColName)

  pixelGroupData <- pixelGroupData1
  rm(pixelGroupData1)

  ## Determine "mixed" -----------------------------------
  if (mixedType == 0) {
    ## 1. sort on pixelGroup and speciesProportion, reverse so 1st row of each pixelGroup is the largest
    ## 2. Keep only first row in each pixelGroup
    pixelGroupData3 <- pixelGroupData[, list(
      speciesCode,
      get(pixelGroupColName),
      speciesProportion
    )]
    setnames(pixelGroupData3, "V2", pixelGroupColName)
    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName]
    setnames(pixelGroupData3, "speciesCode", "leading")
  } else if (mixedType == 1) {
    ## create "mixed" class (2019-05: faster than previous below)
    ## 1. anything with >= vegLeadingProportion is "pure"
    ## 2. sort on pixelGroup and speciesProportion, reverse so that 1st row of each pixelGroup is the largest
    ## 3. Keep only first row in each pixelGroup
    ## 4. change column names and convert pure to mixed ==> mixed <- !pure
    pixelGroupData3 <- pixelGroupData[, list(
      pure = speciesProportion >= vegLeadingProportion,
      speciesCode,
      get(pixelGroupColName),
      speciesProportion
    )]
    setnames(pixelGroupData3, "V3", pixelGroupColName)
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
    mixedType2Condition <- quote(
      Type == "Deciduous" &
        speciesProportion < vegLeadingProportion &
        speciesProportion > 1 - vegLeadingProportion
    )
    pixelGroupData3[, mixed := FALSE]

    pixelGroupData3[eval(mixedType2Condition), mixed := TRUE, by = pixelGroupColName]
    pixelGroupData3[, mixed := any(mixed), by = pixelGroupColName]

    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    set(pixelGroupData3, NULL, "Type", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName] ## sp. w/ highest prop. per pixelGroup
    pixelGroupData3[mixed == TRUE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    set(pixelGroupData3, NULL, "leading", factor(pixelGroupData3[["leading"]]))
  } else {
    stop("invalid mixedType! Must be one of '0', '1' or '2'.")
  }

  cols <- c(pixelGroupColName, "leading")
  xx <- pixelGroupData3[, ..cols][x, on = pixelGroupColName]

  return(xx)
}

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
  if (!nrowCohortData > 0) {
    stop("cohortData is empty")
  }

  if (isTRUE(doAssertion)) {
    message("LandR::vegTypeMapGenerator: NROW(x) == ", nrowCohortData)
  }

  leadingBasedOn
}

# derived from plyr::mapvalues
mapvalues2 <- function(x, from, to) {
  #
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}

#' Reduces the age of cohorts that exceed their `longevity x adjustmentFactor`
#'
#' @template pixelCohortData
#'
#' @param longevity A data.table with the longevity of each species.
#'
#' @param adjustmentFactor A numeric controlling the proportion of species longevity
#' that cohort ages cannot exceed.
#'
#' @returns A `cohortData` data.table with corrected ages.
#'
#' @export
adjustAgeToLongevity <- function(pixelCohortData, longevity, adjustmentFactor) {
  ## Check inputs requirements
  if (!all(c("longevity", "speciesCode") %in% colnames(longevity))) {
    stop("longevity data.frame needs the columns longevity and speciesCode")
  }
  if (!is.numeric(adjustmentFactor) | adjustmentFactor < 0.5 | adjustmentFactor > 1) {
    stop("adjustmentFactor needs to be a number between 0.5 and 1")
  }

  ## Calculate the maximum age accepted for each species
  maxAges <- longevity[, .(speciesCode, maxAge = round(longevity * adjustmentFactor))]
  ## Correct the age for cohorts that exceed that limit
  correctedPixelCohortData <- pixelCohortData[maxAges, on = .(speciesCode)]
  ## Identify species for which some cohorts exceed longevity*adjustmentFactor
  speciesToCorrect <- unique(as.character(correctedPixelCohortData[age > maxAge, speciesCode]))
  for (sp in speciesToCorrect) {
    message(
      "Adjusting ages of some ",
      sp,
      " cohorts that exceed `longevity*P(sim)$adjustmentFactor`."
    )
    sp_row <- which(correctedPixelCohortData[, "speciesCode"] == sp)
    sp_longevity <- longevity[speciesCode == sp, longevity]
    correctedPixelCohortData[sp_row, "age"] <- ageAdjust(
      age = correctedPixelCohortData[sp_row, age],
      adjustmentFactor = adjustmentFactor,
      longevity = sp_longevity,
      decayRange = 20
    )
  }
  correctedPixelCohortData[, maxAge := NULL]
  return(correctedPixelCohortData)
}

#' Applies the smoothed age correction for a single species.
#'
#' @inheritParams adjustAgeToLongevity
#'
#' @param age input age values
#'
#' @param decayRange range along the x-axis where [decayFunction] is applied
#'
#' @keywords internal
ageAdjust <- function(age, adjustmentFactor, longevity, decayRange = 96) {
  ageOrig <- age
  maxAge <- round(adjustmentFactor * longevity)
  maxUnaffectedAge <- maxAge - decayRange - 1

  ageFromMaxUnaffectedAge <- round(age - maxUnaffectedAge)
  whNeedAdjusting <- which(ageFromMaxUnaffectedAge > 0)

  age3 <- decayFunction(x = ageOrig[whNeedAdjusting], Ymax = maxAge, startAt = maxUnaffectedAge)
  ageOrig[whNeedAdjusting] <- age3

  newAge <- round(ageOrig)
  return(newAge)
}

#' Apply the decay function
#'
#' Exponential function starting at `startAt` and with an asymptote at `Ymax`
#'
#' @param x input value along the x-axis
#' @param Ymax maximum value on the y-axis, defining the asymptote
#' @param startAt value along the x-axis where the decay function is applied
#'
#' @keywords internal
decayFunction <- function(x, Ymax, startAt = 120) {
  Ymax - (Ymax - startAt) * exp(-(x - startAt) / (Ymax - startAt))
}
