utils::globalVariables(c(".N", "V1", "V2", "relativeAbundObsrvd", "relativeAbund",
                         "reps", "years"))

#' Assertions
#'
#' - `assert1`: Assert that `ecoregionCodes` that were replaced, were correctly identified;
#'
#' @param cohortData34to36 A `cohortData` `data.table` with only the
#'                         pixels what were LCC 34:36
#' @param rmZeroBiomassQuote  An expression to evaluate, in the form of `quote(B>0)`,
#'    used to select cohorts with biomass.
#'
#' @param classesToReplace Integer vector of classes that are are to be replaced,
#'     e.g., 34, 35, 36 on LCC2005, which are burned young, burned 10 year, and cities.
#'
#' @template cohortData
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assert1 <- function(cohortData34to36, cohortData, rmZeroBiomassQuote, classesToReplace = 34:36,
                    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    allCodesAre34to36 <- all(grepl(paste(paste0(".*_", classesToReplace), collapse = "|"),
                                   as.character(cohortData34to36$initialEcoregionCode)))
    if (!allCodesAre34to36)
      stop("lcc classes were mismanaged; contact developers: code 234")

    ## do they match the codes found for pre selected areas (e.g. areas with biomass)?
    ## if not, is this because the non-matching codes had no biomass to start with?
    ## or is it because they need to be masked?
    if (!is.null(rmZeroBiomassQuote)) {
      onlyExistingCodes <- all(unique(cohortData34to36$ecoregionGroup) %in%
                                 unique(cohortData[eval(rmZeroBiomassQuote), initialEcoregionCode]))
    } else {
      onlyExistingCodes <- all(unique(cohortData34to36$ecoregionGroup) %in%
                                 unique(cohortData$initialEcoregionCode))
    }

    if (!onlyExistingCodes) {
      ## was there any biomass/species data in these pixels?
      if (!is.null(rmZeroBiomassQuote)) {
        temp1 <- unique(cohortData34to36$ecoregionGroup)
        temp2 <- unique(cohortData[eval(rmZeroBiomassQuote), initialEcoregionCode])
        nonMatching <- setdiff(temp1, temp2)

        nonMatchingB <- sum(cohortData[initialEcoregionCode %in% nonMatching, B], na.rm = TRUE)
        if (nonMatchingB)
          stop("There are some ecoregionCodes created post replacement of 34 and 35")
      }

      ## are they pixels that will be masked because they couldn't be converted?
      if (!is.null(rmZeroBiomassQuote)) {
        temp1 <- unique(cohortData34to36$ecoregionGroup)
        temp2 <- unique(cohortData[, initialEcoregionCode])
        nonMatching <- setdiff(temp1, temp2)
        if (!all(is.na(nonMatching)))
          stop("There are some ecoregionCodes created post replacement of 34 and 35")
      }
    }
  }
}

#' - `assert2`: assert that `ecoregionCodes` that were replaced, were correctly identified;
#'
#' @param cohortDataNo34to36 A `cohortData` `data.table` with only the
#'                         pixels what were LCC 34:36
#'
#' @param classesToReplace Integer vector of classes that are are to be replaced,
#'     e.g., 34, 35, 36 on LCC2005, which are burned young, burned 10 year, and cities.
#'
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assert2 <- function(cohortDataNo34to36, classesToReplace = 34:36,
                    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    noCodesAre34to36 <- all(!grepl(paste(paste0(".*_", classesToReplace), collapse = "|"),
                                   as.character(cohortDataNo34to36$ecoregionGroup)))
    if (!noCodesAre34to36)
      stop("lcc classes were mismanaged and some classes were not replaced correctly;
           contact developers: code 235")
  }
}

#' - `assertSppMaxBMaxANPP`: assert that all species have `maxB` and `maxANPP` values
#'   in the landscape;
#'
#' @template speciesEcoregion
#'
#' @template cohortData
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertSppMaxBMaxANPP <- function(speciesEcoregion,
                                 doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    tempDT <- speciesEcoregion[, .(sum(maxB, na.rm = TRUE),
                                   sum(maxANPP, na.rm = TRUE)),
                               by = speciesCode]
    if (any(tempDT$V1 == 0)) {
      noMaxB <- tempDT[V1 == 0, speciesCode]
      stop("The following species have no maxB values throughout the landscape:\n",
           paste(noMaxB, collapse = ", "), "\n",
           "This may be due to missing trait values.")
    }

    if (any(tempDT$V2 == 0)) {
      noMaxANPP <- tempDT[V2 == 0, speciesCode]
      stop("The following species have no maxANPP values throughout the landscape:\n",
           paste(noMaxANPP, collapse = ", "), "\n",
           "This may be due to missing trait values")
    }
  }
}

#' - `assertUniqueCohortData`: assert that `cohortData` has unique lines when subsetting
#' for a given set of columns;
#'
#' @param columns Vector of column names on which to test for unique `cohortData`
#'
#' @export
#' @rdname assertions
assertUniqueCohortData <- function(cohortData, columns,
                                   doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    obj <- cohortData[, .N, by = columns]
    whNEQOne <- which(obj$N != 1)
    test1 <- length(whNEQOne) == 0
    if (!test1)
      stop("There are identical cohorts (based on ", paste(columns, collapse = ", "),
           ") within the same pixelGroup")
  }
}

#' - `assertERGs`: assert that `ecoregionGroups` match across different objects;
#'
#' @template ecoregionMap
#' @template speciesEcoregion
#' @template minRelativeB
#'
#' @export
#' @rdname assertions
assertERGs <- function(ecoregionMap, cohortData, speciesEcoregion, minRelativeB,
                       doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    erg <- list()
    if (!missing(ecoregionMap)) {
      erg[[1]] <- sort(na.omit(unique(factorValues2(ecoregionMap, as.vector(ecoregionMap[]),
                                                    att = "ecoregionGroup"))))
      if (is.character(erg[[1]])) { # this can happen if SpatRaster
        erg[[1]] <- factor(erg[[1]])
      }
    }
    if (!missing(cohortData))
      erg[[2]] <- sort(unique(cohortData$ecoregionGroup))
    if (!missing(speciesEcoregion))
      erg[[3]] <- sort(unique(speciesEcoregion$ecoregionGroup))
    if (!missing(minRelativeB))
      erg[[4]] <- sort(unique(minRelativeB$ecoregionGroup))

    ## first test will detect differences in both factor levels and values
    lens <- vapply(erg, function(x) length(x) > 1, FUN.VALUE = logical(1))
    erg <- erg[lens]
    # test3 <- all(vapply(seq(erg)[-1], function(x)
    # identical(erg[[1]], erg[[x]]), FUN.VALUE = logical(1)))
    test3 <- all(outer(erg, erg, FUN = Vectorize(identical)))

    ## second test only detects differences in values
    erg <- lapply(erg, as.character)
    # test4 <- all(vapply(seq(erg)[-1], function(x)
    #   identical(erg[[1]], erg[[x]]), FUN.VALUE = logical(1)))
    test4 <- all(outer(erg, erg, FUN = Vectorize(identical)))

    if (!test4) {
      message(str(erg, 1))
      stop("speciesEcoregion, cohortData, and ecoregionMap should all have exactly the same",
           "\n  ecoregionGroups. They do not. This needs to be fixed before proceeding.")
    }

    ## only necessary to check this one if the first doesn't fail.
    if (!test3) {
      message(str(erg, 1))
      stop("speciesEcoregion, cohortData, and ecoregionMap should all have exactly the same",
           "\n  ecoregionGroup levels. They do not. This needs to be fixed before proceeding.")
    }
  }
}

#' - `assertColumns`: assert that an object contains a particular set of columns;
#'
#' @param obj A `data.frame`- or `data.table`-like object
#' @param colClasses A named vector of column classes, where the names are the column names
#'
#' @export
#' @rdname assertions
assertColumns <- function(obj, colClasses,
                          doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    ## put them in order of the object
    colNames <- names(colClasses)
    test1Indiv <- colNames %in% colnames(obj)
    test1 <- all(test1Indiv)
    mess1 <- character()
    if (!test1) {
      mess1 <- paste0("obj has missing column(s): ", paste(collapse = ", ", colNames[!test1Indiv]))
    }
    colClasses2 <- colClasses[na.omit(match(names(obj), names(colClasses)))]
    colNames2 <- names(colClasses2)
    test2Indiv <- vapply(seq_len(NCOL(obj[, colNames2, with = FALSE])),
                         function(colNum) {
                           is(obj[, colNames2[colNum], with = FALSE][[1]],
                              colClasses2[colNum])
                         }, FUN.VALUE = logical(1))
    test2 <- all(test2Indiv)
    mess2 <- character()
    if (!test2) {
      wh <- which(!test2Indiv)
      wrongCols <- colClasses2[wh]
      mess2 <- paste("obj column classes need correction.",
                      "It should have the following class(es):",
                      paste(names(wrongCols), "=", wrongCols, collapse = ", "))
    }
    if (!test1 || !test2) {
      stop(paste(mess1, mess2, sep = "\n"))
    }
  }
}

#' - `assertCohortData`: test that `pixelGroups` in `pixelGroupMap` and `cohortData` match;
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @template cohortDefinitionCols
#'
#' @param maxExpectedNumDiverge A numeric, length 1, indicating by how many they
#'   can diverge. Default 1.
#'
#' @param message An optional message to print. This may help identify where this function
#'   was called.
#'
#' @template verbose
#'
#' @export
#' @rdname assertions
assertCohortData <- function(cohortData, pixelGroupMap, maxExpectedNumDiverge = 1,
                             message = "", doAssertion = getOption("LandR.assertions", TRUE),
                             verbose = getOption("LandR.verbose", TRUE),
                             cohortDefinitionCols = c("pixelGroup", "age", "speciesCode")) {
  if (doAssertion) {
    if (!isTRUE("pixelGroup" %in% names(cohortData))) {
      stop("cohortData must have pixelGroup")
    }
    a <- sort(unique(na.omit(as.vector(pixelGroupMap[]))))
    b <- sort(unique(na.omit(cohortData$pixelGroup)))
    ## test1 and test2 can be 1 because there could be pixelGroup of 0, which is OK to not match
    test1 <- sum(!a %in% b)
    test2 <- sum(!b %in% a)

    cohortDataN <- cohortData[, .N, by = cohortDefinitionCols]
    test3 <- which(cohortDataN$N != 1)
    if (test1 > maxExpectedNumDiverge || test2 > maxExpectedNumDiverge) {
      if (nchar(message) > 0) message(message)
      if (verbose) {
        if (test1 > maxExpectedNumDiverge) message("test1 is ", test1,
                                                   " -- too many pixelGroups on pixelGroupMap")
        if (test2 > maxExpectedNumDiverge) message("test2 is ", test2,
                                                   " -- too many pixelGroups in cohortData")
      }
      stop("The sim$pixelGroupMap and cohortData have unmatching pixelGroup.",
           " They must be matching. Please contact the module developers")
    }
    if (length(test3) != 0)
      stop("There are duplicate, identical cohorts: ", print(cohortDataN[test3]))

    if (verbose > 1) {
      message(crayon::green("  -- assertion passed using assertCohortData --"))
    }
  }
}

#' - `assertPixelCohortData`: test that `pixelGroupMap` and `pixelCohortData` match `pixelIndex`;
#'
#' This is the full `pixelCohortData`, not the collapsed one.
#'
#' @template pixelCohortData
#' @template pixelGroupMap
#'
#' @export
#' @rdname assertions
assertPixelCohortData <- function(pixelCohortData, pixelGroupMap,
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    uniquePixelsInCohorts <- pixelGroupMap[][unique(pixelCohortData$pixelIndex)]
    pixelsOnMap <- sum(!is.na(as.vector(pixelGroupMap[])), na.rm = TRUE)
    lenUniquePixelsInCohorts <- length(unique(pixelCohortData$pixelIndex))
    lenBurnedPixels <- sum(as.vector(pixelGroupMap[]) == 0, na.rm = TRUE) # 30927
    pixelsRegeneratedOnZeros <- sum(uniquePixelsInCohorts == 0)
    allPixelsNotInCohortData <- as.vector(pixelGroupMap[])[-unique(pixelCohortData$pixelIndex)]
    numPixelsNoRegen <- sum(allPixelsNotInCohortData == 0, na.rm = TRUE)
    tableB <- sum(!is.na(allPixelsNotInCohortData)) # 25166

    test2 <- identical(as.integer(pixelsOnMap - tableB), lenUniquePixelsInCohorts)
    test3 <- identical(as.integer(pixelsOnMap - (lenBurnedPixels - pixelsRegeneratedOnZeros)),
                       lenUniquePixelsInCohorts)

    uniqueAllPixelsNotInCohortData <- unique(allPixelsNotInCohortData)
    ## after converting to `terra` NAs sometimes appear as NaNs
    test1 <- all(uniqueAllPixelsNotInCohortData %in% c(NA, NaN, 0L))
    if (!test1 || !test2 || !test3) {
      stop("Every value on pixelGroupMap >0 must have a pixelIndex in pixelCohortData.\n",
           "This test is failing, i.e., there are some pixelGroups in pixelGroupMap",
           " that aren't in pixelCohortData.")
    }
  }
}

#' - `assertSpeciesPlotLabels`: Check that each species as a unique label in the
#'   `EN_generic_short` and `Leading` columns of the `sppEquiv` table;
#'
#' @param speciesNames A vector of species names for which the labels will be checked
#'
#' @template sppEquiv
#'
#' @export
#' @rdname assertions
assertSpeciesPlotLabels <- function(speciesNames, sppEquiv,
                                    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    sppLabelsENshort <- equivalentName(speciesNames,  sppEquiv, "EN_generic_short")
    sppLabelsLeading <- equivalentName(speciesNames,  sppEquiv, "Leading")
    if (any(duplicated(sppLabelsENshort)) ||
            any(duplicated(sppLabelsLeading)))
      stop(paste(
      "Two or more species share the same label under the 'EN_generic_short' or",
      "'Leading' columns of 'sppEquiv'.\n",
      "Please provide unique 'EN_generic_short' and 'Leading' for each species."
     ))
  }
}

#' - `assertFireToleranceDif`: Assert that the difference between fire severity and
#'   species fire tolerances ranges between -4 and 4.;
#'
#' @template burnedPixelCohortData
#'
#' @export
#' @rdname assertions
assertFireToleranceDif <- function(burnedPixelCohortData,
                                   doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    test1 <- TRUE
    test2 <- TRUE
    if (min(burnedPixelCohortData$severityToleranceDif, na.rm = TRUE) < -4 ||
        max(burnedPixelCohortData$severityToleranceDif, na.rm = TRUE) > 4) {
      if (min(burnedPixelCohortData$firetolerance, na.rm = TRUE) < 1 ||
          min(burnedPixelCohortData$firetolerance, na.rm = TRUE) < 5)
        test1 <- FALSE
      if (min(burnedPixelCohortData$severity, na.rm = TRUE) < 1 ||
          min(burnedPixelCohortData$severity, na.rm = TRUE) < 5)
        test2 <- FALSE
    }
    if (!test1)
      stop("The difference between severity and species fire tolerance must be [-4,4].
           Fire tolerance has values outside of [1,5], please check your
           species traits table ('species')")
    if (!test2)
      stop("The difference between severity and species fire tolerance must be [-4,4].
           Severity has values outside of [1,5], please debug Biomass_regenerationPM")
    if (!test1 && !test2)
      stop("The difference between severity and species fire tolerance must be [-4,4].
           Severity and fire tolerances have values outside of [1,5], please debug
           Biomass_regenerationPM and check your species traits table ('species')")
  }
}

#' - `assertSpeciesLayers`: assert that species layers exist with species cover in the study area;
#'
#' @template speciesLayers
#' @param thresh the minimum number of pixels where the species must have
#'   `biomass > 0` to be considered present in the study area.
#'   Defaults to 1.
#'
#' @export
#' @rdname assertions
assertSpeciesLayers <- function(speciesLayers, thresh,
                                doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    ## covert to list if not a stack or single/multi layer SpatRaster
    if (is(speciesLayers, "RasterStack")) {
      speciesLayers <- unstack(speciesLayers)
    }
    ## convert to list if only one RasterLayer
    if (is(speciesLayers, "RasterLayer")) {
      speciesLayers <- list(speciesLayers)
    }

    test1 <- vapply(speciesLayers, FUN = function(x) {
      all(is.na(as.vector(values(x))))
    }, FUN.VALUE = logical(1))

    if (all(test1))
      stop("no pixels found were found with species % cover >=", thresh,
           ". Try lowering the threshold.")
  }
}

#' - `assertRstLCChange`: assert that `rstLCChange` is a mask-type raster layer
#'   and matches `rasterToMatch`;
#'
#' @param rstLCChange a raster layer indicating pixels were land-use change occurred as 1s
#' @template rasterToMatch
#'
#' @export
#' @rdname assertions
assertRstLCChange <- function(rstLCChange, rasterToMatch,
                              doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    ## check conformity with RTM
    if (!.compareRas(rstLCChange,
                     rasterToMatch, stopOnError = FALSE)) {
      stop("'rstLCChange' and 'rasterToMatch' differ in
         their properties. Please check")
    }

    ## check if it's a maks
    temp <- setdiff(as.vector(rstLCChange[]), c(1, NA))
    if (length(temp)) {
      stop("rstLCChange should be a 'mask', with 1s in disturbed pixels and NAs everywhere else")
    }
  }
}

#' - `assertSpeciesEcoregionCohortDataMatch`: assert that the `cohortData` and `speciesEcoregion`
#'   have matching classes (specifically, whether all combinations of `ecoregionGroup`
#'   and `speciesCode` are in both objects, no more no less);
#'
#' @template cohortData
#'
#' @template speciesEcoregion
#'
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertSpeciesEcoregionCohortDataMatch <- function(
    cohortData, speciesEcoregion,
    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    a <- setdiff(unique(paste(speciesEcoregion$ecoregionGroup, speciesEcoregion$speciesCode)),
                 unique(paste(cohortData$ecoregionGroup, cohortData$speciesCode)))

    if (length(a) > 0)
      stop("speciesEcoregion has ecoregionGroup x speciesCode values that are not in cohortData: ",
           paste(a, collapse = ", "))
    b <- setdiff(unique(paste(cohortData$ecoregionGroup, cohortData$speciesCode)),
                 unique(paste(speciesEcoregion$ecoregionGroup, speciesEcoregion$speciesCode)))
    if (length(b) > 0)
      stop("cohortData has ecoregionGroup x speciesCode values that are not in speciesEcoregion: ",
           paste(b, collapse = ", "))
  }
}

#' - `assertPixelCohortDataValid`: assert that the `standCohortData` has no `NA` values;
#'
#' @param standCohortData A `data.table` with simulated and observed stand data for validation
#'
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertPixelCohortDataValid <- function(standCohortData,
                                       doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    test <- standCohortData[, vapply(.SD, FUN = function(x) any(is.na(x)), FUN.VALUE = logical(1))]
    if (any(test)) {
      stop("There are NAs in either the observed or simulated data.\n",
           "Please check 'standCohortData'.")
    }

    test2 <- standCohortData[, sum(relativeAbundObsrvd), by = .(rep, year, pixelIndex)]$V1 |>
      unique()
    test3 <- standCohortData[, sum(relativeAbund), by = .(rep, year, pixelIndex)]$V1 |>
      unique()

    ## need to round, because some values "appear" to be 1, but probably have v. small decimals.
    if (length(setdiff(round(test2, 6), c(1, 0)))) {
      stop("Observed relative abundances do not sum to 1 (per pixelIndex/rep/year)")
    }

    if (length(setdiff(round(test3, 6), c(1, 0)))) {
      stop("Simulated relative abundances do not sum to 1 (per pixelIndex/rep/year)")
    }
  }
}

#' - `assertRepsAllCohortData`: assert that the `allCohortData` has the expected years
#'   and reps combinations;
#'
#' @param allCohortData A `data.table` with all simulated `cohortData` to use for
#'                      validation
#' @param reps repetition ids
#' @param years years
#'
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertRepsAllCohortData <- function(allCohortData, reps, years,
                                    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    test1 <- allCohortData[rep == reps[1] & year == years[1]]
    set(test1, NULL, "rep", NULL)
    out <- vapply(reps[-1], FUN = function(x, test1) {
      test2 <- allCohortData[rep == x & year == years[1]]
      set(test2, NULL, "rep", NULL)

      if (!identical(test1, test2)) {
        stop(paste("Simulation starting conditions are not identical between reps",
                   reps[1], "and", x))
      } else {
        TRUE
      }
    }, test1 = test1, FUN.VALUE = logical(1))
  }
}


#' - `assertStandAgeMapAttr`: assert that `standAgeMap` is a `RasterLayer` with
#'   attribute `imputedPixID`;
#'
#' @template standAgeMap
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertStandAgeMapAttr <- function(standAgeMap,
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    if (!inherits(standAgeMap, c("RasterLayer", "SpatRaster"))) {
      stop("standAgeMap should be a RasterLayer or SpatRaster")
    }
    if (is.null(attr(standAgeMap, "imputedPixID"))) {
      stop("standAgeMap should have a 'imputedPixID' attribute.",
           " If no pixel ages were imputed, please set attribute to `integer(0)`")
    } else {
      if (!(is(attr(standAgeMap, "imputedPixID"), "numeric") ||
            is(attr(standAgeMap, "imputedPixID"), "integer"))) {
        stop("standAgeMap attribute 'imputedPixID' should be numeric/integer")
      }
    }
  }
}

#' - `assertCohortDataAttr`: assert that `cohortData` has attribute `imputedPixID`;
#'
#' @template cohortData
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertCohortDataAttr <- function(cohortData,
                                 doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    if (is.null(attr(cohortData, "imputedPixID"))) {
      stop("cohortData should have a 'imputedPixID' attribute")
    } else {
      if (!(is(attr(cohortData, "imputedPixID"), "numeric") ||
            is(attr(cohortData, "imputedPixID"), "integer"))) {
        stop("cohortData attribute 'imputedPixID' should be numeric/integer")
      }
    }
  }
}

#' `assertSppVectors`: assert that vectors of species match.
#' It compares the species listed in the `sppEquiv` table (under the `sppEquivCol` column),
#' with species listed in the `sppNameVector` and `sppColorVect`.
#' The comparison excludes the "Mixed" type, if it exists.
#'
#' @param sppEquiv an table with a column containing species names
#'
#' @template sppNameVector
#'
#' @param sppEquivCol the column name to use from `sppEquiv`.
#'
#' @template sppColorVect
#'
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assertSppVectors <- function(sppEquiv = NULL, sppNameVector = NULL, sppColorVect = NULL,
                             sppEquivCol = NULL,
                             doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    if (is.null(sppEquiv)) {
      stop("Please provide 'sppEquiv' and at least one vector")
    }
    if (is.null(sppEquivCol)) {
      stop("Please provide 'sppEquivCol'")
    }

    if (is.null(sppNameVector) && is.null(sppColorVect)) {
      stop("Please provide 'sppNameVector' and/or 'sppColorVect'")
    }

    sppInSppEquiv <- unique(sppEquiv[[sppEquivCol]])

    if (!is.null(sppNameVector)) {
      test1 <- any(length(union(sppInSppEquiv, sppNameVector)) != length(sppInSppEquiv),
                   length(union(sppInSppEquiv, sppNameVector)) != length(sppNameVector))

      if (test1) {
        stop("Species listed in 'sppEquiv' and 'sppNameVector' differ.")
      }
    }

    if (!is.null(sppColorVect)) {
      sppInColourV <- setdiff(names(sppColorVect), "Mixed")
      test2 <- any(length(union(sppInSppEquiv, sppInColourV)) != length(sppInColourV),
                   length(union(sppInSppEquiv, sppInColourV)) != length(sppInSppEquiv))

      if (test2) {
        stop("Species listed in 'sppEquiv' and 'sppColorVect' differ (excluding 'Mixed' type).")
      }
    }
  }
}

#' Assert post-fire disturbance mortality and regeneration
#'
#' @param cohortDataOrig original `cohortData` (prior to any modification)
#'
#' @param pixelGroupMapOrig original `pixelGroupMap` (prior to any modification)
#'
#' @param cohortDataNew modified `cohortData` output from `updateCohortData()`
#'
#' @param pixelGroupMapNew modified `pixelGroupMap` output from `updateCohortData()`
#'
#' @param postDistPixelCohortData modified `cohortData` output from `updateCohortData()`
#'
#' @param distrbdPixelCohortData `cohortData`-like object containing all dead, surviving and
#' new cohorts (i.e. activated by serotiny/resprouting)
#'
#' @template doAssertion
#'
#' @return NULL
#' @export
assertPostPartialDist <- function(cohortDataOrig, pixelGroupMapOrig,
                                  cohortDataNew, pixelGroupMapNew,
                                  postDistPixelCohortData, distrbdPixelCohortData,
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    oldPCohortData <- addPixels2CohortData(cohortDataOrig, pixelGroupMapOrig, doAssertion = FALSE)
    newPCohortData <- addPixels2CohortData(cohortDataNew, pixelGroupMapNew)

    testPGs <- oldPCohortData[pixelGroup %in% postDistPixelCohortData$pixelGroup, pixelGroup]
    testPIs <- postDistPixelCohortData[pixelGroup %in% testPGs, pixelIndex]

    cols <- c("pixelIndex", "speciesCode", "age")
    test1 <- unique(newPCohortData[pixelIndex %in% testPIs, ..cols])
    test2 <- unique(postDistPixelCohortData[pixelIndex %in% testPIs, ..cols])
    setorderv(test1, cols)
    setorderv(test2, cols)

    if (!identical(test1, test2)) {
      stop("Post-fire disturbances miscalculated:  missing survivor/regenerated cohorts")
    }

    test3 <- unique(distrbdPixelCohortData[B == 0, ..cols])
    test4 <- unique(newPCohortData[, ..cols])
    if (nrow(newPCohortData[test3, on = cols, nomatch = 0L])) {
      stop("Killed cohorts are still in cohortData table")
    }
  }
}

#' Assert raw species table column types
#'
#' @param speciesTableRaw raw species traits `data.table` (see [getSpeciesTable()])
#'
#' @export
#' @rdname assertions
assertSpeciesTableRaw <- function(speciesTableRaw,
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  assertColumns(speciesTableRaw,
                c(LandisCode = "character",
                  Area = "factor",
                  Longevity = "integer",
                  Maturity = "integer",
                  Shade = "numeric",
                  Fire = "integer",
                  SeedEffDist = "integer",
                  SeedMaxDist = "integer",
                  VegProb = "numeric",
                  MinAgeVeg = "integer",
                  MaxAgeVeg = "integer",
                  PostFireRegen = "factor",
                  LeafLongevity = "integer",
                  WoodDecayRate = "numeric",
                  MortalityCurve = "numeric",
                  GrowthCurve = "numeric",
                  LeafLignin = "numeric",
                  HardSoft = "character"))
}

#' Assert species table column types
#'
#' @param speciesTable species traits `data.table` (see [prepSpeciesTable()])
#'
#' @export
#' @rdname assertions
assertSpeciesTable <- function(speciesTable,
                               doAssertion = getOption("LandR.assertions", TRUE)) {
  assertColumns(speciesTable,
                c(species = "character", Area = "factor", longevity = "integer",
                  sexualmature = "integer", shadetolerance = "numeric",
                  firetolerance = "integer", seeddistance_eff = "integer",
                  seeddistance_max = "integer", resproutprob = "numeric",
                  resproutage_min = "integer", resproutage_max = "integer",
                  postfireregen = "factor", leaflongevity = "integer",
                  wooddecayrate = "numeric", mortalityshape = "integer",
                  growthcurve = "numeric", leafLignin = "numeric", hardsoft = "factor"))
}
