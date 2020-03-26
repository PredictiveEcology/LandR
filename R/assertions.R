if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".N"))
}

#' Assertions
#'
#' Assert that \code{ecoregionCodes} that were replaced, were correctly identified.
#'
#' @param cohortData34to36 A \code{cohortData} \code{data.table} with only the
#'                         pixels what were LCC 34:36
#' @param rmZeroBiomassQuote  An expression to evaluate, in the form of \code{quote(B>0)},
#'    used to select cohorts with biomass.
#'
#' @template cohortData
#' @template doAssertion
#'
#' @export
#' @rdname assertions
assert1 <- function(cohortData34to36, cohortData, rmZeroBiomassQuote,
                    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    allCodesAre34to36 <- all(grepl(".*34|.*35|.*36",
                                   as.character(cohortData34to36$initialEcoregionCode)))
    if (!allCodesAre34to36)
      stop("lcc classes were mismanaged; contact developers: code 234")

    ## do they match the codes found for areas with biomass?
    ## if not, is this because the non-matching codes had no biomass to start with?
    onlyExistingCodes <- all(unique(cohortData34to36$ecoregionGroup) %in%
                               unique(cohortData[eval(rmZeroBiomassQuote), initialEcoregionCode]))
    if (!onlyExistingCodes) {
      temp1 <- unique(cohortData34to36$ecoregionGroup)
      temp2 <- unique(cohortData[eval(rmZeroBiomassQuote), initialEcoregionCode])
      nonMatching <- setdiff(temp1, temp2)

      ## was there any biomass/species data in these pixels?
      nonMatchingB <- sum(cohortData[initialEcoregionCode %in% nonMatching, B], na.rm = TRUE)

      if (nonMatchingB)
        stop("There are some ecoregionCodes created post replacement of 34 and 35")
    }
  }
}

#' Assert that \code{cohortData} has unique lines when subsetting for a given set of columns
#'
#' @param columns Vector of column names on which to test for unique \code{cohortData}
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

#' Assert that \code{ecoregionGroups} match across different objects
#'
#' @param ecoregionMap The \code{ecoregionMap}, a raster of all the unique groupings
#' @template speciesEcoregion
#' @param minRelativeB TODO: add description
#'
#' @export
#' @importFrom utils str
#' @rdname assertions
assertERGs <- function(ecoregionMap, cohortData, speciesEcoregion, minRelativeB,
                       doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    erg <- list()
    if (!missing(ecoregionMap))
      erg[[1]] <- sort(na.omit(unique(factorValues2(ecoregionMap, ecoregionMap[],
                                                    att = "ecoregionGroup"))))
    if (!missing(cohortData))
      erg[[2]] <- sort(unique(cohortData$ecoregionGroup))
    if (!missing(speciesEcoregion))
      erg[[3]] <- sort(unique(speciesEcoregion$ecoregionGroup))
    if (!missing(minRelativeB))
      erg[[4]] <- sort(unique(minRelativeB$ecoregionGroup))

    erg <- lapply(erg, as.character)
    lens <- sapply(erg, function(x) length(x) > 1)
    erg <- erg[lens]
    test3 <- all(sapply(seq(erg)[-1], function(x)
      identical(erg[[1]], erg[[x]])))

    if (!test3) {
      message(str(erg, 1))
      stop("speciesEcoregion, cohortData, and ecoregionMap should all have exactly the same",
           "\n  ecoregionGroups. They do not. This needs to be fixed before proceeding.")
    }
  }
}


#' Assert that an object contains a particular set of columns
#'
#' @param obj A data.frame or data.table-like object
#' @param colClasses A named vector of column classes, where the names are the column names
#'
#' @export
#' @rdname assertions
assertColumns <- function(obj, colClasses,
                          doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    colNames <- names(colClasses)
    test1 <- all(colNames %in% colnames(obj))
    test2 <- all(sapply(seq(NCOL(obj[, colNames, with = FALSE])),
                        function(colNum) {
                          is(obj[, colNames[colNum], with = FALSE][[1]],
                             colClasses[colNum])
                        }))
    if (!test1 || !test2)
      stop("obj should be a data.table with at least ", NROW(colNames), " columns: ",
           paste(colNames, collapse = ", "),
           " ... of classes: ", paste(colClasses, collapse = ", "))
  }
}

#' A test that \code{pixelGroupMap} and \code{cohortData} match
#'
#' @template cohortData
#'
#' @template pixelGroupMap
#'
#' @param sim If the \code{simList} is included, then the \code{browser()} call will be more useful.
#'
#' @param maxExpectedNumDiverge A numeric, length 1, indicating by how many they
#'   can diverge. Default 1.
#'
#' @param message An optional message to print. This may help identify where this function
#'   was called.
#'
#' @param verbose Controls message output. Defaults to \code{getOption("LandR.verbose")}
#'
#' @note
#' TODO
#'
#' @export
#' @importFrom crayon green
#' @importFrom stats na.omit
#' @rdname assertions
assertCohortData <- function(cohortData, pixelGroupMap, sim, maxExpectedNumDiverge = 1,
                             message = "", doAssertion = getOption("LandR.assertions", TRUE),
                             verbose = getOption("LandR.verbose", TRUE)) {
  if (doAssertion) {
    if (!isTRUE("pixelGroup" %in% names(cohortData))) {
      stop("cohortData must have pixelGroup")
    }
    a <- sort(unique(na.omit(pixelGroupMap[])))
    b <- sort(unique(na.omit(cohortData$pixelGroup)))
    ## test1 and test2 can be 1 because there could be pixelGroup of 0, which is OK to not match
    test1 <- sum(!a %in% b)
    test2 <- sum(!b %in% a)

    browser(expr = exists("aaaa"))
    cohortDataN <- cohortData[, .N, by = c("pixelGroup", "speciesCode", "age", "B")]
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

#' A test that \code{pixelGroupMap} and \code{pixelCohortData} match \code{pixelIndex}
#'
#' This is the full \code{pixelCohortData}, not the collapsed one.
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
    pixelsOnMap <- sum(!is.na(pixelGroupMap[]), na.rm = TRUE)
    lenUniquePixelsInCohorts <- length(unique(pixelCohortData$pixelIndex))
    lenBurnedPixels <- sum(pixelGroupMap[] == 0, na.rm = TRUE) # 30927
    pixelsRegeneratedOnZeros <- sum(uniquePixelsInCohorts == 0)
    allPixelsNotInCohortData <- pixelGroupMap[][-unique(pixelCohortData$pixelIndex)]
    numPixelsNoRegen <- sum(allPixelsNotInCohortData == 0, na.rm = TRUE)
    tableB <- sum(!is.na(allPixelsNotInCohortData)) # 25166

    test2 <- identical(as.integer(pixelsOnMap - tableB), lenUniquePixelsInCohorts)
    test3 <- identical(as.integer(pixelsOnMap - (lenBurnedPixels - pixelsRegeneratedOnZeros)),
                       lenUniquePixelsInCohorts)

    uniqueAllPixelsNotInCohortData <- unique(allPixelsNotInCohortData)
    test1 <- all(uniqueAllPixelsNotInCohortData %in% c(NA, 0L))
    if (!test1 | !test2 | !test3) {
      stop("Every value on pixelGroupMap >0 must have a pixelIndex in pixelCohortData.\n",
           "This test is failing, i.e., there are some pixelGroups in pixelGroupMap",
           " that aren't in pixelCohortData.")
    }
  }
}

#' Check that each species as a unique label in the 'EN_generic_short' and 'Leading'
#'   columns of the \code{sppEquiv} table.
#'
#' @param speciesNames A vector of species names for which the labels will be checked
#' @template sppEquiv
#'
#' @export
#' @rdname assertions
assertSpeciesPlotLabels <- function(speciesNames, sppEquiv,
                                    doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    sppLabelsENshort <- equivalentName(speciesNames,  sppEquiv, "EN_generic_short")
    sppLabelsLeading <- equivalentName(speciesNames,  sppEquiv, "Leading")
    if (any(duplicated(sppLabelsENshort) |
            duplicated(sppLabelsLeading)))
      stop("2 or more species share the same label under the 'EN_generic_short' or 'Leading' columns of sim$sppEquiv.
          Please provide unique 'EN_generic_short' and 'Leading' for each species")
  }
}

#' Assert that the difference between fire severity and species fire tolerances
#'  ranges between -4 and 4.
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
    if (min(burnedPixelCohortData$severityToleranceDif, na.rm = TRUE) < -4 |
        max(burnedPixelCohortData$severityToleranceDif, na.rm = TRUE) > 4) {
      if (min(burnedPixelCohortData$firetolerance, na.rm = TRUE) < 1 |
          min(burnedPixelCohortData$firetolerance, na.rm = TRUE) < 5 )
        test1 <- FALSE
      if (min(burnedPixelCohortData$severity, na.rm = TRUE) < 1 |
          min(burnedPixelCohortData$severity, na.rm = TRUE) < 5 )
        test2 <- FALSE
    }
    if (!test1)
      stop("The difference between severity and species fire tolerance must be [-4,4].
           Fire tolerance has values outside of [1,5], please check your
           species traits table ('species')")
    if (!test2)
      stop("The difference between severity and species fire tolerance must be [-4,4].
           Severity has values outside of [1,5], please debug Biomass_regenerationPM")
    if (!test1 & !test2)
      stop("The difference between severity and species fire tolerance must be [-4,4].
           Severity and fire tolerances have values outside of [1,5], please debug
           Biomass_regenerationPM and check your species traits table ('species')")
  }
}


#' Assert that species layers exist with species cover in the study area
#'
#' @param speciesLayers A \code{RasterStack} or \code{RasterLayer} that
#'   should contain species cover data in the study area
#' @param thresh the minimum number of pixels where the species must have
#'   \code{biomass > 0} to be considered present in the study area.
#'   Defaults to 1.
#'
#' @export
#' @rdname assertions
assertSpeciesLayers <- function(speciesLayers, thresh,
                                doAssertion = getOption("LandR.assertions", TRUE)) {
  if (doAssertion) {
    ## covert to list if not a stack
    if (class(speciesLayers) != "RasterStack") {
      speciesLayers <- list(speciesLayers)

      test1 <- sapply(speciesLayers, FUN = function(x)
        all(is.na(getValues(x))))
    } else {
      test1 <- sapply(1:nlayers(speciesLayers), FUN = function(x)
        all(is.na(getValues(speciesLayers[[x]]))))
    }

    if (all(test1))
      stop("no pixels found were found with species % cover >=", thresh,
           ". Try lowering the threshold.")
  }
}

#' Assert that \code{rstLCChange} is a mask-type raster layer and matches
#' rasterToMatch
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
    if (!compareRaster(rstLCChange,
                       rasterToMatch, stopiffalse = FALSE)) {
      stop("'rstLCChange' and 'rasterToMatch' differ in
         their properties. Please check")
    }

    ## check if it's a maks
    temp <- setdiff(getValues(rstLCChange), c(1, NA))
    if (length(temp))
      stop("rstLCChange should be a 'mask', with 1s in disturbed pixels and NAs everywhere else")

  }
}


#' Assert that the cohortData speciesEcoregion have matching clases
#'
#' Specifically, whether all combinations of ecoregionGroup and speciesCode are in both objects, no more
#' no less.
#'
#' @param cohortData A cohortData object
#' @param speciesEcoregion A speciesEcoregion object
#' @export
assertSpeciesEcoregionCohortDataMatch <- function(cohortData, speciesEcoregion,
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
