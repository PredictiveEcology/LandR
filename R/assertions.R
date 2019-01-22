if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".N"))
}

#' Assertions
#'
#' @param cohortData34to36 A \code{cohortData} \code{data.table} with only the
#'                         pixels what were lcc 34:36
#'
#' @param cohortData The full \code{cohortData} \code{data.table}
#'
#' @export
#' @rdname assertions
assert1 <- function(cohortData34to36, cohortData) {
  if (getOption("LandR.assertions")) {
    allCodesAre34to36 <- all(grepl(".*34|.*35|.*36",
                                   as.character(cohortData34to36$initialEcoregionCode)))
    if (!allCodesAre34to36)
      stop("lcc classes were mismanaged; contact developers: code 234")

    onlyExistingCodes <- all(unique(cohortData34to36$ecoregionGroup) %in% unique(cohortData$initialEcoregionCode))
    if (!onlyExistingCodes)
      stop("There are some ecoregionCodes created post replacement of 34 and 35")
  }
}

#' @param columns Vector of column names on which to test for unique cohortData
#'
#' @export
#' @rdname assertions
assertUniqueCohortData <- function(cohortData, columns) {
  if (getOption("LandR.assertions")) {
    obj <- cohortData[ , .N, by = columns]
    whNEQOne <- which(obj$N != 1)
    test1 <- length(whNEQOne) == 0
    if (!test1)
      stop("There are identical cohorts (based on ", paste(columns, collapse = ", "),
           ") within the same pixelGroup")
  }
}

#' @param ecoregionMap The \code{ecoregionMap}, a raster of all the unique groupings
#' @param speciesEcoregion The \code{speciesEcoregion} \code{data.table}
#' @param minRelativeB TODO: add description
#'
#' @export
#' @importFrom utils str
#' @rdname assertions
assertERGs <- function(ecoregionMap, cohortData, speciesEcoregion, minRelativeB) {
  if (getOption("LandR.assertions", FALSE)) {
    erg <- list()
    if (!missing(ecoregionMap))
      erg[[1]] <- sort(na.omit(unique(factorValues2(ecoregionMap, ecoregionMap[], att = "ecoregionGroup"))))
    if (!missing(cohortData))
      erg[[2]] <- sort(unique(cohortData$ecoregionGroup))
    if (!missing(speciesEcoregion))
      erg[[3]] <- sort(unique(speciesEcoregion$ecoregionGroup))
    if (!missing(minRelativeB))
      erg[[4]] <- sort(unique(minRelativeB$ecoregionGroup))

    lens <- sapply(erg, function(x) length(x)>1)
    erg <- erg[lens]
    test3 <- all(sapply(seq(erg)[-1], function(x)
      identical(erg[[1]], erg[[x]])))

    if (!test3) {
      print(str(erg, 1))
      stop("speciesEcoregion, cohortData, and ecoregionMap should all have exactly the same",
           "\n  ecoregionGroups. They do not. This needs to be fixed before proceeding.")
    }
  }
}

#' @param obj A data.frame or data.table-like object
#' @param colClasses A named vector of column classes, where the names are the column names
#' @export
#' @rdname assertions
assertColumns <- function(obj, colClasses) {
  if (getOption("LandR.assertions")) {
    colNames <- names(colClasses)
    test1 <- all(colNames %in% colnames(obj))
    test2 <- all(sapply(seq(NCOL(obj[, colNames, with = FALSE])),
                        function(colNum) {
                          is(obj[, colNames[colNum], with = FALSE][[1]],
                             colClasses[colNum])
                        }))
    if (!test1 || !test2)
      stop("obj should be a data.table with at least ",NROW(colNames)," columns: ",
           paste(colNames, collapse = ", "),
           " ... of classes: ", paste(colClasses, collapse = ", "))
  }
}

#' A test that \code{pixelGroupMap} and \code{cohortData} match
#'
#' @inheritParams updateCohortData
#' @param sim If the simList is included, then the browser() call will be more useful
#' @param maxExpectedNumDiverge A numeric, length 1, indicating by how many they
#'   can diverge. Default 1.
#' @param message An optional message to print. This may help identify where this function
#'   was called.
#'
#' @note
#' TODO
#'
#' @export
#' @importFrom crayon green
#' @importFrom stats na.omit
assertCohortData <- function(cohortData, pixelGroupMap, sim, maxExpectedNumDiverge = 1,
                           message = "", verbose = getOption("LandR.verbose", TRUE)) {
  if (getOption("LandR.assertions", FALSE)) {
    a <- sort(unique(na.omit(pixelGroupMap[])))
    b <- sort(unique(na.omit(cohortData$pixelGroup)))
    test1 <- sum(!a %in% b)  # can be 1 because there could be pixelGroup of 0, which is OK to not match
    test2 <- sum(!b %in% a)  # can be 1 because there could be pixelGroup of 0, which is OK to not match

    browser(expr = exists("aaaa"))
    cohortDataN <- cohortData[, .N, by = c("pixelGroup", "speciesCode", "age", "B")]
    test3 <- which(cohortDataN$N != 1)
    if (test1 > maxExpectedNumDiverge || test2 > maxExpectedNumDiverge) {
      if (nchar(message) > 0) message(message)
      if (test1 > maxExpectedNumDiverge) message("test1 is ", test1, " -- too many pixelGroups on pixelGroupMap")
      if (test2 > maxExpectedNumDiverge) message("test2 is ", test2, " -- too many pixelGroups in cohortData")
      stop("The sim$pixelGroupMap and cohortData have unmatching pixelGroup. They must be matching.",
           " Please contact the module developers")
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
#' This is the full pixelCohortData, not the collapsed one
#'
#' @param pixelCohortData The full \code{cohortData} \code{data.table}
#'
#' @inheritParams updateCohortData
#'
#' @export
assertPixelCohortData <- function(pixelCohortData, pixelGroupMap) {
  if (isTRUE(getOption("LandR.assertions"))) {
    uniquePixelsInCohorts <- pixelGroupMap[][unique(pixelCohortData$pixelIndex)]
    pixelsOnMap <- sum(!is.na(pixelGroupMap[]), na.rm = TRUE)
    lenUniquePixelsInCohorts <- length(unique(pixelCohortData$pixelIndex))
    lenBurnedPixels <- sum(pixelGroupMap[] == 0, na.rm = TRUE) # 30927
    pixelsRegeneratedOnZeros <- sum(uniquePixelsInCohorts == 0)
    allPixelsNotInCohortData <- pixelGroupMap[][-unique(pixelCohortData$pixelIndex)]
    numPixelsNoRegen <- sum(allPixelsNotInCohortData == 0, na.rm = TRUE)
    tableB <- table(allPixelsNotInCohortData) # 25166

    test2 <- identical(as.integer(pixelsOnMap - tableB), lenUniquePixelsInCohorts)
    test3 <- identical(as.integer(pixelsOnMap - (lenBurnedPixels - pixelsRegeneratedOnZeros)), lenUniquePixelsInCohorts)

    uniqueAllPixelsNotInCohortData <- unique(allPixelsNotInCohortData)
    test1 <- all(uniqueAllPixelsNotInCohortData %in% c(NA, 0L))
    if (!test1 | !test2 | !test3) {
      stop("Every value on pixelGroupMap greater than 0 must have a pixelIndex in pixelCohortData.",
           " This test is failing, i.e., there are some pixelGroupMaps have pixelGroups, and aren't in pixelCohortData.")
    }
  }
}
