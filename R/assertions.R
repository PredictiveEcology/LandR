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
  if (getOption("LandR.assertions", FALSE)) {
    a <- sort(unique(na.omit(pixelGroupMap[])))
    b <- sort(unique(na.omit(cohortData$pixelGroup)))
    test1 <- sum(!a %in% b)  # can be 1 because there could be pixelGroup of 0, which is OK to not match
    test2 <- sum(!b %in% a)  # can be 1 because there could be pixelGroup of 0, which is OK to not match
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
      stop("There are duplicate, identical cohorts: ", cohortDataN[test3])

    message(crayon::green("  -- assertion passed using testCohortData --"))
  }
}