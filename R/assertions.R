#' Assertions
#'
#' @param cohortData34to36 A cohortData data.table with only the pixels what were lcc 34:36
#' @param cohortData The full cohortData data.table
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

#' @param cohortData The full cohortData data.table
#' @param ecoregionMap The ecoregionMap, a raster of all the unique groupings
#' @param speciesEcoregion The speciesEcoregion data.table
#' @export
#' @rdname assertions
assertERGs <- function(ecoregionMap, cohortData, speciesEcoregion, minRelativeB) {
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

