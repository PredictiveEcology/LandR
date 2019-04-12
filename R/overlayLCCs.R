if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("ecoregionCode", "NAs"))
}

#' Overlay different LCC data sources
#'
#' @param LCCs A named list or named \code{RasterStack} of rasters whose content is Land Cover Class
#' @param forestedList A named list of same length and names as \code{LCCs} indicating
#'   which classes in each LCC raster are 'forested', either permanent or transient
#' @param outputLayer A character string that matches one of the named elements in \code{LCCs}.
#'   This will be the classification system returned.
#' @param NAcondition A character string  of a vectorized logical statement that will be
#'   parsed within the "forestEquivalencies" table.
#'   It should be a set of conditions with \code{== 0}, i.e., non-forested.
#'   Examples:, e.g., \code{"LCC2005 == 0"} or \code{"CC == 0 | LCC2005 == 0"},
#'   where \code{0} is the non-forested pixels based on converting LCCs and forestedList to 1s and 0s.
#' @param NNcondition A character string of a vectorized logical statement that will be parsed
#'   within the "forestEquivalencies" table.
#'   It should be a set of conditions with \code{== 0}, i.e., non-forested.
#'   Examples:, e.g., \code{"LCC2005 == 0"} or \code{"CC == 0 | LCC2005 == 0"},
#'   where \code{0} is the non-forested pixels based on converting LCCs and forestedList to 1s and 0s.
#' @param remapCondition Not yet implemented. This would be for a situation where
#'   2 LCC layers are provided, one has information in a pixel, but not the one
#'   which is \code{outputLayer}, so this needs a reclassify or remap.
#' @param classesToReplace Passed to \code{convertUnwantedLCC}, for the pixels where
#'   \code{NNcondition} is \code{TRUE}
#' @param availableERC_by_Sp Passed to \code{convertUnwantedLCC}, for the pixels where
#'   \code{NNcondition} is \code{TRUE}. If this is \code{NULL}, then it will be
#'   created internally with all pixels with:
#'   \code{data.table(initialEcoregionCode = LCCs[[outputLayer]][])}
#' @param forestEquivalencies A data.frame or NULL. If \code{NULL}, this function will
#'   derive this table automatically from the other arguments. Otherwise, the user must
#'   provide a data.frame with \code{length(LCCs) + 1} columns, and \code{2 ^ length(LCCs)}
#'   rows. Currently not used.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom data.table as.data.table
#' @importFrom raster nlayers
overlayLCCs <- function(LCCs, forestedList, outputLayer,
                        NAcondition, NNcondition, remapCondition,
                        classesToReplace, availableERC_by_Sp,
                        forestEquivalencies = NULL) {
  forestedListFail <- FALSE
  if (!missing(remapCondition)) warning("remapCondition is not yet implemented")
  if (is.null(names(forestedList))) forestedListFail <- TRUE
  if (!identical(sort(names(forestedList)), sort(names(LCCs))))
    forestedListFail <- TRUE

  if (isTRUE(forestedListFail))
    stop("forestedList must be named the same names as LCCs")

  # make sure they are in same order
  theOrder <- match(names(forestedList), names(LCCs))
  if (!identical(theOrder, seq_along(forestedList)))
    forestedList <- forestedList[theOrder]

  if (!is(LCCs, "RasterStack"))
    LCCs <- stack(LCCs)
  if (nlayers(LCCs) > 1) {

    forestedStack <- stack(LCCs)
    forestedStack[] <- 0
    forestedStack <- stack(forestedStack)
    names(forestedStack) <- names(LCCs)

    for(x in names(LCCs)) {
      forestedStack[[x]][LCCs[[x]][] %in% forestedList[[x]]] <- 1
    }

    namesForestedStack <- names(forestedStack)
    names(namesForestedStack) <- namesForestedStack
    dt <- as.data.table(lapply(namesForestedStack, function(x) forestedStack[[x]][]))
    dt[, ecoregionCode := LCCs[[outputLayer]][]]

    # 1. Put in NAs
    if (!missing(NAcondition)) {
      dt[, NAs := eval(parse(text = NAcondition)), with = TRUE]
      dt[NAs == TRUE, ecoregionCode := NA ]
      dt[, NAs := NULL]
    }

    if (is.null(availableERC_by_Sp)) {
      availableERC_by_Sp <- data.table(initialEcoregionCode = LCCs[[outputLayer]][])
    }

    if (!missing(NNcondition)) {

      dt[, NNcondition := eval(parse(text = NNcondition))]
      dt <- cbind(dt, availableERC_by_Sp)
      dt[, pixelIndex := seq(ncell(LCCs[[outputLayer]]))]
      dt1 <- dt[NNcondition == TRUE, c("initialEcoregionCode", "pixelIndex")]
      if (any(classesToReplace %in% dt1$initialEcoregionCode)) {
        # It is possible that there are no pixels with classesToReplace that are fulfill the NNcondition
        a <- convertUnwantedLCC(classesToReplace = classesToReplace,
                         rstLCC = LCCs[[outputLayer]],
                         theUnwantedPixels = dt1$pixelIndex,
                         availableERC_by_Sp = na.omit(dt)[, c("initialEcoregionCode", "pixelIndex")])
        dt <- a[dt, on = "pixelIndex"]
        dt[!is.na(ecoregionGroup), ecoregionCode := ecoregionGroup]
        dt[, ecoregionGroup := NULL]
      }
    }
    # replace all values in the raster
    LCCs[[outputLayer]][] <- dt$ecoregionCode
  }
  return(LCCs[[outputLayer]])
}
