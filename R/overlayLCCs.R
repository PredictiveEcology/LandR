utils::globalVariables(c("ecoregionCode", "NAs", "newLCC"))

#' Overlay different LCC data sources
#'
#' @param LCCs A named list or named \code{RasterStack} of layers whose content
#'   is Land Cover Class.
#' @param forestedList A named list of same length and names as \code{LCCs} indicating
#'   which classes in each LCC raster are 'forested', either permanent or transient
#' @param outputLayer A character string that matches one of the named elements
#'   in \code{LCCs}. This will be the classification system returned.
#' @param NAcondition The condition when a pixel is deemed to be \code{NA}.
#'   Given as a character string of a vectorized logical statement that will be
#'   within the \code{forestEquivalencies} table.
#'   It should be a set of conditions with \code{== 0}, i.e., non-forested.
#'   Examples:, e.g., \code{"LCC2005 == 0"} or \code{"CC == 0 | LCC2005 == 0"},
#'   where \code{0} is the non-forested pixels based on converting LCCs and
#'   \code{forestedList} to \code{1} and \code{0}.
#' @param NNcondition The 'nearest-neighbour' condition; i.e., the condition when
#'   a nearest-neighbour search is done to fill in the pixel with forested type.
#'   Given as a character string of a vectorized logical statement that will be
#'   parsed within the \code{forestEquivalencies} table.
#'   It should be a set of conditions with \code{== 0}, i.e., non-forested.
#'   Examples:, e.g., \code{"LCC2005 == 0"} or \code{"CC == 0 | LCC2005 == 0"},
#'   where \code{0} is the non-forested pixels based on converting LCCs and
#'   \code{forestedList} to \code{1} and \code{0}.
#' @param remapTable \code{data.table}. This would be for a situation where
#'   2 LCC layers are provided, one has information in a pixel, but not the one
#'   which is \code{outputLayer}, so this needs a reclassify or remap.
#' @param classesToReplace Passed to \code{convertUnwantedLCC}, for the pixels where
#'   \code{NNcondition} is \code{TRUE}
#' @param availableERC_by_Sp Passed to \code{convertUnwantedLCC}, for the pixels where
#'   \code{NNcondition} is \code{TRUE}. If this is \code{NULL}, then it will be
#'   created internally with all pixels with:
#'   \code{data.table(initialEcoregionCode = LCCs[[outputLayer]][])}
#' @param forestEquivalencies A \code{data.frame} or \code{NULL}.
#'   If \code{NULL}, this function will derive this table automatically from the
#'   other arguments. Otherwise, the user must provide a \code{data.frame} with
#'   \code{length(LCCs) + 1} columns, and \code{2 ^ length(LCCs)} rows.
#'   Currently not used.
#'
#' @author Eliot McIntire and Alex Chubaty
#' @export
#' @importFrom data.table as.data.table
#' @importFrom raster nlayers stack
overlayLCCs <- function(LCCs, forestedList, outputLayer,
                        NAcondition, NNcondition, remapTable = NULL,
                        classesToReplace, availableERC_by_Sp,
                        forestEquivalencies = NULL) {
  forestedListFail <- FALSE
  if (is.null(names(forestedList))) forestedListFail <- TRUE
  if (!identical(sort(names(forestedList)), sort(names(LCCs))))
    forestedListFail <- TRUE

  if (isTRUE(forestedListFail))
    stop("forestedList must use the same names as LCCs")

  # make sure they are in same order
  theOrder <- match(names(forestedList), names(LCCs))
  if (!identical(theOrder, seq_along(forestedList)))
    forestedList <- forestedList[theOrder]

  if (!is(LCCs, "RasterStack"))
    LCCs <- stack(LCCs)

  if (nlayers(LCCs) > 1) {
    forestedStack <- stack(LCCs)
    forestedStack[] <- 0 ## will convert to a brick, so need to restack
    forestedStack <- stack(forestedStack)
    names(forestedStack) <- names(LCCs)

    for (x in names(LCCs)) {
      forestedStack[[x]][LCCs[[x]][] %in% forestedList[[x]]] <- 1
    }

    namesForestedStack <- names(forestedStack)
    names(namesForestedStack) <- namesForestedStack

    if (is.null(availableERC_by_Sp)) {
      availableERC_by_Sp <- data.table(initialEcoregionCode = LCCs[[outputLayer]][])
    }

    if (is.null(remapTable)) {
      dt <- as.data.table(lapply(namesForestedStack, function(x) forestedStack[[x]][]))
      dt[, ecoregionCode := LCCs[[outputLayer]][]]

      # 1. Put in NAs
      if (!missing(NAcondition)) {
        dt[, NAs := eval(parse(text = NAcondition)), with = TRUE]
        dt[NAs == TRUE, ecoregionCode := NA]
        dt[, NAs := NULL]
      }

      if (!missing(NNcondition)) {
        dt[, NNcondition := eval(parse(text = NNcondition))]
        dt <- cbind(dt, availableERC_by_Sp)
        dt[, pixelIndex := seq(ncell(LCCs[[outputLayer]]))]
        dt1 <- dt[NNcondition == TRUE, c("initialEcoregionCode", "pixelIndex")]
        if (any(classesToReplace %in% dt1$initialEcoregionCode)) {
          ## It's possible that there are no pixels with classesToReplace that fulfill the NNcondition
          a <- convertUnwantedLCC(classesToReplace = classesToReplace,
                                  rstLCC = LCCs[[outputLayer]],
                                  theUnwantedPixels = dt1$pixelIndex,
                                  availableERC_by_Sp = na.omit(dt)[, c("initialEcoregionCode", "pixelIndex")])
          dt <- a[dt, on = "pixelIndex"]
          dt[!is.na(ecoregionGroup), ecoregionCode := ecoregionGroup]
          dt[, ecoregionGroup := NULL]
        }
      }
    } else {
      if (is.data.frame(remapTable))
        as.data.table(remapTable)

      if (!all(names(LCCs) %in% colnames(remapTable)))
        stop("All LCC names must be columns in remapTable")

      namesLCCs <- names(LCCs)
      names(namesLCCs) <- namesLCCs
      dt <- as.data.table(lapply(namesLCCs, function(x) LCCs[[x]][]))
      dt[, ecoregionCode := LCCs[[outputLayer]][]]
      dt <- cbind(dt, availableERC_by_Sp)
      dt[, pixelIndex := seq(ncell(LCCs[[outputLayer]]))]

      setkeyv(dt, names(LCCs))
      setkeyv(remapTable, names(LCCs))
      dt2 <- merge(dt, remapTable, all.x = TRUE) ## TODO: use 'dt' instead of 'dt2'
      dt2[, ':='(initialEcoregionCode = newLCC, ecoregionCode = newLCC)]
      dt2[, newLCC := NULL]

      dt1 <- dt2[initialEcoregionCode %in% classesToReplace, c("initialEcoregionCode", "pixelIndex")]

      ## It's possible that there are no pixels with classesToReplace
      if (any(classesToReplace %in% dt1$initialEcoregionCode)) {
        a <- convertUnwantedLCC(classesToReplace = classesToReplace,
                                rstLCC = LCCs[[outputLayer]],
                                #theUnwantedPixels = dt1$pixelIndex,
                                availableERC_by_Sp = na.omit(dt2)[, c("initialEcoregionCode", "pixelIndex")])
        dt <- a[dt, on = "pixelIndex"]
        dt[!is.na(ecoregionGroup), ecoregionCode := ecoregionGroup]
        dt[, ecoregionGroup := NULL]
      }

      setkey(dt, pixelIndex)
    }

    # replace all values in the raster
    LCCs[[outputLayer]][] <- dt$ecoregionCode
  }
  return(LCCs[[outputLayer]])
}
