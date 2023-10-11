utils::globalVariables(c(
  "ecoregionCode", "NAs", "newLCC"
))

#' Overlay different LCC data sources
#'
#' @param LCCs A named list or named `RasterStack` of layers whose content
#'   is Land Cover Class.
#'
#' @param forestedList A named list of same length and names as `LCCs` indicating
#'   which classes in each LCC raster are 'forested', either permanent or transient
#'
#' @param outputLayer A character string that matches one of the named elements
#'   in `LCCs`. This will be the classification system returned.
#'
#' @param NAcondition The condition when a pixel is deemed to be `NA`.
#'   Given as a character string of a vectorized logical statement that will be
#'   within the `forestEquivalencies` table.
#'   It should be a set of conditions with `== 0`, i.e., non-forested.
#'   Examples:, e.g., `"LCC2005 == 0"` or `"CC == 0 | LCC2005 == 0"`,
#'   where `0` is the non-forested pixels based on converting LCCs and
#'   `forestedList` to `1` and `0`.
#'
#' @param NNcondition The 'nearest-neighbour' condition; i.e., the condition when
#'   a nearest-neighbour search is done to fill in the pixel with forested type.
#'   Given as a character string of a vectorized logical statement that will be
#'   parsed within the `forestEquivalencies` table.
#'   It should be a set of conditions with `== 0`, i.e., non-forested.
#'   Examples:, e.g., `"LCC2005 == 0"` or `"CC == 0 | LCC2005 == 0"`,
#'   where `0` is the non-forested pixels based on converting LCCs and
#'   `forestedList` to `1` and `0`.
#'
#' @param remapTable `data.table`. This would be for a situation where
#'   2 LCC layers are provided, one has information in a pixel, but not the one
#'   which is `outputLayer`, so this needs a reclassify or remap.
#'
#' @param classesToReplace Passed to [convertUnwantedLCC()], for the pixels where
#'   `NNcondition` is `TRUE`
#'
#' @param availableERC_by_Sp Passed to [convertUnwantedLCC()], for the pixels where
#'   `NNcondition` is `TRUE`. If this is `NULL`, then it will be
#'   created internally with all pixels with:
#'   `data.table(initialEcoregionCode = LCCs[[outputLayer]][])`
#'
#' @param forestEquivalencies A `data.frame` or `NULL`.
#'   If `NULL`, this function will derive this table automatically from the
#'   other arguments. Otherwise, the user must provide a `data.frame` with
#'   `length(LCCs) + 1` columns, and `2 ^ length(LCCs)` rows.
#'   Currently not used.
#'
#' @author Eliot McIntire and Alex Chubaty
#' @export
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

  fn <- rasterRead
  if (!is(LCCs, "RasterStack") && !identical(fn, terra::rast)) {
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster pkg is not installed; ",
           "either set options('reproducible.rasterRead' = 'terra::rast'), ",
           "or install.packages('raster').")
    }
    fn <- raster::stack
  }
  LCCs <- fn(LCCs)


  if (length(names(LCCs)) > 1) {
    forestedStack <- fn(LCCs)
    forestedStack[] <- 0 ## will convert to a brick, so need to restack
    forestedStack <- fn(forestedStack)
    names(forestedStack) <- names(LCCs)

    for (x in names(LCCs)) {
      forestedStack[[x]][LCCs[[x]][] %in% forestedList[[x]]] <- 1
    }

    namesForestedStack <- names(forestedStack)
    names(namesForestedStack) <- namesForestedStack

    if (is.null(availableERC_by_Sp)) {
      availableERC_by_Sp <- data.table(initialEcoregionCode = as.vector(LCCs[[outputLayer]][]))
    }

    if (is.null(remapTable)) {
      dt <- as.data.table(lapply(namesForestedStack, function(x) as.vector(forestedStack[[x]][])))
      dt[, ecoregionCode := as.vector(LCCs[[outputLayer]][])]

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
        remapTable <- as.data.table(remapTable)

      if (!all(names(LCCs) %in% colnames(remapTable)))
        stop("All LCC names must be columns in remapTable")

      namesLCCs <- names(LCCs)
      names(namesLCCs) <- namesLCCs
      dt <- as.data.table(lapply(namesLCCs, function(x) as.vector(LCCs[[x]][])))
      dt[, ecoregionCode := as.vector(LCCs[[outputLayer]][])]
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
  if (is(LCCs[[outputLayer]], "SpatRaster"))
    LCCs[[outputLayer]] <- as.int(LCCs[[outputLayer]])

  return(LCCs[[outputLayer]])
}
