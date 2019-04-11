#' Overlay different LCC data sources
#'
#' @param LCCs A named list or named \code{RasterStack} of rasters whose content is Land Cover Class
#' @param forestedList A named list of same length and names as \code{LCCs} indicating
#'   which classes in each LCC raster are 'forested', either permanent or transient
#' @param outputLayer A character string that matches one of the named elements in \code{LCCs}.
#'   This will be the classification system returned.
#' @param NAcombos A character string that will be parsed within the "forestEquivalencies" table.
#'   It should be a set of conditions that == 0, e.g., \code{"LCC2005 == 0"} or \code{"CC == 0 | LCC2005 == 0"},
#'   where \code{0} is the non-forested pixels based on converting LCCs and forestedList to 1s and 0s.
#' @param forestEquivalencies A data.frame or NULL. If \code{NULL}, this function will
#'   derive this table automatically from the other arguments. Otherwise, the user must
#'   provide a data.frame with \code{length(LCCs) + 1} columns, and \code{2 ^ length(LCCs)}
#'   rows. Currently not used.
#'
#' @author Eliot McIntire
#' @export
overlayLCCs <- function(LCCs, forestedList, outputLayer,
                        NAcombos,
                        forestEquivalencies = NULL) {
  forestedListFail <- FALSE
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

    browser()

    dt[, NAs := eval(parse(text = NAcombos)), with = TRUE]

    # Next is convertUnwanted and remap
  }
}
