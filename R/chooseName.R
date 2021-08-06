#' Return equivalent name from a \code{data.frame} of equivalencies
#'
#' This is simply a wrapper around \code{match} or \code{\%in\%} for
#' a specific \code{data.frame} of values.
#'
#' @param value   Vector of values to match in \code{df}.
#' @param df      A \code{data.frame} where every row is a set of equivalent names.
#' @param column  A character string or numeric of length 1, indicating the column
#'                in \code{df} to return names from.
#' @param multi   Logical. If \code{TRUE}, then all matches will be returned. Default
#'                is \code{FALSE} for backwards compatibility. This may result in more
#'                elements than were in \code{value}, as each \code{value} may be matched
#'                by more than one entry, returning more than one result each from \code{column}.
#'                If a species is found in the \code{df}, but there is no corresponding
#'                value in \code{column}, it will return \code{NA}
#' @param searchColumn Optionally, provide the name of a column in \code{df} that results
#'                must be found in. The return value will still be from \code{column}
#'
#' @export
#' @rdname equivalentName
equivalentName <- function(value, df, column, multi = FALSE, searchColumn = NULL) {
  out <- equivalentNameAsList(value, df, multi)
  likelyMatch <- equivalentNameColumn(value, df, column, multi = multi, searchColumn = searchColumn)
  df[[column]][out[[likelyMatch]]]
}

#' \code{equivalentNameColumn} will provide the name of the column that best matches.
#' @export
#' @rdname equivalentName
equivalentNameColumn <- function(value, df, column, multi = FALSE, searchColumn = NULL) {
  out <- equivalentNameAsList(value, df, multi)
  likelyMatch <- if (is.null(searchColumn)) {
    names(which.max(unlist(lapply(out, function(x) sum(!is.na(x))))))
  } else {
    searchColumn
  }
  likelyMatch
}

equivalentNameAsList <- function(value, df, multi) {
  out <- lapply(df, function(x) {
    if (isTRUE(multi)) {
      which(x %in% as.character(value))
    } else {
      match(as.character(value), x)
    }
  })
}


#' Check and expand sppEquiv
#'
#' This will expand a \code{sppEquiv} object that is only a vector or only a one-column
#' data.table into a many column data.table, if the columns that are present do not
#' contain \code{ensureColumns}.
#'
#' @param sppEquiv A character vector or data.table with named column(s). If this
#'   data.table does not have columns named \code{ensureColumns}, then it will
#'   attempt to merge this data.table with LandR::sppEquivalencies_CA to get
#'   \code{ensureColumns}.
#' @param ensureColumns A character vector of column names that must be in \code{sppEquiv}.
#'   If these are not present, then the function will attempt to merge with
#'   LandR::sppEquivalencies_CA, so the column name(s) of \code{sppEquiv} must match
#'   column names in LandR::sppEquivalencies_CA.
#' @export
#' @importFrom data.table setDT setnames data.table
#' @return
#' A data.table with potentially all columns in LandR::sppEquivalencies_CA.
#'
sppEquivCheck <- function(sppEquiv, ensureColumns = NULL) {
  if (is.null(dim(sppEquiv))) {
    sppEquiv <- equivalentName(sppEquiv, LandR::sppEquivalencies_CA,
                               column = P(sim)$sppEquivCol, multi = TRUE)
    sppEquiv <- data.table(tmp = sppEquiv)
    setnames(sppEquiv, new = P(sim)$sppEquivCol)
  }
  if (is(sppEquiv, "data.frame")) {
    if (!is.data.table(sppEquiv)) setDT(sppEquiv)
    if (!all(c(ensureColumns) %in% colnames(sppEquiv) )) {
      if (all(colnames(sppEquiv) %in% colnames(LandR::sppEquivalencies_CA))) {
        sppEquiv <- LandR::sppEquivalencies_CA[
          sppEquiv, on = intersect(colnames(sppEquiv), colnames(LandR::sppEquivalencies_CA))]
      } else {
        stop("Please provide sppEquiv as a data.table with at least one column of species names, ",
             "that shares a column name and species naming convention as in LandR::sppEquivalencies_CA")
      }
    }
  } else {
    stop("sppEquiv must be a data.table")
  }
  sppEquiv[]

}
