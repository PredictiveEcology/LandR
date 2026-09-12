#' Return equivalent name from a `data.frame` of equivalencies
#'
#' This is simply a wrapper around `match` or `\%in\%` for
#' a specific `data.frame` of values.
#'
#' @param value   Vector of values to match in `df`.
#' @param df      A `data.frame` where every row is a set of equivalent names.
#' @param column  A character string or numeric of length 1, indicating the column
#'                in `df` to return names from.
#' @param multi   Logical. If `TRUE`, then all matches will be returned.
#'                Default `FALSE` for backwards compatibility.
#'                This may result in more elements than were in `value`,
#'                as each `value` may be matched by more than one entry,
#'                returning more than one result each from `column`.
#'                If a species is found in the `df`, but there is no corresponding
#'                value in `column`, it will return `NA`.
#' @param searchColumn Optionally, provide the name of a column in `df` that results
#'                must be found in. The return value will still be from `column`
#'
#' @export
#' @rdname equivalentName
equivalentName <- function(value, df, column, multi = FALSE, searchColumn = NULL) {
  out <- equivalentNameAsList(value, df, multi)
  likelyMatch <- equivalentNameColumn(value, df, column, multi = multi, searchColumn = searchColumn)
  df[[column]][out[[likelyMatch]]]
}

#' `equivalentNameColumn` will provide the name of the column that best matches.
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

#' Check and expand `sppEquiv`
#'
#' This will expand a `sppEquiv` object that is only a vector or only a one-column
#' `data.table` into a many column `data.table`, if the columns that are present do not
#' contain `ensureColumns`.
#'
#' @param sppEquiv A character vector or `data.table` with named column(s).
#'   If this `data.table` does not have columns named `ensureColumns`,
#'   then it will attempt to merge this data.table with `sppEquivalencies_CA`
#'   to get `ensureColumns`.
#' @param ensureColumns A character vector of column names that must be in `sppEquiv`.
#'   If these are not present, then the function will attempt to merge with
#'   `sppEquivalencies_CA`, so the column name(s) of `sppEquiv` must match
#'   column names in `sppEquivalencies_CA`.
#' @param sppEquivCol Optional. Column in `sppEquivalencies_CA` to use for equivalent names
#'   when `sppEquiv` not provided (i.e., when `sppEquivalencies_CA` is used instead).
#'
#' @return A `data.table` with potentially all columns in `sppEquivalencies_CA`.
#'
#' @export
sppEquivCheck <- function(sppEquiv, ensureColumns = NULL, sppEquivCol = NULL) {
  sppEquivalencies_CA <- get(data("sppEquivalencies_CA", package = "LandR",
                                  envir = environment()), inherits = FALSE)
  if (is.null(dim(sppEquiv))) {
    stopifnot(!is.null(sppEquivCol))
    sppEquiv <- equivalentName(sppEquiv, sppEquivalencies_CA, column = sppEquivCol, multi = TRUE)
    sppEquiv <- data.table(tmp = sppEquiv)
    setnames(sppEquiv, new = sppEquivCol)
  }
  if (is(sppEquiv, "data.frame")) {
    if (!is.data.table(sppEquiv)) setDT(sppEquiv)
    if (!all(c(ensureColumns) %in% colnames(sppEquiv))) {
      if (all(colnames(sppEquiv) %in% colnames(sppEquivalencies_CA))) {
        sppEquiv <- sppEquivalencies_CA[
          sppEquiv,
          on = intersect(colnames(sppEquiv), colnames(sppEquivalencies_CA))
        ]
      } else {
        stop(
          "Please provide 'sppEquiv' as a data.table with at least one column of species names ",
          "that shares a column name and species naming convention as in LandR::sppEquivalencies_CA"
        )
      }
    }
  } else {
    stop("sppEquiv must be a data.table")
  }
  sppEquiv[]
}
