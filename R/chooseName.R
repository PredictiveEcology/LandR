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
equivalentName <- function(value, df, column, multi = FALSE, searchColumn = NULL) {
  out <- lapply(df, function(x) {
    if (isTRUE(multi)) {
      which(x %in% as.character(value))
    } else {
      match(as.character(value), x)
    }
  })
  likelyMatch <- if (is.null(searchColumn)) {
    which.max(unlist(lapply(out, function(x) sum(!is.na(x)))))
  } else {
    searchColumn
  }
  df[[column]][out[[names(likelyMatch)]]]
}
