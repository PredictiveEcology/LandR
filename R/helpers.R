#' Convert a character string to sentence case
#'
#' @param x A character string.
#'
#' @export
toSentenceCase <- function(x) {
  newNames <- tolower(x)
  substr(newNames, 1, 1) <- toupper(substr(newNames, 1, 1))
  newNames
}

#' Determine name of current user
#'
#' @param username Character.
#'
#' @return Logical indicating whether the current user matches \code{username};
#' or, if \code{username} is NULL, a character string indicating the current user.
#'
#' @export
user <- function(username = NULL) {
  if (is.null(username)) {
    Sys.info()[["user"]]
  } else {
    identical(username, Sys.info()[["user"]])
  }
}
