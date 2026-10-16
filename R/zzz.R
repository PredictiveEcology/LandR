.onLoad <- function(libname, pkgname) {
  ## set options using the approach used by devtools
  opts <- options()
  opts.LandR <- LandROptions()

  toset <- !(names(opts.LandR) %in% names(opts))
  if (any(toset)) options(opts.LandR[toset])

  invisible()
}

#' The `LandR` package environment
#'
#' Environment used internally to store internal package objects and methods.
#'
#' @keywords internal
#' @rdname pkgEnv
.pkgEnv <- new.env(parent = emptyenv())
