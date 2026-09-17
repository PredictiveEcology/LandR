.onLoad <- function(libname, pkgname) {
  ## set options using the approach used by devtools
  opts <- options()
  reproCachePath <- getOption("reproducible.cachePath")
  opts.LandR <- list( #nolint
    LandR.assertions = TRUE,
    LandR.verbose = 1,
    ## fetch SCANFI from the arbutus mirror when no reproducible.urlRemap is set; see
    ## ?scanfiUrlRemap
    LandR.scanfiMirror = TRUE
  )
  toset <- !(names(opts.LandR) %in% names(opts))
  if (any(toset)) options(opts.LandR[toset])
  .setScanfiMirror()

  invisible()
}

#' The `LandR` package environment
#'
#' Environment used internally to store internal package objects and methods.
#'
#' @keywords internal
#' @rdname pkgEnv
.pkgEnv <- new.env(parent = emptyenv())
