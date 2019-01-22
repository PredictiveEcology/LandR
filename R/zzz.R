## be sure to update the 'Package Options' section of the package help file
##   in R/spades-core-package.R
##
.onLoad <- function(libname, pkgname) {
  ## set options using the approach used by devtools
  opts <- options()
  reproCachePath <- getOption("reproducible.cachePath")
  opts.LandR <- list( # nolint
    LandR.assertions = TRUE,
    LandR.verbose = 1
  )
  toset <- !(names(opts.LandR) %in% names(opts))
  if (any(toset)) options(opts.LandR[toset])

  invisible()
}
