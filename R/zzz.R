.onLoad <- function(libname, pkgname) {
  ## set options using the approach used by devtools
  opts <- options()
  reproCachePath <- getOption("reproducible.cachePath")
  opts.LandR <- list( #nolint
    LandR.assertions = TRUE,
    LandR.mergeHybridSpruce = "engelmann",
    LandR.verbose = 1,
    ## The share of a stand held by one type above which it stops being called mixed. Every
    ## function that needs it reads a nested pair:
    ##
    ##   getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", <default>))
    ##
    ## so `NTEMS.mixedwoodProp` is one knob that moves all of them together, while an unset
    ## `NTEMS.mixedwoodProp` leaves each family on the default it has always had. That is why
    ## `NTEMS.mixedwoodProp` is deliberately NOT set here: setting it would consume the outer
    ## slot and the fallthrough could never happen.
    ##
    ## The two inner defaults differ only by history, not by concept -- both are the same
    ## purity threshold on a biomass-like share. They are kept apart so that adding the
    ## options changes no existing result.
    LandR.vegLeadingProportion = 0.8,   ## vegTypeMapGenerator(), vegTypeGenerator(), plotVTM()
    ## lccMapGenerator(): a land-cover legend code. NTEMS/EOSD (Wulder & Nelson 2003) call a
    ## stand coniferous or broadleaf at 75% or more of total basal area, mixed wood below that.
    LandR.lccLeadingProportion = 0.75
  )
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
