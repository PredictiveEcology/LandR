#' `LandR` options
#'
#' These provide top-level settings for `LandR` functions and the `LandR` `SpaDES` modules.
#' To see defaults, run `LandROptions()`. See Details below.
#'
#' @details
#' Below are options that can be set with `options("LandR.xxx" = newValue)`, where `xxx` is
#' one of the values below, and `newValue` is a new value to give the option. Sometimes these
#' options can be placed in the user's `.Rprofile` file so they persist between sessions.
#'
#' \describe{
#'   \item{`assertions`}{
#'     Default: `TRUE`. If `TRUE`, additional code checks are run during function calls;
#'     see [assertions].
#'   }
#'   \item{`lccLeadingProportion`}{
#'     Default: `0.75`. The same threshold as `vegLeadingProportion` below, for
#'     [lccMapGenerator()], which works from land-cover legend codes. NTEMS/EOSD
#'     (Wulder & Nelson 2003) call a stand coniferous or broadleaf at 75% or more of total
#'     basal area, and mixed wood below that.
#'   }
#'   \item{`mergeHybridSpruce`}{
#'     Default: `"engelmann"`. Which species the hybrid white x Engelmann spruce
#'     (`Pice_eng_gla`) is merged into in the `sppEquiv` returned by [speciesInStudyArea()]:
#'     `"engelmann"` (`Pice_eng`), `"white"` (`Pice_gla`), or `NA` for no merging.
#'   }
#'   \item{`scanfiMirror`}{
#'     Default: `TRUE`. When `LandR` is loaded and no `reproducible.urlRemap` is set, fetch
#'     SCANFI from the PredictiveEcology arbutus mirror instead of Google Drive, so the ids
#'     that 404 for anonymous users resolve and no login is needed; see [scanfiUrlRemap()].
#'   }
#'   \item{`vegLeadingProportion`}{
#'     Default: `0.8`. The share of a stand held by one type, above which the stand stops
#'     being called mixed. Used by [vegTypeMapGenerator()], [vegTypeGenerator()] and
#'     [plotVTM()].
#'   }
#'   \item{`verbose`}{
#'     Default: `1`. The default `verbose` argument of functions that report their progress,
#'     e.g., [updateCohortData()] and [LANDISDisp()]. Higher numbers give more messages.
#'   }
#' }
#'
#' One option does not use the `LandR.` prefix, because it is shared with the NTEMS-derived
#' products the leading-species thresholds come from:
#'
#' \describe{
#'   \item{`NTEMS.mixedwoodProp`}{
#' Default: `NULL`, i.e., left unset. Every function that needs a purity threshold reads a
#' nested pair of options,
#' `getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", <default>))`,
#' so setting `NTEMS.mixedwoodProp` moves `vegLeadingProportion` and `lccLeadingProportion`
#' together, while leaving it unset keeps each of them on its own default. Its `NULL` default
#' here documents the option without setting it: `options()` ignores a `NULL`, so the
#' fallthrough still happens.
#'   }
#' }
#'
#' @return
#' A named list of the `LandR` options and their default values, as set by `.onLoad()` when
#' the package is loaded. Options the user has already set are not overwritten.
#'
#' @aliases opts.LandR
#' @export
#' @rdname LandROptions
LandROptions <- function() {
  list(
    LandR.assertions = TRUE,
    ## lccMapGenerator(): a land-cover legend code, so a different historical default than
    ## LandR.vegLeadingProportion. The two are the same purity threshold on a biomass-like
    ## share, and are kept apart so that adding these options changed no existing result.
    LandR.lccLeadingProportion = 0.75,
    LandR.mergeHybridSpruce = "engelmann",
    ## fetch SCANFI from the arbutus mirror when no reproducible.urlRemap is set; see
    ## ?scanfiUrlRemap
    LandR.scanfiMirror = TRUE,
    LandR.vegLeadingProportion = 0.8, ## vegTypeMapGenerator(), vegTypeGenerator(), plotVTM()
    LandR.verbose = 1,
    ## A NULL default documents the option without setting it: `options()` ignores a NULL, so
    ## the `getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", ...))`
    ## fallthrough still reaches the inner default. Setting it here would consume the outer
    ## slot and the fallthrough could never happen.
    NTEMS.mixedwoodProp = NULL
  )
}
