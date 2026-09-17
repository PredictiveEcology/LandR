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
#'   \item{`leadingSpeciesProp`}{
#'     Default: unset, which means it takes the value of `mixedwoodProp`. The share ONE
#'     species must hold for the stand to be called its own rather than mixed. A user who
#'     wants "just a majority" can set it to `0.51` without moving the mixedwood definition.
#'     Read with [leadingSpeciesProp()]; used by [vegTypeMapGenerator()] and
#'     [vegTypeGenerator()] for `mixedType` other than 2.
#'   }
#'   \item{`subsetDataSize`}{
#'     Default: `500L`. The maximum number of rows [subsetDT()] keeps per group when
#'     subsampling data for model fitting, and the default of the `LandR` modules'
#'     `subsetData*Model` parameters. Read with [subsetDataSize()]. It was 50 for years, set
#'     when these fits were expensive; at that size repeated runs of the same simulation gave
#'     visibly different `maxB`.
#'   }
#'   \item{`mergeHybridSpruce`}{
#'     Default: `"engelmann"`. Which species the hybrid white x Engelmann spruce
#'     (`Pice_eng_gla`) is merged into in the `sppEquiv` returned by [speciesInStudyArea()]:
#'     `"engelmann"` (`Pice_eng`), `"white"` (`Pice_gla`), or `NA` for no merging.
#'   }
#'   \item{`mixedwoodProp`}{
#'     Default: `0.75`. The share of a stand held by a GROUP -- all conifers, or all
#'     broadleaves -- at or above which the stand stops being mixedwood. This is the
#'     definition the national products use: NTEMS/EOSD (Wulder & Nelson 2003) and the NFI
#'     photo plot dictionary call a stand coniferous or broadleaf at 75% or more of total
#'     basal area (photo plots: total tree volume), mixed wood when neither group reaches it.
#'     Deciduous conifers count as conifers (`Larix` is `Type == "Conifer"` in
#'     [sppEquivalencies_CA]). Read with [mixedwoodProp()]; used by [vegTypeMapGenerator()]
#'     and [vegTypeGenerator()] with `mixedType = 2`, [lccMapGenerator()] and [plotVTM()].
#'   }
#'   \item{`verbose`}{
#'     Default: `1`. The default `verbose` argument of functions that report their progress,
#'     e.g., [updateCohortData()] and [LANDISDisp()]. Higher numbers give more messages.
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
    ## A NULL default documents the option without setting it: `options()` ignores a NULL, so
    ## leadingSpeciesProp() still falls through to LandR.mixedwoodProp. Setting it here would
    ## fix it at load, and moving LandR.mixedwoodProp afterwards would no longer move it.
    LandR.leadingSpeciesProp = NULL,
    LandR.mergeHybridSpruce = "engelmann",
    ## THE one place a leading/mixedwood threshold is written down. Everything reads it through
    ## mixedwoodProp() / leadingSpeciesProp(); no function carries its own default.
    LandR.mixedwoodProp = 0.75,
    ## THE one place the subsample size is written down; see subsetDataSize().
    LandR.subsetDataSize = 500L,
    LandR.verbose = 1
  )
}
