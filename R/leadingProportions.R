#' Leading-species and mixedwood proportion thresholds
#'
#' Two proportions decide how a stand's composition is labelled, and they answer different
#' questions. Both are read from options so a project sets them once.
#'
#' `mixedwoodProp()` is the **group** threshold: the share of the stand held by all conifers, or
#' by all broadleaves, at or above which the stand stops being mixedwood. It is the definition
#' the national products are written on -- NTEMS/EOSD (Wulder & Nelson 2003) and the NFI photo
#' plot dictionary call a stand coniferous or broadleaf at 75% or more of total basal area (NFI
#' photo plots: total tree volume), and mixed wood when neither group reaches it.
#'
#' `leadingSpeciesProp()` is the **single-species** threshold: the share one species must hold to
#' be called leading. It is a different question and a user may reasonably want a different
#' answer -- 0.51, "just a majority", is a sensible setting -- so it has its own option. When it
#' is not set it takes the mixedwood value, which is what makes the two agree by default.
#'
#' Set them with, e.g., `options(LandR.mixedwoodProp = 0.7)` or
#' `options(LandR.leadingSpeciesProp = 0.51)`. The shipped default for `LandR.mixedwoodProp` lives
#' in [LandROptions()] and is the only place the number itself is written down, so changing it
#' later changes every function at once.
#'
#' Deciduous conifers are conifers here. `Larix` is `Type == "Conifer"` in
#' [sppEquivalencies_CA], so it never counts toward the broadleaf share.
#'
#' @return A numeric proportion between 0 and 1.
#'
#' @seealso [vegTypeMapGenerator()], [lccMapGenerator()]
#'
#' @export
#' @rdname leadingProportions
#' @aliases leadingProportions
mixedwoodProp <- function() {
  getOption("LandR.mixedwoodProp", LandROptions()[["LandR.mixedwoodProp"]])
}

#' @export
#' @rdname leadingProportions
leadingSpeciesProp <- function() {
  getOption("LandR.leadingSpeciesProp", mixedwoodProp())
}

## The threshold a caller that did not supply one should use, given what it is asking.
##
## `mixedType = 2` asks the mixedwood question (conifer group vs broadleaf group), so it takes
## `mixedwoodProp()`. Every other mixedType asks about a single species and takes
## `leadingSpeciesProp()`. Because `mixedType = 2` is the default in every caller, a user who
## sets `LandR.leadingSpeciesProp` and changes nothing else would silently see no effect; that
## is worth a warning rather than a surprise.
.leadingProp <- function(vegLeadingProportion = NULL, mixedType = 2) {
  if (!is.null(vegLeadingProportion)) {
    return(vegLeadingProportion)   ## an explicit argument always wins, and is never ambiguous
  }

  if (isTRUE(mixedType == 2)) {
    set <- getOption("LandR.leadingSpeciesProp")
    if (!is.null(set) && !isTRUE(all.equal(set, mixedwoodProp()))) {
      warning(
        "'LandR.leadingSpeciesProp' is set to ", set, ", but mixedType = 2 asks whether the ",
        "stand is mixedwood (all conifers vs all broadleaves), so it uses ",
        "'LandR.mixedwoodProp' (", mixedwoodProp(), ") instead. Set 'LandR.mixedwoodProp' to ",
        "change this call, or use mixedType = 1 for a single-species leading rule."
      )
    }
    mixedwoodProp()
  } else {
    leadingSpeciesProp()
  }
}

## Is this stand mixedwood? The broadleaf share is summed over the WHOLE broadleaf group and
## compared with the threshold; the conifer share is its complement. Called once per
## pixelGroup, so `speciesProportion` and `Type` are that group's rows.
##
## Deciduous conifers are conifers: `Type` comes from `sppEquiv`, where `Larix` is "Conifer",
## so tamarack biomass counts toward the conifer share and never makes a stand mixedwood.
## `Type %in% "Deciduous"` (not `==`) so a species missing from `sppEquiv`, whose Type is NA
## after the merge, is simply not counted as broadleaf rather than poisoning the sum.
##
## Previously this tested each deciduous species SEPARATELY -- a stand was mixedwood if any ONE
## deciduous species sat inside the band. Three broadleaf species at 0.15 each is 45% broadleaf,
## which is mixedwood by the group definition, but no single species reached the band, so the
## stand was called pure conifer.
.isMixedwood <- function(speciesProportion, Type, vegLeadingProportion) {
  broadleaf <- sum(speciesProportion[Type %in% "Deciduous"], na.rm = TRUE)
  broadleaf < vegLeadingProportion & broadleaf > 1 - vegLeadingProportion
}
