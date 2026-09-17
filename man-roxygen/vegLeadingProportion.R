#' @param vegLeadingProportion Numeric between 0-1, the share of a stand above which it stops
#'               being called mixed. Which share depends on what is being asked, and the two
#'               questions have their own options (see [leadingProportions]):
#'               with `mixedType = 2` it is the whole broadleaf group against the whole conifer
#'               group, and defaults to `mixedwoodProp()` (`LandR.mixedwoodProp`, the NTEMS/EOSD
#'               and NFI value: coniferous or broadleaf at 75% or more of total basal area,
#'               mixed wood when neither reaches it); otherwise it is one species' share and
#'               defaults to `leadingSpeciesProp()` (`LandR.leadingSpeciesProp`, which takes the
#'               mixedwood value unless it is set). Deciduous conifers count as conifers.
