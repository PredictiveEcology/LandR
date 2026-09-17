#' @param vegLeadingProportion Numeric between 0-1, determining the relative biomass
#'               threshold a species needs to pass to be considered "leading".
#'               The default is a nested pair of options,
#'               `getOption("NTEMS.mixedwoodProp", getOption("LandR.vegLeadingProportion", 0.8))`
#'               -- set `NTEMS.mixedwoodProp` to move every such threshold at once, or the
#'               inner option to move only the vegetation-typing functions.
#'               [lccMapGenerator()] nests the same outer option over
#'               `LandR.lccLeadingProportion` (0.75), the NTEMS/EOSD value: Wulder & Nelson
#'               (2003) call a stand coniferous or broadleaf at 75% or more of total basal
#'               area, and mixed wood when neither reaches 75%.
