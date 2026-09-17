#' How many rows to keep per group when subsampling for model fitting
#'
#' [subsetDT()] keeps at most this many rows per group (e.g., per `ecoregionGroup` x
#' `speciesCode`) when fitting the statistical models, and the `LandR` modules use it as the
#' default of their `subsetData*Model` parameters. It is a function, and the number behind it
#' lives once in [LandROptions()], so a project can move it in one place with
#' `options(LandR.subsetDataSize = ...)` and every module follows.
#'
#' The default is 500. It was 50 for years, chosen when the fits were computationally
#' expensive; that subsample was small enough that repeated runs of the same simulation gave
#' visibly different `maxB` (a median coefficient of variation of 10% across
#' ecoregion x species, up to 62%, on a 60 km boreal test window). Ten times as many rows is
#' well within what current machines fit comfortably.
#'
#' @return A single number: the maximum rows per group to keep.
#'
#' @seealso [subsetDT()]
#'
#' @export
subsetDataSize <- function() {
  getOption("LandR.subsetDataSize", LandROptions()[["LandR.subsetDataSize"]])
}
