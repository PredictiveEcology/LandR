#' @param newPixelCohortData must be a complete `cohortData` object with newly created
#'        cohorts. They do not have to have `pixelGroup` values yet;
#'        they can be overlapping with `cohortData`, (i.e., they can be regenerated on empty
#'        pixels or on already occupied pixels).
#'        Must contain the columns: `pixelIndex`, `speciesCode`, `ecoregionGroup`.
#'        The remaining 4 (see `cohortData`) will be created with `0`s.
