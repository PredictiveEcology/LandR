#' @param initialB the initial biomass of new cohorts. Defaults to 10.
#' If `NULL` or `NA`, initial cohort biomass is calculated as in LANDIS-II Biomass Succession
#' Extension v3.2.1 (Scheller & Miranda, 2015):
#'
#' ```r
#' initialB = asInteger(pmin(maxANPP, asInteger(pmax(1, maxANPP \* exp(-1.6 \* sumB / maxB_eco)))))
#' ```
#'
#' where `maxANPP` and `maxB_eco` are the maximum ANPP and B parameters of the species
#' in question within the pixel's `ecolocation`, and `sumB` is the total stand biomass
#' excluding cohorts with ages less than `successionTimestep`.
