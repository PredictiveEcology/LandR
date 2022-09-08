#' @param speciesTable The raw species traits table, with \emph{at least} the following columns:
#'  \code{speciesCode} an character representation of species, \code{Area} a character of
#'  geographical area within a species range for which trait values are relevant,
#'  \code{longevity} species longevity in years, \code{sexualmature} age in years at sexual maturity,
#'  \code{shadetolerance} numeric value between 1-5 of relative shade tolerance (with respect to oteher species),
#'  \code{seeddistance_eff} the "effective" seed dispersal distance, \code{seeddistance_max} a
#'  numeric with the maximum seed dispersal distance, and \code{mortalityshape} and \code{growthcureve} to growth
#'  curve shape parameters. Other columns (e.g. fire-related traits) can also be included depeding on LandR modules in use.
#'   Please see the LANDIS-II Biomass Succession Extension v3.2.1 manual (Scheller and Miranda 2015) for further detail.
