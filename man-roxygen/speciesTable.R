#' @param speciesTable A species traits table, with **at least** the following columns:
#'  - `speciesCode` an character representation of species;
#'  - `Area` a character of geographical area within a species range for which trait values are relevant;
#'  - `longevity` species longevity in years, `sexualmature` age in years at sexual maturity;
#'  - `shadetolerance` numeric value between 1-5 of relative shade tolerance (with respect to other species);
#'  - `seeddistance_eff` the "effective" seed dispersal distance;
#'  - `seeddistance_max` a numeric with the maximum seed dispersal distance;
#'  - `mortalityshape` and `growthcurve`: growth curve shape parameters.
#'  Other columns (e.g. fire-related traits) can also be included depending on LandR modules in use.
#'  Please see the LANDIS-II Biomass Succession Extension v3.2.1 manual (Scheller and Miranda 2015)
#'  for further detail.
