#' @param speciesEcoregion A `data.table` with `species`-`ecoregion`-specific species trait values.
#' Ecoregion refers to "ecolocation", a categorical variable grouping sites with similar biophysical
#' characteristics. The table should have at least the following columns: `speciesCode` and
#' `ecoregionGroup`, character representation of species and ecoregion groups respectively,
#' `maxB` the maximum biomass for the species in a given 'ecoregion', `maxANPP` the maximum
#' aboveground net primary productivity and `SEP` the species establishment probability.
#' May contain columns `inflationFactor`  (used to adjust `maxB`) and `mANPPproportion`
#' (used to calculate `maxANPP`).
