#' @param speciesEcoregion A \code{data.table} with \code{species}-{ecoregion}-specific species
#' trait values. Ecoregion refers to "ecolocation", a categorical variable grouping sites with similar
#' biophysical characteristics. The table should have at least the following columns: \code{speciesCode} and
#' \code{ecoregionGroup}, character representation of species and ecoregion groups respectively,
#' \code{maxB} the maximum biomass for the species in a given 'ecoregion', \code{maxANPP} the maximum
#' aboveground net primary productivity and \code{SEP} the species establishment probability.
