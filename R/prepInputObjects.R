if (getRversion() >= "3.1.0") {
  utils::globalVariables(c())
}


#' Check if all species in have trait values
#'
#' @param speciesLayers stack of species layers rasters
#' @param species a \code{data.table} that has species traits such as longevity, shade tolerance, etc.
#' @param sppColorVect A named vector of colors to use for plotting. The names must conform with \codes{names(speciesLayers)}
#'   and should also contain a color for 'Mixed'.
#'
#' @return
#' A \code{list} with the \code{speciesLayers} and \code{sppColorVect}
#'   containing only the species that have trait values in \code{species}
#'
#' @export
#' @importFrom stats complete.cases
#' @importFrom crayon blue

checkSpeciesTraits <- function(speciesLayers, species, sppColorVect) {
  missTraits <- setdiff(names(speciesLayers), species$species)
  missTraits <- c(missTraits, setdiff(species$species,
                                      species[complete.cases(species), species]))
  if (length(missTraits)) {
    message(blue("The following species in 'speciesLayers' have missing traits",
                 "and will be excluded:\n", paste(missTraits, collapse = " "),
                 "\n If this is wrong check if species synonyms are included in 'sppEquiv'"))
    speciesLayers <- speciesLayers[[which(!names(speciesLayers) %in% missTraits)]]
    sppColorVect <- sppColorVect[c(names(speciesLayers), "Mixed")]
  }

  return(list(speciesLayers = speciesLayers, sppColorVect = sppColorVect))
}

#' Make pixelTable from biomass, age, land-cover and species cover data
#'
#' @param speciesLayers stack of species layers rasters
#' @param species a \code{data.table} that has species traits such as longevity, shade tolerance, etc.
#' @param standAgeMap raster of stand age
#' @param ecoregionFiles A list with two objects: the \code{ecoregionMap} and a table summarizing
#'   it's information per pixelID. See \code{ecoregionProducer}
#' @param biomassMap raster of total stand biomass
#' @param rasterToMatch raster mask of study area
#' @param LCC2005 raster of land-cover class
#' @param pixelGroupAgeClass When assigning pixelGroup membership, this defines the resolution
#'   of ages that will be considered 'the same pixelGroup', e.g., if it is 10, then 6 and 14
#'   will be the same
#'
#' #' @return
#' A \code{data.table} as many rows as non-NA pixels in \code{rasterToMath} and
#'  the columns containing pixel data from the input raster layers.
#'
#' @importFrom pemisc factorValues2
#' @importFrom crayon blue
#' @importFrom raster ncell
#' @importFrom data.table data.table

makePixelTable <- function(speciesLayers, species, standAgeMap, ecoregionFiles,
                           biomassMap, rasterToMatch, LCC2005, pixelGroupAgeClass) {
  message(blue("Round age to nearest P(sim)$pixelGroupAgeClass, which is",
               pixelGroupAgeClass))
  coverMatrix <- matrix(asInteger(speciesLayers[]),
                        ncol = length(names(speciesLayers)))
  colnames(coverMatrix) <- names(speciesLayers)

  pixelTable <- data.table(age = asInteger(ceiling(asInteger(standAgeMap[]) /
                                                     pixelGroupAgeClass) * pixelGroupAgeClass),
                           logAge = log(standAgeMap[]),
                           initialEcoregionCode = factor(factorValues2(ecoregionFiles$ecoregionMap,
                                                                       ecoregionFiles$ecoregionMap[],
                                                                       att = 5)),
                           totalBiomass = asInteger(biomassMap[] * 100), # change units
                           cover = coverMatrix,
                           pixelIndex = seq(ncell(rasterToMatch)),
                           lcc = LCC2005[],
                           rasterToMatch = rasterToMatch[])

  coverColNames <- paste0("cover.", species$species)
  pixelTable1 <- na.omit(pixelTable, cols = c("rasterToMatch"))
  pixelTable2 <- na.omit(pixelTable, cols = c("rasterToMatch", "initialEcoregionCode"))
  pixelTable <- na.omit(pixelTable2, cols = c(coverColNames))

  if (NROW(pixelTable1) != NROW(pixelTable))
    warning("Setting pixels to NA where there is NA in sim$speciesLayers'. Vegetation succession",
            "\n  parameters will only be calculated where there is data for species cover.",
            "\n  Check if sim$rasterToMatch shoudn't also only have data where there is cover data,",
            "\n  as this may affect other modules.")
  if (NROW(pixelTable2) != NROW(pixelTable))
    warning("Setting pixels to NA where there is NA in dummy 'ecoregionMap'")

  message(blue("rm NAs, leaving", magenta(NROW(pixelTable)), "pixels with data"))
  message(blue("This is the summary of the input data for age, ecoregionGroup, biomass, speciesLayers:"))
  print(summary(pixelTable))

  return(pixelTable)
}
