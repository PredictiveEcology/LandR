if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("resproutprob", "postfireregen", "speciesCode",
                           "establishprob", "ecoregionGroup", "speciesCode",
                           "longevity", "maxB", "maxANPP", "cover", "lcc"))
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


#' Make speciesEcoregion using statistically estimated maxB, maxANPP and establishment probabilities

#' @param cohortDataNoBiomass
#' @param cohortDataShort
#' @param cohortDataShortNoCover
#' @param species a \code{data.table} that has species traits such as longevity, shade tolerance, etc.
#' @param speciesEcoregion A \code{speciesEcoregion} table.
#' @param modelCover statistical model of species presence/absence
#' @param modelBiomass statistical model of species biomass
#' @param successionTimestep The time between successive seed dispersal events.
#' @param currentYear \code{time(sim)}
#'
#' #' @return
#' A \code{speciesEcoregion} table with added columns for parameters
#'   \code{maxB}, \code{maxANPP} and \code{establishprob}
#'
#' @importFrom data.table unique rbindlist

makeSpeciesEcoregion <- function(cohortDataNoBiomass, cohortDataShort, cohortDataShortNoCover, species,
                                 speciesEcoregion, modelCover, modelBiomass, successionTimestep, currentYear) {
  ## Create speciesEcoregion table
  joinOn <- c("ecoregionGroup", "speciesCode")
  speciesEcoregion <- unique(cohortDataNoBiomass, by = joinOn)
  speciesEcoregion[, c("B", "logAge", "cover") := NULL]
  species[, speciesCode := as.factor(species)]
  speciesEcoregion <- species[, .(speciesCode, longevity)][speciesEcoregion, on = "speciesCode"]
  speciesEcoregion[ , ecoregionGroup := factor(as.character(ecoregionGroup))]

  #################################################
  ## establishProb
  establishprobBySuccessionTimestep <- 1 - (1 - modelCover$pred)^successionTimestep
  cohortDataShort[, establishprob := establishprobBySuccessionTimestep]
  cohortDataShort <- species[, .(resproutprob, postfireregen, speciesCode)][cohortDataShort, on = "speciesCode"]
  cohortDataShort[, establishprob := pmax(0, pmin(1, (establishprob * (1 - resproutprob))))]

  cohortDataShort <- rbindlist(list(cohortDataShort, cohortDataShortNoCover),
                               use.names = TRUE, fill = TRUE)
  cohortDataShort[is.na(establishprob), establishprob := 0]

  # Join cohortDataShort with establishprob predictions to speciesEcoregion
  speciesEcoregion <- cohortDataShort[, .(ecoregionGroup, speciesCode, establishprob)][
    speciesEcoregion, on = joinOn]

  #################################################
  # maxB
  # Set age to the age of longevity and cover to 100%
  speciesEcoregion[, `:=`(logAge = log(longevity), cover = 100)]
  speciesEcoregion[ , maxB := asInteger(predict(modelBiomass$mod,
                                                newdata = speciesEcoregion,
                                                type = "response"))]
  speciesEcoregion[maxB < 0L, maxB := 0L] # fix negative predictions

  ########################################################################
  # maxANPP
  message(blue("Add maxANPP to speciesEcoregion -- currently --> maxB/30"))
  speciesEcoregion[ , maxANPP := asInteger(maxB / 30)]

  ########################################################################
  # Clean up unneeded columns
  speciesEcoregion[ , `:=`(logAge = NULL, cover = NULL, longevity = NULL, #pixelIndex = NULL,
                           lcc = NULL)]

  speciesEcoregion[ , year := currentYear]
  return(speciesEcoregion)
}


