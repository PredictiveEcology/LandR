utils::globalVariables(c(
  "cover", "ecoregionGroup", "establishprob", "lcc", "longevity", "maxB", "maxANPP",
  "postfireregen", "resproutprob", "speciesCode"
))

#' Check if all species in have trait values
#'
#' @param speciesLayers stack of species layers rasters
#' @param species a \code{data.table} with species traits such as longevity, shade tolerance, etc.
#' @param sppColorVect A named vector of colours to use for plotting.
#'                     The names must conform with \code{names(speciesLayers)} and should also
#'                     contain a colour for 'Mixed'.
#'
#' @return
#' A \code{list} with the \code{speciesLayers} and \code{sppColorVect}
#'   containing only the species that have trait values in \code{species}
#'
#' @export
#' @importFrom crayon blue
#' @importFrom stats complete.cases
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

#' Make \code{pixelTable} from biomass, age, land-cover and species cover data
#'
#' @param speciesLayers stack of species layers rasters
#' @param species a \code{data.table} with species traits such as longevity, shade tolerance, etc.
#' @param standAgeMap raster of stand age
#' @param ecoregionFiles A list with two objects: the \code{ecoregionMap} and a table summarizing
#'   its information per \code{pixelID.} See \code{ecoregionProducer}.
#' @param biomassMap raster of total stand biomass
#' @template rasterToMatch
#' @param rstLCC raster of land-cover class
#' @param printSummary Logical. If \code{TRUE}, the default, a print out of the
#'   \code{summary(pixelTable)} will occur.
#' @template doAssertion
#'
#' @return
#' A \code{data.table} as many rows as non-NA pixels in \code{rasterToMath} and
#'  the columns containing pixel data from the input raster layers.
#'
#' @export
#' @importFrom crayon blue
#' @importFrom data.table data.table
#' @importFrom pemisc factorValues2
#' @importFrom raster ncell
makePixelTable <- function(speciesLayers, species, standAgeMap, ecoregionFiles,
                           biomassMap, rasterToMatch, rstLCC, #pixelGroupAgeClass = 1,
                           printSummary = TRUE,
                           doAssertion = getOption("LandR.assertions", TRUE)) {
  if (missing(rasterToMatch)) {
    rasterToMatch <- raster(speciesLayers[[1]])
    rasterToMatch[] <- 0
    rasterToMatch[is.na(speciesLayers[[1]])] <- NA
  }

  if (missing(ecoregionFiles)) {
    ecoregionFiles <- list()
    ecoregionFiles$ecoregionMap <- raster(rasterToMatch)
    rtmNotNA <- which(!is.na(rasterToMatch[]))
    ecoregionFiles$ecoregionMap[rtmNotNA] <- seq_along(rtmNotNA)
    initialEcoregionCodeVals <- ecoregionFiles$ecoregionMap[]
  } else {
    initialEcoregionCodeVals <- factorValues2(
      ecoregionFiles$ecoregionMap,
      ecoregionFiles$ecoregionMap[],
      att = 5)
  }

  if (missing(species)) {
    species <- list()
    species$species <- names(speciesLayers)
  }

  # message(blue("Round age to nearest pixelGroupAgeClass, which is", pixelGroupAgeClass))
  coverMatrix <- matrix(asInteger(speciesLayers[]), ncol = length(names(speciesLayers)))
  colnames(coverMatrix) <- names(speciesLayers)

  # faster to use as.factor, which is fine for a numeric.
  iec <- if (is.numeric(initialEcoregionCodeVals)) {
    as.factor(initialEcoregionCodeVals)
  } else {
    factor(initialEcoregionCodeVals)
  }
  pixelTable <- data.table(initialEcoregionCode = iec,
                           cover = coverMatrix,
                           pixelIndex = seq(ncell(rasterToMatch)),
                           rasterToMatch = rasterToMatch[]
  )
  if (!missing(standAgeMap)) {
    set(pixelTable, NULL, "age", asInteger(standAgeMap[]))
    set(pixelTable, NULL, "logAge", .logFloor(standAgeMap[]))
  }

  if (!missing(biomassMap)) {
    set(pixelTable, NULL, "totalBiomass", asInteger(biomassMap[] * 100) ) # change units)
  }

  if (!missing(rstLCC)) {
    set(pixelTable, NULL, "lcc", rstLCC[])
  }

  #pixelTable <- data.table(#age = asInteger(ceiling(asInteger(standAgeMap[]) /
                          #                           pixelGroupAgeClass) * pixelGroupAgeClass),
                          # logAge = .logFloor(standAgeMap[]),
                          # initialEcoregionCode = factor(initialEcoregionCodeVals),
                          # totalBiomass = asInteger(biomassMap[] * 100), # change units
                          # cover = coverMatrix,
                          # pixelIndex = seq(ncell(rasterToMatch)),
                          # lcc = rstLCC[],
                          # rasterToMatch = rasterToMatch[])

  # Remove NAs from pixelTable
  ## 1) If in rasterToMatch
  pixelTable1 <- na.omit(pixelTable, cols = c("rasterToMatch"))
  ## 2) If in rasterToMatch and initialEcoregionCode
  pixelTable2 <- na.omit(pixelTable, cols = c("rasterToMatch", "initialEcoregionCode"))
  ## 3) For species that we have traits for (i.e., where species$species exists in the speciesLayers)
  coverColNames <- paste0("cover.", species$species)
  pixelTable <- na.omit(pixelTable2, cols = c(coverColNames))

  if (NROW(pixelTable1) != NROW(pixelTable))
    warning("Setting pixels to NA where there is NA in sim$speciesLayers'. Vegetation succession",
            " parameters will only be calculated where there is data for species cover.",
            "\n  Check if rasterToMatch shoudn't also only have data where there is cover data,",
            " as this may affect other modules.")
  if (NROW(pixelTable2) != NROW(pixelTable))
    warning("Setting pixels to NA where there is NA in dummy 'ecoregionMap'")

  message(blue("rm NAs, leaving", magenta(NROW(pixelTable)), "pixels with data"))
  message(blue("This is the summary of the input data for age, ecoregionGroup, biomass, speciesLayers:"))
  if (isTRUE(printSummary)) print(summary(pixelTable))

  return(pixelTable)
}

#' Create \code{speciesEcoregion}
#'
#' Use statistically estimated \code{maxB}, \code{maxANPP} and establishment probabilities
#' to generate \code{specieEcoregion} table.
#'
#' See Details.
#'
#' @param cohortDataBiomass a subset of \code{cohortData}
#' @param cohortDataShort a subset of \code{cohortData}
#' @param cohortDataShortNoCover a subset of \code{cohortData}
#' @param species a \code{data.table} of species traits, e.g., longevity, shade tolerance, etc.
#' @param modelCover statistical model of species presence/absence
#' @param modelBiomass statistical model of species biomass
#' @param successionTimestep The time between successive seed dispersal events.
#' @param currentYear \code{time(sim)}
#'
#' @section \code{establishprob}:
#' This section takes the cover as estimated from the mature tree cover and
#' partitions it between resprouting and seeds Unfortunately, establishment by
#' seed is not independent of resprouting, i.e., some pixels would have both
#' Since we don't know the level of independence, we can't correctly assess how
#' much to discount the two. If there is resprouting > 0, then this is the
#' partitioning:
#' \code{establishprob = f(establishprob + resproutprob + jointEstablishProbResproutProb)}
#' If \code{jointEstablishProbResproutProb} is 0, then these are independent events
#' and the total cover probability can be partitioned easily between seeds and
#' resprout. This is unlikely ever to be the case. We are picking 50% overlap as
#' a number that is better than 0 (totally independent probabilities, meaning no
#' pixel has both seeds and resprout potential) and  100% overlap (totally
#' dependent probabilities, i.e., every pixel where there is seeds will also be
#' a pixel with resprouting) This is expressed with the "* 0.5" in the code.
#'
#' #' @return
#' A \code{speciesEcoregion} \code{data.table} with added columns for parameters
#'   \code{maxB}, \code{maxANPP} and \code{establishprob}
#'
#' @export
#' @importFrom data.table rbindlist
makeSpeciesEcoregion <- function(cohortDataBiomass, cohortDataShort, cohortDataShortNoCover,
                                 species, modelCover, modelBiomass, successionTimestep, currentYear) {
  ## Create speciesEcoregion table
  joinOn <- c("ecoregionGroup", "speciesCode")
  speciesEcoregion <- unique(cohortDataBiomass, by = joinOn)
  speciesEcoregion[, c("B", "logAge", "cover") := NULL]
  species[, speciesCode := as.factor(species)]
  speciesEcoregion <- species[, .(speciesCode, longevity)][speciesEcoregion, on = "speciesCode"]
  speciesEcoregion[ , ecoregionGroup := factor(as.character(ecoregionGroup))]

  #################################################
  ## establishProb
  predictedCoverVals <- if (is(modelCover, "numeric")) {
    modelCover
  } else {
    predict(modelCover$mod, newdata = cohortDataShort,
                                type = "response")
  }
  establishprobBySuccessionTimestep <- 1 - (1 - predictedCoverVals)^successionTimestep
  cohortDataShort[, establishprob := establishprobBySuccessionTimestep]
  cohortDataShort <- species[, .(resproutprob, postfireregen, speciesCode)][cohortDataShort,
                                                                            on = "speciesCode"]

  # Partitioning between seed and resprout. See documentation about the "* 0.5"
  cohortDataShort[, establishprob := pmax(0, pmin(1, (establishprob * (1 - resproutprob * 0.5))))]

  cohortDataShort <- rbindlist(list(cohortDataShort, cohortDataShortNoCover),
                               use.names = TRUE, fill = TRUE)
  cohortDataShort[is.na(establishprob), establishprob := 0]

  # Join cohortDataShort with establishprob predictions to speciesEcoregion
  speciesEcoregion <- cohortDataShort[, .(ecoregionGroup, speciesCode, establishprob)][
    speciesEcoregion, on = joinOn]

  #################################################
  # maxB
  # Set age to the age of longevity and cover to 100%
  speciesEcoregion[, `:=`(logAge = .logFloor(longevity), cover = 100)]
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
  speciesEcoregion[ , `:=`(logAge = NULL, cover = NULL, longevity = NULL,  lcc = NULL)]

  speciesEcoregion[ , year := currentYear]
  return(speciesEcoregion)
}

#' Create \code{biomassMap}
#'
#' This is a function that creates the \code{biomassMap} raster used  for simulations in
#' \code{Biomass_core} module, using estimated data based on \code{rawBiomassMap} contained in
#' \code{pixelCohortData}.
#'
#' @template pixelCohortData
#' @template rasterToMatch
#'
#' @return The \code{biomassMap}, a raster of total stand biomass per pixel.
#'
#' @export
#' @importFrom raster raster
makeBiomassMap <-  function(pixelCohortData, rasterToMatch) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  biomassMap <- raster(rasterToMatch)
  # suppress this message call no non-missing arguments to min;
  # returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(biomassMap[pixelData$pixelIndex] <- pixelData$totalBiomass)

  return(biomassMap)
}

#' Create \code{minRelativeB} table
#'
#' The table contains expert-based values for minimum relative biomass of each shade tolerance
#' class (the minimum relative biomass a cohort with a given shade tolerance should have to be able
#' to germinate), in each unique ecoregion group.
#' All ecoregion groups currently have the same values.
#'
#' @template pixelCohortData
#'
#' @return a data.frame of min relative biomass values per ecoregion group.
#'
#' @export
makeMinRelativeB <- function(pixelCohortData) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  ## D. Cyr's values result in too many cohorts in more moisture-limited forests of Western Canada.
  ## https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/blob/master/LandisInputs/BSW/biomass-succession-main-inputs_BSW_Baseline.txt
  ##
  ## Adjusted values for western forests:
  minRelativeB <- data.frame(ecoregionGroup = as.factor(levels(pixelData$ecoregionGroup)),
                             X1 = 0.15, ## 0.2
                             X2 = 0.25, ## 0.4
                             X3 = 0.50, ## 0.5
                             X4 = 0.75, ## 0.7
                             X5 = 0.85) ## 0.9

  return(minRelativeB)
}

#' Create \code{makePixelGroupMap}
#'
#' Create the \code{makePixelGroupMap} raster containing \code{pixelGroups} in \code{pixelCohortData}.
#'
#' @template pixelCohortData
#' @template rasterToMatch
#'
#' @return a raster with pixel groups
#'
#' @export
#' @importFrom raster raster
makePixelGroupMap <- function(pixelCohortData, rasterToMatch) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  pixelGroupMap <- raster(rasterToMatch)

  ## suppress this message call no non-missing arguments to min;
  ## returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(pixelGroupMap[pixelData$pixelIndex] <- as.integer(pixelData$pixelGroup))

  return(pixelGroupMap)
}

#' Create \code{standAgeMap}
#'
#' Create the \code{standAgeMap} raster containing age estimates for \code{pixelCohortData}.
#' A separate \code{prepInputs} call will source NDFB data used to update ages of recently burned
#' pixels.
#'
#' @param ... additional arguments passed to \code{prepInputs}
#' @param ageURL url where age map is downloaded
#' @param ageFun passed to 'fun' arg of \code{prepInputs} of stand age map
#' @param maskWithRTM passed to \code{prepInputs} of stand age map
#' @param method passed to \code{prepInputs} of stand age map
#' @param datatype passed to \code{prepInputs} of stand age map
#' @param destinationPath directory where  age and fire data will be downloaded
#' @param filename2 passed to \code{prepInputs} of stand age map
#' @param fireURL url to download fire polygons used to update age map
#' @param fireFun passed to \code{prepInputs} of fire data
#' @template rasterToMatch
#' @param fireField field used to rasterize fire polys
#' @param startTime date of first fire year.
#'
#' @export
#' @importFrom httr with_config
#' @importFrom raster crs
#' @importFrom reproducible Cache prepInputs
prepInputsStandAgeMap <- function(..., ageURL = NULL,
                                  ageFun = "raster::raster",
                                  maskWithRTM = TRUE,
                                  method = "bilinear",
                                  datatype = "INT2U",
                                  destinationPath = NULL,
                                  filename2 = NULL,
                                  fireURL = NULL,
                                  fireFun = "sf::st_read",
                                  rasterToMatch, fireField = "YEAR",
                                  startTime) {
  if (is.null(ageURL))
    ageURL <- paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                     "canada-forests-attributes_attributs-forests-canada/",
                     "2001-attributes_attributs-2001/",
                     "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")

  if (is.null(fireURL))
    fireURL <- paste0("https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/",
                      "current_version/NFDB_poly.zip")

  with_config(config = config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
    standAgeMap <- Cache(prepInputs, ...,
                         maskWithRTM = maskWithRTM,
                         method = method,
                         datatype = datatype,
                         filename2 = filename2,
                         destinationPath = destinationPath,
                         url = ageURL,
                         fun = ageFun,
                         rasterToMatch = rasterToMatch)
  })
  standAgeMap[] <- asInteger(standAgeMap[])

  if (!(is.null(fireURL) || is.na(fireURL))) {
    fireYear <- Cache(prepInputsFireYear, ...,
                      url = fireURL,
                      fun = fireFun,
                      destinationPath = destinationPath,
                      rasterToMatch = rasterToMatch)

    toChange <- !is.na(fireYear[]) & fireYear[] <= asInteger(startTime)
    standAgeMap[] <- asInteger(standAgeMap[])
    standAgeMap[toChange] <- asInteger(startTime) - asInteger(fireYear[][toChange])
  }
  standAgeMap
}

#' Create a raster of fire polygons
#'
#' Wrapper on \code{prepInputs} that will rasterize fire polygons.
#'
#' @param ... Additional arguments passed to \code{prepInputs}
#' @template rasterToMatch
#' @param field field used to rasterize fire polys
#' @param earliestYear the earliest fire date to allow
#'
#' @export
#' @importFrom fasterize fasterize
#' @importFrom raster crs
#' @importFrom reproducible Cache prepInputs
#' @importFrom sf st_cast st_transform
prepInputsFireYear <- function(..., rasterToMatch, field = 'YEAR', earliestYear = 1950) {
  a <- Cache(prepInputs, ...)
  gg <- st_cast(a, "MULTIPOLYGON") # collapse them into a single multipolygon
  d <- st_transform(gg, crs(rasterToMatch))
  fireRas <- fasterize(d, raster = rasterToMatch, field = field)
  fireRas[!is.na(getValues(fireRas)) & getValues(fireRas) < earliestYear] <- NA
  return(fireRas)
}
