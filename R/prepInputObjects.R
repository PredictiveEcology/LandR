utils::globalVariables(c(
  "cover", "ecoregionGroup", "establishprob", "lcc", "longevity", "maxB", "maxANPP",
  "postfireregen", "resproutprob", "speciesCode", "logAge"
))

#' Check if all species in have trait values
#'
#' @template speciesLayers
#' @template species
#' @template sppColorVect
#'
#' @return
#' A `list` with the `speciesLayers` and `sppColorVect`
#'   containing only the species that have trait values in `species`
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

#' Make `pixelTable` from biomass, age, land-cover and species cover data
#'
#' @template speciesLayers
#' @template standAgeMap
#' @param ecoregionFiles A list with two objects: the `ecoregionMap` and a table summarizing
#'   its information per `pixelID.` See `ecoregionProducer`.
#' @param biomassMap raster of total stand biomass
#' @template rasterToMatch
#' @template rstLCC
#' @param printSummary Logical. If `TRUE`, the default, a print out of the
#'   `summary(pixelTable)` will occur.
#' @template doAssertion
#'
#' @return
#' A `data.table` as many rows as non-NA pixels in `rasterToMath` and
#'  the columns containing pixel data from the input raster layers.
#'
#' @export
#' @importFrom crayon blue
#' @importFrom data.table data.table
#' @importFrom pemisc factorValues2
#' @importFrom raster ncell
makePixelTable <- function(speciesLayers, standAgeMap, ecoregionFiles,
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
  ## 3) For species that we have traits for
  coverColNames <- paste0("cover.", names(speciesLayers))
  pixelTable <- na.omit(pixelTable2, cols = c(coverColNames))

  if (NROW(pixelTable1) != NROW(pixelTable))
    message("Setting pixels to NA where there is NA in sim$speciesLayers. Vegetation succession",
            " parameters will only be calculated where there is data for species cover.",
            "\n  Check if rasterToMatch shoudn't also only have data where there is cover data,",
            " as this may affect other modules.")
  if (NROW(pixelTable2) != NROW(pixelTable))
    message("Setting pixels to NA where there is NA in 'ecoregionMap'")

  message(blue("rm NAs, leaving", magenta(NROW(pixelTable)), "pixels with data"))
  message(blue("This is the summary of the input data for age, ecoregionGroup, biomass, speciesLayers:"))
  if (isTRUE(printSummary)) print(summary(pixelTable))

  return(pixelTable)
}

#' Create `speciesEcoregion`
#'
#' Use statistically estimated `maxB`, `maxANPP` and establishment probabilities
#' to generate `specieEcoregion` table.
#'
#' See Details.
#'
#' @param cohortDataBiomass a subset of `cohortData`
#' @param cohortDataShort a subset of `cohortData`
#' @param cohortDataShortNoCover a subset of `cohortData`
#' @template species
#' @param modelCover statistical model of species presence/absence
#' @param modelBiomass statistical model of species biomass
#' @param successionTimestep The time between successive seed dispersal events.
#' @param currentYear `time(sim)`
#'
#' @section `establishprob`:
#' This section takes the cover as estimated from the mature tree cover and
#' partitions it between resprouting and seeds Unfortunately, establishment by
#' seed is not independent of resprouting, i.e., some pixels would have both
#' Since we don't know the level of independence, we can't correctly assess how
#' much to discount the two. If there is resprouting > 0, then this is the
#' partitioning:
#' `establishprob = f(establishprob + resproutprob + jointEstablishProbResproutProb)`
#' If `jointEstablishProbResproutProb` is 0, then these are independent events
#' and the total cover probability can be partitioned easily between seeds and
#' resprout. This is unlikely ever to be the case. We are picking 50% overlap as
#' a number that is better than 0 (totally independent probabilities, meaning no
#' pixel has both seeds and resprout potential) and  100% overlap (totally
#' dependent probabilities, i.e., every pixel where there is seeds will also be
#' a pixel with resprouting) This is expressed with the "* 0.5" in the code.
#'
#' #' @return
#' A `speciesEcoregion` `data.table` with added columns for parameters
#'   `maxB`, `maxANPP` and `establishprob`
#'
#' @export
#' @importFrom data.table rbindlist
makeSpeciesEcoregion <- function(cohortDataBiomass, cohortDataShort, cohortDataShortNoCover,
                                 species, modelCover, modelBiomass, successionTimestep, currentYear) {
  if (!is.null(modelBiomass$scaledVarsModelB)) {
    if (!is(modelBiomass$scaledVarsModelB, "list"))
      stop("modelBiomass$scaledVarsModelB must be a list")

    if (!all(names(modelBiomass$scaledVarsModelB) %in% c("cover", "logAge")))
      stop("modelBiomass$scaledVarsModelB must be a list with 'cover' and 'logAge' entries")
  }

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
    predict(modelCover$mod, newdata = cohortDataShort, type = "response")
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

  ## rescale if need be (modelBiomass may have been fitted on scaled variables)
  if (!is.null(modelBiomass$scaledVarsModelB)) {
    speciesEcoregion2 <- copy(speciesEcoregion)
    speciesEcoregion2[, `:=`(logAge = scale(logAge,
                                            center = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:center"),
                                            scale = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:scale")),
                             cover = scale(cover,
                                           center = attr(modelBiomass$scaledVarsModelB$cover, "scaled:center"),
                                           scale = attr(modelBiomass$scaledVarsModelB$cover, "scaled:scale")))]
    speciesEcoregion2[ , maxB := asInteger(predict(modelBiomass$mod,
                                                   newdata = speciesEcoregion2,
                                                   type = "response"))]
    speciesEcoregion[, maxB := speciesEcoregion2$maxB]
  } else {
    speciesEcoregion[ , maxB := asInteger(predict(modelBiomass$mod,
                                                  newdata = speciesEcoregion,
                                                  type = "response"))]
  }

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

#' Create `biomassMap`
#'
#' This is a function that creates the `biomassMap` raster used  for simulations in
#' `Biomass_core` module, using estimated data based on `rawBiomassMap` contained in
#' `pixelCohortData`.
#'
#' @template pixelCohortData
#' @template rasterToMatch
#'
#' @return The `biomassMap`, a raster of total stand biomass per pixel.
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

#' Create `minRelativeB` table
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
                             minRelativeBDefaults()
                             # X1 = 0.15, ## 0.2
                             # X2 = 0.25, ## 0.4
                             # X3 = 0.50, ## 0.5
                             # X4 = 0.75, ## 0.7
                             # X5 = 0.85  ## 0.9
  )

  return(minRelativeB)
}

#' minRelativeB defaults for Western Boreal Forest Canada
#'
#' @export
minRelativeBDefaults <- function() data.frame(X1 = 0.15, ## 0.2
                                              X2 = 0.25, ## 0.4
                                              X3 = 0.50, ## 0.5
                                              X4 = 0.75, ## 0.7
                                              X5 = 0.85)

#' Create `makePixelGroupMap`
#'
#' Create the `makePixelGroupMap` raster containing `pixelGroups` in `pixelCohortData`.
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

#' Create `standAgeMap`
#'
#' Create the `standAgeMap` raster containing age estimates for `pixelCohortData`.
#' A separate `prepInputs` call will source Canadian National Fire Data Base
#' data to update ages of recently burned pixels. To suppress this, pass NULL/NA `fireURL`
#'
#' @param ... additional arguments passed to `prepInputs`
#' @param ageURL url where age map is downloaded
#' @param ageFun passed to 'fun' arg of `prepInputs` of stand age map
#' @param maskWithRTM passed to `prepInputs` of stand age map
#' @param method passed to `prepInputs` of stand age map
#' @param datatype passed to `prepInputs` of stand age map
#' @param filename2 passed to `prepInputs` of stand age map
#' @param firePerimeters fire raster layer fire year values.
#' @param fireURL url to download fire polygons used to update age map. If NULL or NA age
#'   imputation is bypassed. Requires passing `rasterToMatch`. Only used if `firePerimeters`
#'   is missing.
#' @param fireFun passed to `prepInputs` of fire data. Only used if `firePerimeters`
#'   is missing.
#' @param fireField field used to rasterize fire polys. Only used if `firePerimeters`
#'   is missing.
#' @template destinationPath
#' @template rasterToMatch
#' @template startTime
#'
#' @return a raster layer stand age map corrected for fires, with an attribute vector of pixel IDs
#'  for which ages were corrected. If no corrections were applied the attribute vector is `integer(0)`.
#'
#' @export
#' @importFrom raster crs
#' @importFrom reproducible Cache prepInputs
#'
#' @examples
#' library(SpaDES.tools)
#' library(raster)
#' library(reproducible)
#' randomPoly <- randomStudyArea(size = 1e7)
#' randomPoly
#' ras2match <- raster(res = 250, ext = extent(randomPoly), crs = crs(randomPoly))
#' ras2match <- rasterize(randomPoly, ras2match)
#' tempDir <- tempdir()
#'
#' ## NOT USING FIRE PERIMETERS TO CORRECT STAND AGE
#' ## rasterToMatch does not need to be provided, but can be for masking/cropping.
#' standAge <- prepInputsStandAgeMap(destinationPath = tempDir,
#'                                   rasterToMatch = ras2match,
#'                                   fireURL = NA)   ## or NULL
#' attr(standAge, "imputedPixID")
#'
#' ## USING FIRE PERIMETERS TO CORRECT STAND AGE
#' ## ideaally, get the firePerimenters layer first
#' firePerimeters <- Cache(prepInputsFireYear,
#'                         url = paste0("https://cwfis.cfs.nrcan.gc.ca/downloads",
#'                         "/nfdb/fire_poly/current_version/NFDB_poly.zip"),
#'                         fun = "sf::st_read",
#'                         destinationPath = tempDir,
#'                         rasterToMatch = ras2match)
#'
#' standAge <- prepInputsStandAgeMap(destinationPath = tempDir,
#'                                   firePerimeters = firePerimeters,
#'                                   rasterToMatch = ras2match)
#' attr(standAge, "imputedPixID")
#'
#' ## not providing firePerimeters is still possible, but will be deprecated
#' ## in this case 'rasterToMatch' MUST be provided
#' standAge <- prepInputsStandAgeMap(destinationPath = tempDir,
#'                                   rasterToMatch = ras2match)
#' attr(standAge, "imputedPixID")
#'
prepInputsStandAgeMap <- function(..., ageURL = NULL,
                                  ageFun = "raster::raster",
                                  maskWithRTM = TRUE,
                                  method = "bilinear",
                                  datatype = "INT2U",
                                  destinationPath = NULL,
                                  filename2 = NULL,
                                  firePerimeters = NULL,
                                  fireURL = paste0("https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/",
                                                   "fire_poly/current_version/NFDB_poly.zip"),
                                  fireFun = "sf::st_read",
                                  fireField = "YEAR",
                                  rasterToMatch = NULL,
                                  startTime) {

  if (is.null(ageURL)) {
    ageURL <- paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                     "canada-forests-attributes_attributs-forests-canada/",
                     "2001-attributes_attributs-2001/",
                     "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")
  }

  getFires <- if (is.null(firePerimeters) &&
                  (isFALSE(is.null(fireURL)) && isFALSE(is.na(fireURL)))) {
    TRUE
  } else {
    FALSE
  }

  if (is.null(rasterToMatch)) {
    maskWithRTM <- FALSE
  }

  standAgeMap <- Cache(
    prepInputs, ...,
    maskWithRTM = maskWithRTM,
    method = method,
    datatype = datatype,
    filename2 = filename2,
    destinationPath = destinationPath,
    url = ageURL,
    fun = ageFun,
    rasterToMatch = rasterToMatch
  )
  standAgeMap[] <- asInteger(standAgeMap[])

  imputedPixID <- integer(0)
  if (getFires) {
    if (isFALSE(is.null(rasterToMatch))) {
      firePerimeters <- Cache(prepInputsFireYear, ...,
                              url = fireURL,
                              fun = fireFun,
                              fireField = fireField,
                              destinationPath = destinationPath,
                              rasterToMatch = rasterToMatch)
    } else {
      message("No 'rasterToMatch' or 'firePerimeters' supplied; ages will NOT be adjusted using fire data.")
    }
  }

  if (isFALSE(is.null(firePerimeters))) {
    standAgeMap <- replaceAgeInFires(standAgeMap, firePerimeters, startTime)
    imputedPixID <- attr(standAgeMap, "imputedPixID")
  }

  attr(standAgeMap, "imputedPixID") <- imputedPixID
  return(standAgeMap)
}

#' Create `rawBiomassMap`
#'
#' Create the `rawBiomassMap` raster containing biomass estimates for
#' `pixelCohortData`.
#' Wrapper on `prepInputs` that will rasterize fire polygons.
#'
#' @template studyAreaName
#' @template cacheTags
#' @param ... arguments passed to `prepInputs` and `Cache`. If the following arguments
#'   are not provided, the following values will be used:
#'   \itemize{
#'     \item{url: by default, the 2001 kNN stand biomass map is downloaded from
#'       the NRCan National Forest Inventory}
#'     \item{useSAcrs: FALSE}
#'     \item{method: "bilinear"}
#'     \item{datatype: "INT2U"}
#'     \item{filename2: `suffix("rawBiomassMap.tif", paste0("_", studyAreaName))`}
#'     \item{overwrite: TRUE}
#'     \item{userTags: `c(cacheTags, "rawBiomassMap")`}
#'     \item{omitArgs: `c("destinationPath", "targetFile", "userTags", "stable")`}
#'   }
#'
#' @return a rawBiomassMap raster
#' @export
#' @importFrom reproducible Cache prepInputs
prepRawBiomassMap <- function(studyAreaName, cacheTags, ...) {
  Args <- list(...)
  if (is.null(Args$url)) {
    Args$url <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                       "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
                       "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
  }
  if (is.null(Args$useSAcrs)) {
    Args$useSAcrs <- FALSE
  }
  if (is.null(Args$method)) {
    Args$method <- "bilinear"
  }
  if (is.null(Args$datatype)) {
    Args$datatype <- "INT2U"
  }
  if (is.null(Args$filename2)) {
    Args$filename2 <- .suffix("rawBiomassMap.tif", paste0("_", studyAreaName))
  }
  if (is.null(Args$overwrite)) {
    Args$overwrite <- TRUE
  }
  if (is.null(Args$userTags)) {
    Args$userTags <- c(cacheTags, "rawBiomassMap")
  }
  if (is.null(Args$omitArgs)) {
    Args$omitArgs <- c("destinationPath", "targetFile", "userTags", "stable")
  }
  Args$FUN <- prepInputs

  # httr::with_config(config = httr::config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
  #necessary for KNN
  rawBiomassMap <- do.call(Cache, args = Args)
  # })
  return(rawBiomassMap)
}

#' Create a raster of fire perimeters
#'
#' @param ... Additional arguments passed to `prepInputs`
#' @template rasterToMatch
#' @param fireField field used to rasterize fire polys
#' @param earliestYear the earliest fire date to allow
#'
#' @return a raster layer of fire perimeters with fire year values.
#'
#' @export
#' @importFrom fasterize fasterize
#' @importFrom raster crs
#' @importFrom reproducible Cache prepInputs
#' @importFrom sf st_cast st_transform
#' @importFrom magrittr %>%
#'
#' @examples
#' library(SpaDES.tools)
#' library(raster)
#' library(reproducible)
#' randomPoly <- randomStudyArea()
#' randomPoly
#' ras2match <- raster(res = 10, ext = extent(randomPoly), crs = crs(randomPoly))
#' ras2match <- rasterize(randomPoly, ras2match)
#' tempDir <- tempdir()
#' cacheRepo <- file.path(tempDir, "cache")
#'
#' ## ideally, get the firePerimenters layer first
#' firePerimeters <- Cache(prepInputsFireYear,
#'                         url = paste0("https://cwfis.cfs.nrcan.gc.ca/downloads",
#'                         "/nfdb/fire_poly/current_version/NFDB_poly.zip"),
#'                         fun = "sf::st_read",
#'                         destinationPath = tempDir,
#'                         rasterToMatch = ras2match)
#'
prepInputsFireYear <- function(..., rasterToMatch, fireField = "YEAR", earliestYear = 1950) {
  dots <- list(...)
  a <- if (is.null(dots$fun)) {
    Cache(prepInputs, rasterToMatch = rasterToMatch, ...) %>%
      st_as_sf(.)
  } else {
    if (grepl("st_read", dots$fun)) {
      Cache(prepInputs, ...)
    } else {
      Cache(prepInputs, rasterToMatch = rasterToMatch, ...) %>%
        st_as_sf(.)
    }
  }

  if (nrow(a) > 0) {
    gg <- st_cast(a, "MULTIPOLYGON") # collapse them into a single multipolygon
    d <- st_transform(gg, crs(rasterToMatch))
    if (!is(d[[fireField]], "numeric")) {
      warning("Chosen fireField will be coerced to numeric")
      d[[fireField]] <- as.numeric(as.factor(d[[fireField]]))
    }
    fireRas <- fasterize(d, raster = rasterToMatch, field = fireField)
    fireRas[!is.na(getValues(fireRas)) & getValues(fireRas) < earliestYear] <- NA
    return(fireRas)
  } else {
    return(NULL)
  }
}

#' Replace stand age with time since last fire
#'
#' @param standAgeMap a raster layer stand age map
#' @param firePerimeters the earliest fire date to allow
#' @template startTime
#'
#' @return a raster layer stand age map corrected for fires, with an attribute vector of pixel IDs
#'  for which ages were corrected. If no corrections were applied the attribute vector is `integer(0)`.
#'
#' @export
#' @importFrom fasterize fasterize
#' @importFrom raster crs
#' @importFrom reproducible Cache prepInputs
#' @importFrom sf st_cast st_transform
#' @importFrom magrittr %>%
#'
#' @examples
#' library(SpaDES.tools)
#' library(raster)
#' library(reproducible)
#' randomPoly <- randomStudyArea(size = 1e7)
#' randomPoly
#' ras2match <- raster(res = 250, ext = extent(randomPoly), crs = crs(randomPoly))
#' ras2match <- rasterize(randomPoly, ras2match)
#' tempDir <- tempdir()
#'
#' standAge <- prepInputsStandAgeMap(destinationPath = tempDir,
#'                                   rasterToMatch = ras2match,
#'                                   fireURL = NA)   ## or NULL
#' attr(standAge, "imputedPixID")
#'
#' firePerimeters <- Cache(prepInputsFireYear,
#'                         url = paste0("https://cwfis.cfs.nrcan.gc.ca/downloads",
#'                         "/nfdb/fire_poly/current_version/NFDB_poly.zip"),
#'                         fun = "sf::st_read",
#'                         destinationPath = tempDir,
#'                         rasterToMatch = ras2match)
#' standAge <- replaceAgeInFires(standAge, firePerimeters)
#' attr(standAge, "imputedPixID")
#'
replaceAgeInFires <- function(standAgeMap, firePerimeters, startTime) {
  if (missing(startTime)) {
    message("'startTime' is missing, the most recent fire year will be used.")
    startTime <- max(firePerimeters[], na.rm = TRUE)
  }

  toChange <- !is.na(firePerimeters[]) & firePerimeters[] <= asInteger(startTime)
  standAgeMap[] <- asInteger(standAgeMap[])
  standAgeMap[toChange] <- asInteger(startTime) - asInteger(firePerimeters[][toChange])
  imputedPixID <- which(toChange)

  attr(standAgeMap, "imputedPixID") <- imputedPixID
  return(standAgeMap)
}

#' Create `rasterToMatch` and `rasterToMatchLarge`
#'
#' `rasterToMatch` and `rasterToMatchLarge` raster layers are created
#'   from `studyArea` and `studyAreaLarge` polygons (respectively)
#'   using a template raster (often `rawBiomassMap`)
#'
#' @template studyArea
#' @param studyAreaLarge same as `studyArea`, but larger and completely
#'   covering it.
#' @template rasterToMatch
#' @template rasterToMatchLarge
#' @template destinationPath
#' @param templateRas a template raster used to make `rasterToMatch`
#'   and/or `rasterToMatchLarge`. Must match `studyAreaLarge`.
#' @template studyAreaName
#' @template cacheTags
#'
#' @export
#'
#' @importFrom reproducible Cache postProcessTerra fixErrors writeOutputs .suffix
#' @importFrom raster compareRaster
prepRasterToMatch <- function(studyArea, studyAreaLarge,
                              rasterToMatch, rasterToMatchLarge,
                              destinationPath,
                              templateRas, studyAreaName, cacheTags) {

  if (is.null(rasterToMatch) || is.null(rasterToMatchLarge)) {
    ## if we need rasterToMatch/rasterToMatchLarge, that means a) we don't have it,
    ## but b) we will have templateRas

    if (is.null(rasterToMatchLarge) && !is.null(rasterToMatch)) {
      rasterToMatchLarge <- rasterToMatch
    } else if (is.null(rasterToMatchLarge) && is.null(rasterToMatch)) {
      warning(paste0("rasterToMatch and rasterToMatchLarge are missing. Both will be created \n",
                     "from templateRas and studyArea/studyAreaLarge.\n
                     If this is wrong, provide both rasters"))

      if (is.null(templateRas)) {
        stop(paste("Please provide a template raster to make rasterToMatch(Large).",
                   "An option is to use 'rawBiomassMap'"))
      }
      if (!compareRaster(templateRas, studyAreaLarge, stopiffalse = FALSE)) {
        ## note that extents/origin may never align if the resolution and projection do not allow for it
        templateRas <- Cache(postProcessTerra,
                             templateRas,
                             studyArea = studyAreaLarge,
                             useSAcrs = TRUE,
                             overwrite = TRUE,
                             userTags = c("postRTMtemplate"))
        templateRas <- fixErrors(templateRas)
      }
      rasterToMatchLarge <- templateRas
    }

    if (!anyNA(rasterToMatchLarge[])) {
      whZeros <- rasterToMatchLarge[] == 0
      if (sum(whZeros) > 0) {# means there are zeros instead of NAs for RTML --> change
        rasterToMatchLarge[whZeros] <- NA
        message("There were no NAs on the rasterToMatchLarge, but there were zeros; converting these zeros to NA")
      }
    }

    RTMvals <- rasterToMatchLarge[]
    rasterToMatchLarge[!is.na(RTMvals)] <- 1

    rasterToMatchLarge <- Cache(
      writeOutputs,
      rasterToMatchLarge,
      filename2 = .suffix(file.path(destinationPath, "rasterToMatchLarge.tif"),
                          paste0("_", studyAreaName)),
      datatype = "INT2U",
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatchLarge"),
      omitArgs = c("userTags")
    )
    if (is.null(rasterToMatch)) {
      rasterToMatch <- Cache(postProcessTerra,
                             from = rasterToMatchLarge,
                             studyArea = studyArea,
                             # rasterToMatch = rasterToMatchLarge,   ## Ceres: this messes up the extent. if we are doing this it means BOTH RTMs come from biomassMap, so no need for RTMLarge here.
                             useSAcrs = FALSE,
                             # maskWithRTM = FALSE,   ## mask with SA
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = .suffix(file.path(destinationPath, "rasterToMatch.tif"),
                                                 paste0("_", studyAreaName)),
                             overwrite = TRUE,
                             # useCache = "overwrite",
                             userTags = c(cacheTags, "rasterToMatch"),
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable", "filename2",
                                          "overwrite"))
    }
    ## covert to 'mask'
    if (!anyNA(rasterToMatch[])) {
      whZeros <- rasterToMatch[] == 0
      if (sum(whZeros) > 0) {# means there are zeros instead of NAs for RTML --> change
        rasterToMatch[whZeros] <- NA
        message("There were no NAs on the RTM, but there were zeros; converting these zeros to NA")
      }
    }

    RTMvals <- rasterToMatch[]
    rasterToMatch[!is.na(RTMvals)] <- 1
  }

  return(list(rasterToMatch = rasterToMatch, rasterToMatchLarge = rasterToMatchLarge))
}

