utils::globalVariables(c())

#' Create dummy inputs for test simulations
#'
#' `ecoregionMap`is a raster of all the unique groupings.
#'
#' @template rasterToMatch
#'
#' @return a `RasterLayer` object or, in the case of `makeDummyEcoregionFiles`, a list.
#'
#' @export
#' @rdname dummy-inputs
makeDummyEcoregionMap <- function(rasterToMatch) {
  ecoregionMap <- randomPolygons(
    ras = rasterToMatch,
    res = res(rasterToMatch),
    numTypes = 2
  )
  ecoregionMap <- mask(ecoregionMap, rasterToMatch)
  return(ecoregionMap)
}

#' @details
#' `rawBiomassMap` is a raster of "raw" total stand biomass per pixel,
#'      with values between 100 and 20000 g/m^2.
#'
#' @template rasterToMatch
#'
#' @export
#' @rdname dummy-inputs
makeDummyRawBiomassMap <- function(rasterToMatch) {
  rawBiomassMap <- neutralLandscapeMap(rasterToMatch,
    roughness = 0.65,
    rand_dev = 200,
    rescale = FALSE,
    verbose = FALSE
  )

  rawBiomassMap[] <- round(abs(as.vector(rawBiomassMap[])))
  return(rawBiomassMap)
}

#' @details
#' `standAgeMap` is a raster of stand age per pixel (where biomass exists)
#'      with values between 1 and 300 years.
#'
#' @param rawBiomassMap a `rawBiomassMap` (e.g. the one used
#'     throughout the simulation)
#'
#' @export
#' @rdname dummy-inputs
makeDummyStandAgeMap <- function(rawBiomassMap) {
  standAgeMap <- rawBiomassMap
  standAgeMap[] <- asInteger(rescale(as.vector(rawBiomassMap[]), c(1, 300)))
  return(standAgeMap)
}

#' @details
#' `rstLCC` is a raster land-cover class per pixel, with values between 1 and 5 that have no
#'      correspondence to any real land-cover classes.
#'
#' @template rasterToMatch
#'
#' @export
#' @rdname dummy-inputs
makeDummyRstLCC <- function(rasterToMatch) {
  rstLCC <- randomPolygons(
    ras = rasterToMatch,
    res = res(rasterToMatch),
    numTypes = 5
  )
  rstLCC <- mask(rstLCC, rasterToMatch)
  return(rstLCC)
}

#' @details
#' `ecoregionFiles` uses dummy versions of `ecoregionMap` and `rstLCC`
#' to create a list with two objects: the `ecoregionMap` and a table summarizing its
#' information per `pixelID`.
#' See `ecoregionProducer` (it uses `ecoregionProducer` internally).
#'
#' @template ecoregionMap
#' @template rstLCC
#' @template rasterToMatch
#'
#' @export
#' @rdname dummy-inputs
makeDummyEcoregionFiles <- function(ecoregionMap, rstLCC, rasterToMatch) {
  ecoregionstatus <- data.table(
    active = "yes",
    ecoregion = unique(as.vector(ecoregionMap[]))
  )
  ecoregionstatus <- ecoregionstatus[complete.cases(ecoregionstatus)]

  ecoregionFiles <- Cache(ecoregionProducer,
    ecoregionMaps = list(ecoregionMap, rstLCC),
    ecoregionName = "ECODISTRIC",
    rasterToMatch = rasterToMatch,
    userTags = "ecoregionFiles"
  )
  return(ecoregionFiles)
}

#' Rescale function (as in `scales::rescale`)
#'
#' This is a simple function copied from the `scales` package (almost the same).
#'
#' @param x a `numeric` vector
#' @param to a `numeric` vector of length 2. The new range of values.
#' @export
rescale <- function(x, to) {
  ## This is a simple function copied from the scales package.
  ## (package too heavy to use one simple function)
  from <- range(x, na.rm = TRUE, finite = TRUE)
  if (diff(range(to)) %==% 0 || diff(range(from)) %==% 0) {
    return(mean(to))
  }
  (x - from[1]) / diff(from) * diff(to) + to[1]
}
