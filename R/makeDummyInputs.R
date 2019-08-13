if (getRversion() >= "3.1.0") {
  utils::globalVariables(c())
}


#' Makes a dummy \code{ecoregionMap}
#'
#' This is a function that makes a dummy \code{ecoregionMap} raster for
#' test simulations in LBMR.
#'
#' @param rasterToMatch a \code{rasterToMatch} (e.g. the one used
#'     throughout the simulation)
#'
#' @return
#' The \code{ecoregionMap}, a raster of all the unique groupings.
#'
#' @export
#' @importFrom raster mask res
#' @importFrom SpaDES.tools randomPolygons

makeDummyEcoregionMap <- function(rasterToMatch) {
  ecoregionMap <- randomPolygons(ras = rasterToMatch,
                                 res = res(rasterToMatch),
                                 numTypes = 2)
  ecoregionMap <- mask(ecoregionMap, rasterToMatch)
  return(ecoregionMap)
}


#' Makes a dummy \code{rawBiomassMap}
#'
#' This is a function that makes a dummy \code{rawBiomassMap} raster for
#' test simulations in LBMR.
#'
#' @param rasterToMatch a \code{rasterToMatch} (e.g. the one used
#'     throughout the simulation)
#'
#' @return
#' The \code{rawBiomassMap}, a raster of "raw" total stand biomass per pixel,
#'      with values between 100 and 20000 g/m2
#'
#' @export
#' @importFrom raster mask setValues getValues
#' @importFrom SpaDES.tools gaussMap
#' @importFrom scales rescale

makeDummyRawBiomassMap <- function(rasterToMatch) {
  rawBiomassMap <- gaussMap(rasterToMatch)
  rawBiomassMap <- setValues(rawBiomassMap,
                             rescale(getValues(rawBiomassMap), c(100, 20000)))
  rawBiomassMap <- mask(rawBiomassMap, rasterToMatch)
  return(rawBiomassMap)
}

#' Makes a dummy \code{standAgeMap}
#'
#' This is a function that makes a dummy \code{standAgeMap} raster for
#' test simulations in LBMR.
#'
#' @param rawBiomassMap a \code{rawBiomassMap} (e.g. the one used
#'     throughout the simulation)
#'
#' @return
#' The \code{standAgeMap}, a raster of stand age per pixel (where biomass exists)
#'      with values between 1 and 300 years
#'
#' @export
#' @importFrom raster setValues

makeDummyStandAgeMap <- function(rawBiomassMap) {
  standAgeMap <- setValues(rawBiomassMap, asInteger(rescale(getValues(biomassMap), c(1, 300))))
  return(standAgeMap)
}


#' Makes a dummy \code{rstLCC}
#'
#' This is a function that makes a dummy \code{rstLCC} raster for
#' test simulations in LBMR.
#'
#' @param rasterToMatch a \code{rasterToMatch} (e.g. the one used
#'     throughout the simulation)
#'
#' @return
#' The \code{rstLCC}, a raster land-cover class per pixel,
#'      with values between 1 and 5 that have no
#'      correspondence to any real land-cover classes.
#'
#' @export
#' @importFrom raster mask res
#' @importFrom SpaDES.tools randomPolygons

makeDummyRstLCC <- function(rasterToMatch) {
  rstLCC <- randomPolygons(ras = rasterToMatch,
                           res = res(rasterToMatch),
                           numTypes = 5)
  rstLCC <- mask(rstLCC, rasterToMatch)
  return(rstLCC)
}


#' Makes dummy \code{ecoregionFiles} list
#'
#' This is a function that makes a list of dummy \code{ecoregionFiles},
#'  using dummy versions of \code{ecoregionMap} and \code{rstLCC}
#'  (see \code{makeDummyEcoregionMap} and \code{makeDummyRstLCC})
#'  for test simulations in LBMR. It uses \code{ecoregionProducer} internally
#'
#' @param ecoregionMap a raster of all the unique groupings.
#' @param rstLCC a raster land-cover class per pixel,
#' @param rasterToMatch a \code{rasterToMatch} (e.g. the one used
#'     throughout the simulation)
#'
#' @return
#' The \code{ecoregionFiles}, A list with two objects: the \code{ecoregionMap}
#'   and a table summarizing its information per pixelID. See \code{ecoregionProducer}
#'
#' @export
#' @importFrom raster setValues
#' @importFrom data.table data.table

makeDummyEcoregionFiles <- function(ecoregionMap, rstLCC, rasterToMatch) {
  ecoregionstatus <- data.table(active = "yes",
                                ecoregion = unique(ecoregionMap[]))
  ecoregionstatus <- ecoregionstatus[complete.cases(ecoregionstatus)]

  ecoregionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(ecoregionMap, rstLCC),
                          ecoregionName = "ECODISTRIC",
                          ecoregionActiveStatus = ecoregionstatus,
                          rasterToMatch = rasterToMatch,
                          userTags = "ecoregionFiles")
  return(ecoregionFiles)
}




