utils::globalVariables(c(
  "active", "ecoregion", "ecoregion_lcc", "ecoregionGroup", "ID", "landcover", "mapcode"
))

#' Make \code{ecoregionMap} and \code{ecoregion} table
#'
#' This function combines an ecoregion map and a land cover map (e.g. ecodistricts and LCC)
#' and creates a map and table of containing their combined values and pixel IDs.
#' Used internally in LandR modules to prepare maps for to make \code{cohortData}.
#'
#' @param ecoregionMaps a \code{list} with two rasters, one with ecoregions (e.g. ecodistricts)
#' and another with land cover (e.g. LCC)
#' @param ecoregionName the name describing the type of ecoregions in first map
#' (e.g. \code{"ecoDistrict"}) if passing a polygon file
#' @template rasterToMatch
#'
#' @return
#' A list with two objects: the \code{ecoregionMap} and a table summarizing
#' its information per \code{pixelID}
#'
#' @export
#' @importFrom data.table as.data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom raster getValues levels raster
#' @importFrom sf st_as_sf
#' @importFrom reproducible paddedFloatToChar
ecoregionProducer <- function(ecoregionMaps, ecoregionName = NULL, rasterToMatch) {
  # change the coordinate reference for all spatialpolygons
  message("ecoregionProducer 1: ", Sys.time())
  #ecoregionMapInStudy <- raster::intersect(ecoregionMapFull, fixErrors(aggregate(studyArea)))

  # Alternative
  rstEcoregion <- list()
  rtmNAs <- is.na(rasterToMatch[]) | rasterToMatch[] == 0
  for (erm in seq(ecoregionMaps)) {
    if (!is(ecoregionMaps[[erm]], "Raster")) {
      message("ecoregionProducer fastRasterize: ", Sys.time())
      rstEcoregion[[erm]] <- fasterize::fasterize(sf::st_as_sf(ecoregionMaps[[erm]]),
                                                  raster(rasterToMatch),
                                                  field = ecoregionName)
    } else {
      rstEcoregion[[erm]] <- ecoregionMaps[[erm]]
    }
    rstEcoregion[[erm]][rtmNAs] <- NA
  }
  rstEcoregionNAs <- do.call(`|`, lapply(rstEcoregion, function(x) is.na(x[]) | x[] == 0))
  NAs <- rtmNAs | rstEcoregionNAs
  a <- lapply(rstEcoregion, function(x) getValues(x)[!NAs])
  b <- as.data.table(a)
  b[, (names(b)) := lapply(.SD, function(x) paddedFloatToChar(x, max(nchar(x), na.rm = TRUE)))]

  # Take the first 2 columns, whatever their names, in case they are given something
  ecoregionValues <- factor(paste(b[[1]], b[[2]], sep = "_"))

  rstEcoregion <- raster(rstEcoregion[[1]])
  ecoregionFactorLevels <- levels(ecoregionValues)

  rstEcoregion[!NAs] <- as.integer(ecoregionValues)
  levels(rstEcoregion) <- data.frame(ID = seq(ecoregionFactorLevels),
                                     mapcode = seq(ecoregionFactorLevels),
                                     ecoregion = gsub("_.*", "", ecoregionFactorLevels),
                                     landcover = gsub(".*_", "", ecoregionFactorLevels),
                                     ecoregion_lcc = ecoregionFactorLevels,
                                     stringsAsFactors = TRUE)

  ecoregionTable <- as.data.table(raster::levels(rstEcoregion)[[1]])
  message("ecoregionProducer mapvalues: ", Sys.time())
  ecoregionTable <- ecoregionTable[,.(active = "yes", mapcode, ecoregion, landcover, ecoregion_lcc)]

  return(list(ecoregionMap = rstEcoregion, ecoregion = ecoregionTable))
}

#' Make the \code{ecoregion} table
#'
#' This function creates a table containing pixel-wise ecoregion codes and whether they
#'   are "active" (have biomass > 0) or not for simulation. Unlike \code{ecoregionProducer},
#'   this function creates the \code{ecoregion} table from pixel information contained in
#'   \code{pixelCohortData}
#'
#' @template pixelCohortData
#' @template speciesEcoregion
#'
#' @return
#' A \code{data.table} with ecoregion codes and their active status per \code{pixelID}.
#'
#' @export
#' @importFrom data.table data.table
makeEcoregionDT <- function(pixelCohortData, speciesEcoregion) {
  ## make a table of available ecoregions
  ecoregion <- data.table(
    active = "yes",
    ecoregionGroup = factor(as.character(unique(pixelCohortData$ecoregionGroup)))
  )
  # Some ecoregions have NO BIOMASS -- so they are not active
  ecoregion[!ecoregionGroup %in% unique(speciesEcoregion$ecoregionGroup), active := "no"]

  return(ecoregion)
}

#' Make the \code{ecoregionMap} raster
#'
#' Creates a raster of ecoregion codes per pixel.
#' Unlike \code{ecoregionProducer}, this fills the raster with pixel information contained in
#' \code{pixelCohortData}.
#'
#' @template pixelCohortData
#' @param ecoregionFiles A list with two objects: the \code{ecoregionMap} and a table summarizing
#'   its information per \code{pixelID}.
#'
#' @return A raster with ecoregion codes.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom raster levels raster
makeEcoregionMap <- function(ecoregionFiles, pixelCohortData) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  ecoregionMap <-  raster(ecoregionFiles$ecoregionMap)

  # suppress this message call no non-missing arguments to min; returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(ecoregionMap[pixelData$pixelIndex] <- as.integer(pixelData$ecoregionGroup))
  levels(ecoregionMap) <- data.frame(ID = seq(levels(pixelData$ecoregionGroup)),
                                     ecoregion = gsub("_.*", "", levels(pixelData$ecoregionGroup)),
                                     ecoregionGroup = levels(pixelData$ecoregionGroup),
                                     stringsAsFactors = TRUE)
  return(ecoregionMap)
}
