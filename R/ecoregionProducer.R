if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(
    "mapcode"
  ))
}

#' Make \code{ecoregionMap} and \code{ecoregion} table
#'
#' This function combines an eco-region map and a land cover map (e.g. eco-districts and LCC)
#' and creates a map and table of containing their combined values and pixel IDs.
#' Used internally in LandR modules to prepare maps for to make \code{cohortData}.
#'
#' @param ecoregionMaps a \code{list} with two rasters, one with eco-regions (e.g. eco-districts) and
#' another with land cover (e.g. LCC)
#' @param ecoregionName the name describing the type of eco-regions in first map (e.g. "ecoDistrict")
#' @param ecoregionActiveStatus A two column \code{data.table} detailing with eco-regions are to be considered
#'  active for the simulations. Columns should be named 'active' (with 'yes' or 'no' values) and 'ecoregion'.
#' @param rasterToMatch a \code{rasterToMatch} (e.g. the one used throughout the simulation)
#'
#' @return
#' A list with two objects: the \code{ecoregionMap} and a table summarizing
#' it's information per pixelID
#'
#' @export
#' @importFrom data.table as.data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom raster getValues levels raster
#' @importFrom sf st_as_sf
#' @importFrom SpaDES.core paddedFloatToChar

ecoregionProducer <- function(ecoregionMaps, ecoregionName,
                              ecoregionActiveStatus, rasterToMatch) {
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
  a <- lapply(rstEcoregion, function(x) getValues(x)[!NAs] )
  b <- as.data.table(a)
  b[, (names(b)) := lapply(.SD, function(x) paddedFloatToChar(x, max(nchar(x), na.rm = TRUE)))]

  # Take the first 2 columns, whatever their names, in case they are given something
  ecoregionValues <- factor(paste(b[[1]], b[[2]], sep = "_"))

  rstEcoregion <- raster(rstEcoregion[[1]])
  ecoregionFactorLevels <- levels(ecoregionValues)

  rstEcoregion[!NAs] <- as.integer(ecoregionValues)
  levels(rstEcoregion) <- data.frame(ID = seq(ecoregionFactorLevels),
                                     mapcode = seq(ecoregionFactorLevels),
                                     ecoRegion = gsub("_.*", "", ecoregionFactorLevels),
                                     landCover = gsub(".*_", "", ecoregionFactorLevels),
                                     ecoregion = ecoregionFactorLevels)
  #ecoregionFactorValues <- na.omit(unique(rstEcoregion[]))

  #ecoregionTable <- raster::levels(rstEcoregion)[[1]]

  if (FALSE) {
    data.table(
      mapcode = seq_along(ecoregionFactorLevels[!is.na(ecoregionFactorLevels)]),
      ecoregion = ecoregionFactorLevels[!is.na(ecoregionFactorLevels)]
    )
    message("ecoregionProducer mapvalues: ", Sys.time())
    # rstEcoregion[] <- plyr::mapvalues(rstEcoregion[], from = ecoregionTable$ecoregion, to = ecoregionTable$mapcode)
    ecoregionActiveStatus[, ecoregion := as.factor(ecoregion)]
    ecoregionTable <- ecoregionTable[!is.na(mapcode),][, ecoregion := as.character(ecoregion)]
    message("ecoregionProducer dplyr_leftjoin: ", Sys.time())
    ecoregionTable <- dplyr::left_join(ecoregionTable,
                                       ecoregionActiveStatus,
                                       by = "ecoregion") %>%
      data.table()
    ecoregionTable[is.na(active), active := "no"]

    ecoregionTable <- ecoregionTable[,.(active, mapcode, ecoregion)]
  }

  ecoregionTable <- as.data.table(raster::levels(rstEcoregion)[[1]])
  message("ecoregionProducer mapvalues: ", Sys.time())
  ecoregionTable <- ecoregionTable[,.(active = "yes", mapcode, ecoregion)]

  return(list(ecoregionMap = rstEcoregion,
              ecoregion = ecoregionTable))
}
