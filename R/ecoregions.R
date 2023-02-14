utils::globalVariables(c(
  "active", "ecoregion", "ecoregion_lcc", "ecoregionGroup", "ID", "landcover", "mapcode"
))

#' Make `ecoregionMap` and `ecoregion` table
#'
#' This function combines an ecoregion map and a land cover map (e.g. ecodistricts and LCC)
#' and creates a map and table of containing their combined values and pixel IDs.
#' Used internally in LandR modules to prepare maps for to make `cohortData`.
#'
#' @param ecoregionMaps a `list` with two rasters, one with ecoregions (e.g. ecodistricts)
#' and another with land cover (e.g. LCC)
#' @param ecoregionName the name describing the type of ecoregions in first map
#' (e.g. `"ecoDistrict"`) if passing a polygon file
#' @template rasterToMatch
#'
#' @return
#' A list with two objects: the `ecoregionMap` and a table summarizing
#' its information per `pixelID`
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
  rtmNAs <- is.na(as.vector(rasterToMatch[])) | as.vector(rasterToMatch[]) == 0
  for (erm in seq(ecoregionMaps)) {
    if (!inherits(ecoregionMaps[[erm]], c("Raster", "SpatRaster"))) {
      message("ecoregionProducer fastRasterize: ", Sys.time())
      rstEcoregion[[erm]] <- fasterize::fasterize(sf::st_as_sf(ecoregionMaps[[erm]]),
                                                  raster(rasterToMatch),
                                                  field = ecoregionName)
    } else {
      rstEcoregion[[erm]] <- ecoregionMaps[[erm]]
    }
    rstEcoregion[[erm]][rtmNAs] <- NA
  }
  rstEcoregionNAs <- do.call(`|`, lapply(rstEcoregion, function(x) is.na(as.vector(x[])) | as.vector(x[]) == 0))
  NAs <- rtmNAs | rstEcoregionNAs
  a <- lapply(rstEcoregion, function(x) as.vector(x[])[!NAs])
  b <- as.data.table(a)
  b[, (names(b)) := lapply(.SD, function(x) paddedFloatToChar(x, max(nchar(x), na.rm = TRUE)))]

  # Take the first 2 columns, whatever their names, in case they are given something
  ecoregionValues <- factor(paste(b[[1]], b[[2]], sep = "_"))

  rstEcoregion <- rasterRead(rstEcoregion[[1]])
  ecoregionFactorLevels <- levels(ecoregionValues)

  rstEcoregion[!NAs] <- as.integer(ecoregionValues)
  levs <- data.frame(ID = seq(ecoregionFactorLevels),
                                     mapcode = seq(ecoregionFactorLevels),
                                     ecoregion = gsub("_.*", "", ecoregionFactorLevels),
                                     landcover = gsub(".*_", "", ecoregionFactorLevels),
                                     ecoregion_lcc = ecoregionFactorLevels,
                                     stringsAsFactors = TRUE)
  levels(rstEcoregion) <- levs

  ecoregionTable <- as.data.table(levs)
  message("ecoregionProducer mapvalues: ", Sys.time())
  ecoregionTable <- ecoregionTable[,.(active = "yes", mapcode, ecoregion, landcover, ecoregion_lcc)]

  return(list(ecoregionMap = rstEcoregion, ecoregion = ecoregionTable))
}

#' Make the `ecoregion` table
#'
#' This function creates a table containing pixel-wise ecoregion codes and whether they
#'   are "active" (have biomass > 0) or not for simulation. Unlike `ecoregionProducer`,
#'   this function creates the `ecoregion` table from pixel information contained in
#'   `pixelCohortData`
#'
#' @template pixelCohortData
#' @template speciesEcoregion
#'
#' @return
#' A `data.table` with ecoregion codes and their active status per `pixelID`.
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

#' Make the `ecoregionMap` raster
#'
#' Creates a raster of ecoregion codes per pixel.
#' Unlike `ecoregionProducer`, this fills the raster with pixel information contained in
#' `pixelCohortData`.
#'
#' @template pixelCohortData
#' @param ecoregionFiles A list with two objects: the `ecoregionMap` and a table summarizing
#'   its information per `pixelID`.
#'
#' @return A raster with ecoregion codes.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom raster levels raster
makeEcoregionMap <- function(ecoregionFiles, pixelCohortData) {
  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  ecoregionMap <-  eval(parse(text = getOption("reproducible.rasterRead", "terra::rast")))(ecoregionFiles$ecoregionMap)

  # suppress this message call no non-missing arguments to min; returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(ecoregionMap[pixelData$pixelIndex] <- as.integer(pixelData$ecoregionGroup))
  levels(ecoregionMap) <- data.frame(ID = seq(levels(pixelData$ecoregionGroup)),
                                     ecoregion = gsub("_.*", "", levels(pixelData$ecoregionGroup)),
                                     ecoregionGroup = levels(pixelData$ecoregionGroup),
                                     stringsAsFactors = TRUE)
  return(ecoregionMap)
}


#' Create Stacks of the speciesEcoregion content
#'
#' This will output a list of RasterStack objects. Each `RasterStack` show
#' raster maps of one of the columns listed in `columns` and each
#' `RasterLayer` will be one species.
#' @importFrom data.table data.table setDTthreads
#' @rawNamespace import(data.table, except = getNamespaceExports("data.table"))
#' @importFrom pemisc factorValues2
#' @importFrom raster stack raster
#' @template ecoregionMap
#' @template speciesEcoregion
#' @param columns The columns to use in the `speciesEcoregion` data.table.
#'   Default is `c("establishprob", "maxB", "maxANPP")`
speciesEcoregionStack <- function(ecoregionMap, speciesEcoregion,
                                  columns = c("establishprob", "maxB", "maxANPP")) {
  # stack of SEP
  # Require(c("data.table", "PredictiveEcology/pemisc", "raster"))
  # bm2011 <- biomassMaps2011
  # speciesEcoregion <- bm2011$speciesEcoregion
  orig <- data.table::setDTthreads(2)
  on.exit(data.table::setDTthreads(orig))
  whNonNAs <- which(!is.na(ecoregionMap[]))
  fv <- factorValues2(ecoregionMap,
                      ecoregionMap[][whNonNAs],
                      att = "ecoregionGroup")
  fvdt <- data.table(ecoregionGroup = as.character(fv), pixelID = whNonNAs)
  se2 <- fvdt[speciesEcoregion, on = "ecoregionGroup", allow.cartesian = TRUE]
  seList <- split(se2, by = "speciesCode")
  rasTemplate <- raster(ecoregionMap)
  names(columns) <- columns
  spp <- names(seList)
  stks <- lapply(columns, dtList = seList, rasTemplate = rasTemplate,
                 spp = spp,
                 function(column, dtList, rasTemplate, spp = spp) {
    createStack(dtList = dtList, rasTemplate = rasTemplate,
          column = column, spp = spp)
  })
}

createStack <- function(dtList, rasTemplate, column = "estblishprob", spp) {
  i <- 0
    # mclapply(mc.cores = length(dtList),
    outList <- lapply(
      dtList, rasTemplate = rasTemplate, column = column, spp = spp,
      function(dt, rasTemplate, column, spp) {
        i <<- i + 1
        print(paste(column, " ", spp[i]))
        rasTemplate[dt$pixelID] <- dt[[column]]
        print("... Done!")
        rasTemplate
      })
    raster::stack(outList)

}
