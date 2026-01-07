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
#' and another with land cover (e.g. LCC).
#'
#' @param ecoregionName the name describing the type of ecoregions in first map
#' (e.g. `"ecoDistrict"`) if passing a polygon file.
#' @param ecoregionTable A data.table that has 2 columns, `ecoregionName` (a factor)
#'   and `ID` a factor of `paddedFloatToChar(1:length(unique(ecoregionName)),
#'                         padL = max(nchar(length(unique(ecoregionName)))))`. This
#'   represents all the possible values that are available; this will be joined
#'   to `ecoregionMaps[[1]]` values
#'
#' @template rasterToMatch
#'
#' @return
#' A list with two objects: the `ecoregionMap` and a table summarizing
#' its information per `pixelID`
#'
#' @export
ecoregionProducer <- function(ecoregionMaps, ecoregionName = NULL, rasterToMatch,
                              ecoregionTable) {
  .requireNamespace("fasterize", stopOnFALSE = TRUE)

  ## change the coordinate reference for all spatialpolygons
  message("ecoregionProducer 1: ", Sys.time())
  # ecoregionMapInStudy <- intersect(ecoregionMapFull, fixErrors(aggregate(studyArea)))

  ## alternative
  rstEcoregion <- list()
  rtmNAs <- is.na(as.vector(rasterToMatch[])) | as.vector(rasterToMatch[]) == 0
  for (erm in seq(ecoregionMaps)) {
    if (!inherits(ecoregionMaps[[erm]], c("Raster", "SpatRaster"))) {
      message("ecoregionProducer fastRasterize: ", Sys.time())
      rstEcoregion[[erm]] <- fasterize::fasterize(sf::st_as_sf(ecoregionMaps[[erm]]),
        raster(rasterToMatch),
        field = ecoregionName
      )
    } else {
      rstEcoregion[[erm]] <- ecoregionMaps[[erm]]
    }
    rstEcoregion[[erm]][rtmNAs] <- NA
  }
  rstEcoregionNAs <- do.call(`|`, lapply(rstEcoregion, function(x) {
    is.na(as.vector(x[])) | as.vector(x[]) == 0
  }))
  NAs <- rtmNAs | rstEcoregionNAs
  # The next line fails when there are missing levels of the first one,
  #   if they don't have all the digits of the whole, i.e.,
  #   rstEcoregion[[1]] in one case has only levels 1:9
  #   but the ecoregionTable has 1:11, so the join outside this function
  #   fails
  # a <- lapply(rstEcoregion, function(x) as.vector(x[])[!NAs])
  # b[, (names(b)) := lapply(.SD, function(x) paddedFloatToChar(x, max(nchar(x), na.rm = TRUE)))]
  # New Jan 2, 2026 by Eliot
  a <- lapply(rstEcoregion, function(x) {
    if (is.factor(x)) raster::factorValues(x, values(x, mat = FALSE)[!NAs])
    else as.vector(x[])[!NAs]
    })
  b <- as.data.table(a)
  b <- ecoregionTable[b, on = "ecoregionName"] # join that gets all the correct values
  set(b, NULL, "ecoregionName", NULL)

  ## take the first 2 columns, whatever their names, in case they are given something
  ecoregionValues <- factor(paste(b[[1]], b[[2]], sep = "_"))

  rstEcoregion <- rasterRead(rstEcoregion[[1]])
  ecoregionFactorLevels <- levels(ecoregionValues)

  rstEcoregion[!NAs] <- as.integer(ecoregionValues)
  levs <- data.frame(
    ID = seq(ecoregionFactorLevels),
    mapcode = seq(ecoregionFactorLevels),
    ecoregion = gsub("_.*", "", ecoregionFactorLevels),
    landcover = gsub(".*_", "", ecoregionFactorLevels),
    ecoregion_lcc = ecoregionFactorLevels,
    stringsAsFactors = TRUE
  )
  levels(rstEcoregion) <- levs

  ecoregionTable <- as.data.table(levs)
  message("ecoregionProducer mapvalues: ", Sys.time())
  ecoregionTable <- ecoregionTable[, .(active = "yes", mapcode, ecoregion, landcover, ecoregion_lcc)]

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
makeEcoregionMap <- function(ecoregionFiles, pixelCohortData) {

  truePixelData <- as.data.table(ecoregionFiles$ecoregionMap, cells = TRUE)
  setnames(truePixelData, old = "cell", new = "pixelIndex")
  truePixelData[, mapcode := as.integer(mapcode)] #for join
  truePixelData <- truePixelData[ecoregionFiles$ecoregion, on = c("mapcode")]
  #keep only ecoregions for which we have data
  #but keep all observations of that ecoregion, regardless of whether it is currently filled
  truePixelData <- truePixelData[ecoregionGroup %in% pixelCohortData$ecoregionGroup]
  truePixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))]

  ecoregionMap <- rasterRead(ecoregionFiles$ecoregionMap)

  ## suppress this message call no non-missing arguments to min;
  ## returning Inf min(x@data@values, na.rm = TRUE)
  suppressWarnings(ecoregionMap[truePixelData$pixelIndex] <- as.integer(truePixelData$ecoregionGroup))

  factorDT <- unique(truePixelData[, .(ecoregionGroup, landcover, ecoregionName)])
  factorDT[, ID := seq(levels(ecoregionGroup))]
  factorDT[, ecoregion := gsub("_.*", "", ecoregionGroup)]
  setcolorder(factorDT, c("ID", "ecoregionGroup", "ecoregionName", "ecoregion", "landcover"))

  levels(ecoregionMap) <- factorDT

  return(ecoregionMap)
}

#' Create Stacks of the `speciesEcoregion` content
#'
#' Each `RasterStack` show raster maps of one of the columns listed in `columns`
#' and each `RasterLayer` will be one species.
#'
#' @template ecoregionMap
#'
#' @template speciesEcoregion
#'
#' @param columns The columns to use in the `speciesEcoregion` table.
#'                Default is `c("establishprob", "maxB", "maxANPP")`
#'
#' @returns list of `RasterStack` or `SpatRaster` objects
#'
speciesEcoregionStack <- function(ecoregionMap, speciesEcoregion,
                                  columns = c("establishprob", "maxB", "maxANPP")) {
  ## stack of SEP
  # Require(c("data.table", "PredictiveEcology/pemisc", "raster"))
  # bm2011 <- biomassMaps2011
  # speciesEcoregion <- bm2011$speciesEcoregion
  orig <- data.table::setDTthreads(2)
  on.exit(data.table::setDTthreads(orig), add = TRUE)
  whNonNAs <- which(!is.na(ecoregionMap[]))
  fv <- factorValues2(ecoregionMap,
    ecoregionMap[][whNonNAs],
    att = "ecoregionGroup"
  )
  fvdt <- data.table(ecoregionGroup = as.character(fv), pixelID = whNonNAs)
  se2 <- fvdt[speciesEcoregion, on = "ecoregionGroup", allow.cartesian = TRUE]
  seList <- split(se2, by = "speciesCode")
  rasTemplate <- rasterRead(ecoregionMap)
  names(columns) <- columns
  spp <- names(seList)
  stks <- lapply(columns,
    dtList = seList, rasTemplate = rasTemplate, spp = spp,
    function(column, dtList, rasTemplate, spp = spp) {
      createStack(dtList = dtList, rasTemplate = rasTemplate, column = column, spp = spp)
    }
  )
}

createStack <- function(dtList, rasTemplate, column = "estblishprob", spp) {
  i <- 0
  outList <- lapply(
    dtList,
    rasTemplate = rasTemplate, column = column, spp = spp,
      FUN = function(dt, rasTemplate, column, spp) {
        i <<- i + 1
        print(paste(column, " ", spp[i]))
        rasTemplate[dt$pixelID] <- dt[[column]]
        print("... Done!")
        rasTemplate
      }
  )

  .stack(outList)
}
