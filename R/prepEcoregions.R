#' Prepare ecoregions objects
#'
#' DESCRIPTION NEEDED
#'
#' @param ecoregionRst an optional raster object that could be passed to \code{sim}, representing ecoregions
#' @param ecoregionLayer a spatial polygons object representing ecoregions
#' @param ecoregionLayerField optional. The field in \code{ecoregionLayer} that represents ecoregions.
#' @param rasterToMatchLarge the \code{rasterToMatchLarge} object from \code{sim}
#' @param rstLCCAdj \code{RasterLayer} representing land cover adjusted for non-forest classes
#' @param cacheTags \code{UserTags} to pass to cache
#' @param pixelsToRm a vector of pixels to remove
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table data.table
#' @importFrom fasterize fasterize
#' @importFrom raster getValues levels raster
#' @importFrom reproducible Cache fixErrors
#' @importFrom SpaDES.core paddedFloatToChar
#' @importFrom sf st_as_sf
#' @export
prepEcoregions <- function(ecoregionRst = NULL, ecoregionLayer, ecoregionLayerField = NULL,
                           rasterToMatchLarge, rstLCCAdj, pixelsToRm, cacheTags) {

  appendEcoregionFactor <- FALSE ## whether or not to add the ecoregionClasses to the data

  if (is.null(ecoregionRst)) {
    ecoregionLayer <- fixErrors(ecoregionLayer)
    ecoregionMapSF <- sf::st_as_sf(ecoregionLayer)

    if (is.null(ecoregionLayerField)) {
      if (!is.null(ecoregionMapSF$ECODISTRIC)) {
        ecoregionMapSF$ecoregionLayerField <- as.factor(ecoregionMapSF$ECODISTRIC)
      } else {
        ecoregionMapSF$ecoregionLayerField <- as.numeric(row.names(ecoregionMapSF))
      }
    } else {
      ecoDT <- as.data.table(ecoregionMapSF)
      ecoregionField <- ecoregionLayerField
      ecoDT[, ecoregionLayerField := ecoDT[, get(ecoregionField)]]
      ecoregionMapSF$ecoregionLayerField <- as.factor(ecoDT$ecoregionLayerField)
      rm(ecoDT)
    }
    ecoregionRst <- fasterize::fasterize(ecoregionMapSF, raster = rasterToMatchLarge,
                                         field = "ecoregionLayerField")
    rm(ecoregionLayer)
    if (is.factor(ecoregionMapSF$ecoregionLayerField)) {
      appendEcoregionFactor <- TRUE
      #Preserve factor values
      uniqVals <- unique(ecoregionMapSF$ecoregionLayerField)
      df <- data.frame(ID = seq_len(length(uniqVals)),
                       ecoregionName = uniqVals,
                       stringsAsFactors = FALSE)
      levels(ecoregionRst) <- df #this will preserve the factors
    }
  } else {
    if (!length(ecoregionRst@data@attributes) == 0) {
      # Not sure this is what you intended. The is_empty was making the attribute table return empty
      appendEcoregionFactor <- TRUE
    }
  }

  ecoregionRst[pixelsToRm] <- NA

  message(blue("Make initial ecoregionGroups ", Sys.time()))

  assertthat::assert_that(isTRUE(compareRaster(ecoregionRst, rstLCCAdj,
                                               res = TRUE, orig = TRUE, stopiffalse = FALSE)))

  ecoregionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(ecoregionRst, rstLCCAdj),
                          rasterToMatch = rasterToMatchLarge,
                          userTags = c(cacheTags, "ecoregionFiles", "stable"),
                          omitArgs = c("userTags"))

  if (appendEcoregionFactor) {
    ecoregionTable <- as.data.table(ecoregionRst@data@attributes[[1]])
    ecoregionTable[, ID := as.factor(paddedFloatToChar(ID, max(nchar(ID))))]
    ecoregionFiles$ecoregion <- ecoregionFiles$ecoregion[ecoregionTable, on = c("ecoregion" = "ID")] %>%
      na.omit(.)
    setnames(ecoregionFiles$ecoregion, old = "ecoregion_lcc", new = "ecoregionGroup")
  }

  return(ecoregionFiles)
}
