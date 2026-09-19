#' Prepare ecoregions objects
#'
#' @param ecoregionRst an optional raster object that could be passed to `sim`,
#'        representing ecoregions
#'
#' @param ecoregionLayer an `sf` polygons object representing ecoregions
#'
#' @param ecoregionLayerField optional. The field in `ecoregionLayer` that represents ecoregions.
#'
#' @template rasterToMatchLarge
#'
#' @param rstLCCAdj `RasterLayer` representing land cover adjusted for non-forest classes
#'
#' @param cacheTags `UserTags` to pass to cache
#'
#' @param pixelsToRm a vector of pixels to remove
#'
#' @export
prepEcoregions <- function(ecoregionRst = NULL, ecoregionLayer, ecoregionLayerField = NULL,
                           rasterToMatchLarge, rstLCCAdj, pixelsToRm, cacheTags) {
  appendEcoregionFactor <- FALSE ## whether or not to add the ecoregionClasses to the data

  if (is.null(ecoregionRst)) {
    ecoregionLayer <- fixErrors(ecoregionLayer)
    ecoregionMapSF <- sf::st_as_sf(ecoregionLayer) |>
      sf::st_transform(crs = sf::st_crs(rasterToMatchLarge))

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
      ecoregionMapSF[["ecoregionLayerField"]] <- as.factor(ecoDT$ecoregionLayerField)
      rm(ecoDT)
    }

    ## terra::rasterize creates a factor raster from a factor field, but uses "0" as the first value
    ## we will instead create integer field starting at 1.
    ecoregionMapSF$ecoregionLayerFieldInt <- as.integer(ecoregionMapSF$ecoregionLayerField)
    ecoregionRst <- terra::rasterize(ecoregionMapSF, rasterToMatchLarge, touches = TRUE,
                              field = "ecoregionLayerFieldInt")

    rm(ecoregionLayer)
    if (is.factor(ecoregionMapSF$ecoregionLayerField)) {
      appendEcoregionFactor <- TRUE
      ## Preserve factor values
      uniqVals <- unique(ecoregionMapSF$ecoregionLayerField)
      uniqIDs <- unique(ecoregionMapSF$ecoregionLayerFieldInt)
      df <- data.frame(
        ID = uniqIDs,
        ecoregionName = uniqVals,
        stringsAsFactors = FALSE
      )
      levels(ecoregionRst) <- df ## this will preserve the factors

      ecoregionTable <- as.data.table(df)
      ecoregionTable[, ID := as.factor(paddedFloatToChar(ID, max(nchar(ID))))]
    }
  } else {
    if (inherits(ecoregionRst, "RasterLayer")) {
      if (!length(ecoregionRst@data@attributes) == 0) {
        # Not sure this is what you intended. The is_empty was making the attribute table return empty
        appendEcoregionFactor <- TRUE
        ecoregionTable <- as.data.table(ecoregionRst@data@attributes[[1]])
        ecoregionTable[, ID := as.factor(paddedFloatToChar(ID, max(nchar(ID))))]
      }
    } else if (inherits(ecoregionRst, "SpatRaster")) {
      if (!is.null(levels(ecoregionRst))) {
        appendEcoregionFactor <- TRUE
        ecoregionTable <- as.data.table(levels(ecoregionRst))
        setnames(ecoregionTable, c("ID", "ecoregionName"))
        ecoregionTable[, ID := as.factor(ID)]
      }
    } else {
      stop("problem with ecoregionRst -- it is not a RasterLayer or a SpatRaster")
    }
  }

  if (!is.null(pixelsToRm)) {
    ecoregionRst[pixelsToRm] <- NA
  }

  message(cli::col_blue("Make initial ecoregionGroups ", Sys.time()))

  if (!isTRUE(.compareRas(ecoregionRst, rstLCCAdj, res = TRUE, stopOnError = FALSE))) {
    stop("problem with rasters ecoregionRst and rstLCCAdj -- they don't have same metadata")
  }
  ecoregionFiles <- ecoregionProducer(
    ecoregionMaps = list(ecoregionRst, rstLCCAdj),
    rasterToMatch = rasterToMatchLarge,
    ecoregionTable = ecoregionTable) |>
    Cache(
      # ecoregionTable = ecoregionTable,
      userTags = c(cacheTags, "ecoregionFiles", "stable"),
      omitArgs = c("userTags")
    )

  if (appendEcoregionFactor) {
    # These can be mismatched in the case where there are more factor levels in ecoregionTable than
    #   ecoregionFiles$ecoregion, but only if they have more digits, e.g., paddedFloatToChar may
    #   have created only 1:9 in ecoregionTable, but there may be 1:11 in ecoregionFiles
    maxNcharERT <- max(nchar(as.character(ecoregionTable[["ID"]])))
    erfEChar <- as.character(ecoregionFiles$ecoregion[["ecoregion"]])
    if (maxNcharERT != max(nchar(erfEChar))) {
      aa <- as.integer(erfEChar)
      aa <- factor(paddedFloatToChar(aa, padL = maxNcharERT))
      ecoregionFiles$ecoregion[["ecoregion"]] <- aa
    }
    ecoregionFiles$ecoregion <- ecoregionFiles$ecoregion[ecoregionTable, on = c("ecoregion" = "ID")] |>
      na.omit()
    setnames(ecoregionFiles$ecoregion, old = "ecoregion_lcc", new = "ecoregionGroup")
  }

  return(ecoregionFiles)
}
