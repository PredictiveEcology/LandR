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
    if (!inherits(ecoregionRst, c("RasterLayer", "SpatRaster"))) {
      stop("problem with ecoregionRst -- it is not a RasterLayer or a SpatRaster")
    }
    ## A supplied categorical raster, from terra or raster, is normalized the same way: the
    ## raster's category table and `ecoregionTable` both become (ID, ecoregionName), with IDs
    ## padded as in the polygon branch above. ecoregionProducer() reads the labels from the
    ## raster's `ecoregionName` and joins them to `ecoregionTable` on that name, so the two must
    ## agree; the branches used to differ (the SpatRaster one renamed only the table and did not
    ## pad its IDs; the RasterLayer one renamed nothing).
    ecoregionCats <- .ecoregionCategories(ecoregionRst)
    if (!is.null(ecoregionCats)) {
      appendEcoregionFactor <- TRUE
      levels(ecoregionRst) <- if (inherits(ecoregionRst, "Raster")) list(ecoregionCats) else ecoregionCats
      ecoregionTable <- as.data.table(ecoregionCats)
      ecoregionTable[, ID := as.factor(paddedFloatToChar(ID, max(nchar(ID))))]
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

## (ID, ecoregionName) for a categorical ecoregion raster from terra or raster, taking the label
## column `.labelColumn()` picks; NULL when the raster is not categorical. (A non-categorical
## SpatRaster still has non-NULL `levels()`, so that test is not used.)
.ecoregionCategories <- function(r) {
  isCat <- if (inherits(r, "Raster")) raster::is.factor(r) else isTRUE(terra::is.factor(r)[1])
  if (!isCat) {
    return(NULL)
  }
  tb <- .categoryTable(r)
  data.frame(ID = tb[[1]], ecoregionName = as.character(tb[[.labelColumn(r, tb)]]), stringsAsFactors = FALSE)
}
