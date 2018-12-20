if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", "B", "HQ", "leading", "LQ", "mixed",
                           "pixelGroup", "speciesCode", "speciesGroupB",
                           "speciesProportion","SPP", "totalB"))
}

#' Define flammability map
#'
#' @param LandCoverClassifiedMap A \code{Raster} that represents land cover
#' (e.g., Land Cover Classified map from 2005 or 2010 from the Canadian Forest Service).
#'
#' @param nonFlammClasses numeric vector defining which classes in \code{LandCoverClassifiedMap}.
#'
#' @param mask A raster to use as a mask (see \code{\link[raster]{mask}}).
#'
#' @param filename2 See \code{\link[reproducible]{postProcess}}. Default \code{NULL}.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom quickPlot setColors<-
#' @importFrom raster maxValue minValue ratify reclassify writeRaster
#' @export
defineFlammable <- function(LandCoverClassifiedMap = NULL,
                            nonFlammClasses = c(36L, 37L, 38L, 39L),
                            mask = NULL, filename2 = NULL) {
  if (!is.null(mask))
    if (!is(mask, "Raster")) stop("mask must be a raster layer")
  if (!is(LandCoverClassifiedMap, "RasterLayer")) {
    stop("Need a classified land cover map. Currently only accepts 'LCC2005'")
  }
  if (!is.integer(LandCoverClassifiedMap[]))
    stop("LandCoverCLassifiedMap must be an integer")
  if (is.null(nonFlammClasses))
    stop("Need nonFlammClasses, which are the classes that cannot burn in",
         "the LandCoverClassifiedMap")

  oldClass <- minValue(LandCoverClassifiedMap):maxValue(LandCoverClassifiedMap)
  newClass <- ifelse(oldClass %in% nonFlammClasses, 0L, 1L) ## NOTE: 0 codes for NON-flammable
  flammableTable <- cbind(oldClass, newClass)
  rstFlammable <- ratify(reclassify(LandCoverClassifiedMap, flammableTable))
  if (!is.null(filename2))
    rstFlammable <- writeRaster(rstFlammable, filename = filename2, overwrite = TRUE)

  setColors(rstFlammable, n = 2) <- colorRampPalette(c("blue", "red"))(2)
  if (!is.null(mask)) rstFlammable[is.na(mask[])] <- NA_integer_

  rstFlammable[] <- as.integer(rstFlammable[])
  rstFlammable
}

#' Simple prepInputs for LCC2005 or LCC2010
#'
#' A wrapper around prepInputs for the Canadian Land Cover Classification product(s)
#'
#' @inheritParams reproducible::cropInputs
#' @inheritParams reproducible::postProcess
#' @inheritParams reproducible::prepInputs
#'
#' @param year Numeric, either 2005 or 2010 (not yet implemented)
#'
#' @export
#' @importFrom reproducible asPath prepInputs
prepInputsLCC <- function(year = 2005,
                          destinationPath = asPath("."),
                          studyArea = NULL,
                          rasterToMatch = NULL,
                          filename2 = NULL, ...) {

  dots <- list(...)
  if (is.null(dots$url)) {
    if (identical(as.integer(year), 2005L)) {
      url <- "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip" #nolint
      filename <- asPath("LCC2005_V1_4a.tif")
    } else {
      stop("Don't have the url for LCC2010 yet. Plese pass it using 'url'.")
    }
  }

  Cache(prepInputs, targetFile = filename,
        archive = asPath("LandCoverOfCanada2005_V1_4.zip"),
        url = url,
        destinationPath = asPath(destinationPath),
        studyArea = studyArea,
        rasterToMatch = rasterToMatch,
        method = "bilinear",
        datatype = "INT2U",
        filename2 = filename2, ...)
}

#' Make a vegetation type map from a stack of species abundances
#'
#' @param speciesStack A \code{RasterStack} of species abundances.
#'                     This must be one \code{RasterLayer} per species.
#' @param vegLeadingProportion The threshold as a proportion of the total abundance
#'        that a species must have to be considered a "pure" stand of that type.
#'        If no species reaches this proportion, then the pixel will be 'Mixed'.
#' @param mixed Logical. If \code{TRUE}, then a mixed pixel value will be identified and given
#'        (see \code{vegLeadingProportion} argument)
#' @return A factor raster
#'
#' @export
#' @importFrom quickPlot numLayers
#' @importFrom raster levels maxValue raster
#' @importFrom SpaDES.tools inRange
makeVegTypeMap <- function(speciesStack, vegLeadingProportion, mixed = TRUE) {
  if (!inRange(vegLeadingProportion, 0, 1))
    stop("vegLeadingProportion must be a proportion")

  sumVegPct <- sum(speciesStack) ## TODO: how is the sum >100 ?


  if (isTRUE(mixed)) {
    ## create "mixed" layer, which is given a value slightly higher than any
    ## other layer. if it is deemed a mixed pixel
    ## All layers must be below vegLeadingProportion to be called Mixed
    ##  This check turns stack to binary: 1 if < vegLeadingProportion; 0 if more than.
    ##  Then, sum should be numLayers of all are below vegLeadingProportion
    whMixed <- which(sum(speciesStack < (100 * vegLeadingProportion))[] ==
                       numLayers(speciesStack))
    MixedRas <- speciesStack[[1]]
    MixedRas[!is.na(speciesStack[[1]][])] <- 0
    MixedRas[whMixed] <- max(maxValue(speciesStack)) * 1.01

    speciesStack$Mixed <- MixedRas
  }

  a <- speciesStack[]
  nas <- is.na(a[,1])
  maxes <- apply(a[!nas,], 1, function(x) {
    whMax <- which(x == max(x, na.rm = TRUE))
    if (length(whMax) > 1) {
      whMax <- sample(whMax, size = 1)
    }
    return(whMax)
  })

  vegTypeMap <- raster(speciesStack[[1]])

  vegTypeMap[!nas] <- maxes

  layerNames <- names(speciesStack)
  names(layerNames) <- layerNames
  levels(vegTypeMap) <- data.frame(ID = seq(layerNames), Species = names(layerNames))
  vegTypeMap
}

#' Generate vegetation type map
#'
#' @param cohortdata     A different \code{data.table} of stuff
#' @param pixelGroupMap  A \code{raster}
#' @param vegLeadingProportion Numeric between 0-1.
#' @param colors A named vector of color codes. The names MUST match the names of species
#'               in \code{cohortdata$speciesCode}, plus an optional "Mixed" color
#' @param unitTest A logical indicating whether some internal tests should be run to
#'                 ensure the function is running correctly.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom data.table data.table setkey
#' @importFrom raster getValues projection projection<- setValues
#' @importFrom SpaDES.tools rasterizeReduced
vegTypeMapGenerator <- function(cohortdata, pixelGroupMap, vegLeadingProportion,
                                colors, unitTest = getOption("pemisc.unitTest", FALSE)) {
  # shortcohortdata <- setkey(cohortdata, speciesCode)[setkey(species[, .(speciesCode, speciesGroup)],
  #                                                           speciesCode), nomatch = 0]
  pixelGroupData <- cohortdata[, list(totalB = sum(B, na.rm = TRUE),
                                      speciesCode, B), by = pixelGroup]

  pixelGroupData <- pixelGroupData[, .(speciesGroupB = sum(B, na.rm = TRUE),
                                         totalB = mean(totalB, na.rm = TRUE)),
                                     by = c("pixelGroup", "speciesCode")]
  pixelGroupData[, speciesProportion := speciesGroupB / totalB]

  pixelGroupData2 <- pixelGroupData[, list(mixed = all(speciesProportion < vegLeadingProportion),
                        leading = speciesCode[which.max(speciesProportion)]), by = "pixelGroup"]
  pixelGroupData2[mixed == TRUE, leading := "Mixed"]

  vegTypeMap <- rasterizeReduced(pixelGroupData2, pixelGroupMap, "leading", "pixelGroup")

  if (isTRUE(unitTest)) {
    # TEST THE MAP
    ids <- sample(which(!is.na(vegTypeMap[])), 100, replace = FALSE)
    pgs <- pixelGroupMap[][ids]
    dups <- duplicated(pgs)
    pgs <- pgs[!dups]
    leadingTest <- factorValues2(vegTypeMap, vegTypeMap[ids[!dups]], att = 2)
    names(pgs) <- leadingTest
    pgTest <- pixelGroupData[pixelGroup %in% pgs]
    pgTest2 <- pgTest[, list(mixed = all(speciesProportion < vegLeadingProportion),
                             leading = speciesCode[which.max(speciesProportion)]),
                      by = "pixelGroup"]
    out <- pgTest2[mixed == TRUE, leading := "Mixed"]
    length(unique(out$pixelGroup))
    length(pgs %in% unique(pgTest2$pixelGroup))
    pgs2 <- pgs[pgs %in% unique(pgTest2$pixelGroup)]
    if (!all(setkey(pgTest2, pixelGroup)$leading == names(pgs)[order(pgs)]))
      stop("The vegTypeMap is incorrect. Please debug pemisc::vegTypeMapGenerator")
  }

  # if (FALSE) {
  #
  #   species[species == "Pinu_ban" | species == "Pinu_con" | species == "Pinu_sp",
  #           speciesGroup := "PINU"]
  #   species[species == "Betu_pap" | species == "Popu_bal" | species == "Popu_tre" |
  #             species == "Lari_lar", speciesGroup := "DECI"]
  #   species[species == "Pice_mar", speciesGroup := "PICE_MAR"]
  #   species[species == "Pice_gla", speciesGroup := "PICE_GLA"]
  #   species[species == "Abie_sp", speciesGroup := "ABIE"]
  #
  #   shortcohortdata <- setkey(cohortdata, speciesCode)[setkey(species[, .(speciesCode, speciesGroup)],
  #                                                             speciesCode), nomatch = 0]
  #   shortcohortdata[, totalB := sum(B, na.rm = TRUE), by = pixelGroup]
  #   shortcohortdata <- shortcohortdata[, .(speciesGroupB = sum(B, na.rm = TRUE),
  #                                          totalB = mean(totalB, na.rm = TRUE)),
  #                                      by = c("pixelGroup", "speciesGroup")]
  #   shortcohortdata[, speciesProportion := speciesGroupB / totalB]
  #
  #   speciesLeading <- NULL
  #   Factor <- NULL #nolint
  #   ID <- NULL #nolint
  #   pixelGroup <- NULL
  #   speciesProportion <- NULL
  #   speciesGroup <- NULL
  #   speciesCode <- NULL
  #   totalB <- NULL
  #   B <- NULL #nolint
  #   speciesGroupB <- NULL
  #
  #   shortcohortdata[speciesGroup == "PINU" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 1] # pine leading
  #   shortcohortdata[speciesGroup == "DECI" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 2] # deciduous leading
  #   shortcohortdata[speciesGroup == "PICE_MAR" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 3] # black spruce leading
  #   shortcohortdata[speciesGroup == "PICE_GLA" & speciesProportion > vegLeadingProportion,
  #                   speciesLeading := 4] # white spruce leading
  #   shortcohortdata[is.na(speciesLeading), speciesLeading := 0]
  #   shortcohortdata[, speciesLeading := max(speciesLeading, na.rm = TRUE), by = pixelGroup]
  #   shortcohortdata <- unique(shortcohortdata[, .(pixelGroup, speciesLeading)], by = "pixelGroup")
  #   shortcohortdata[speciesLeading == 0, speciesLeading := 5] # 5 is mixed forests
  #   attritable <- data.table(ID = sort(unique(shortcohortdata$speciesLeading)))
  #   attritable[ID == 1, Factor := "Pine leading"]          #nolint
  #   attritable[ID == 2, Factor := "Deciduous leading"]     #nolint
  #   attritable[ID == 3, Factor := "Black spruce leading"]  #nolint
  #   attritable[ID == 4, Factor := "White spruce leading"]  #nolint
  #   attritable[ID == 5, Factor := "Mixed"]
  #   vegTypeMap <- rasterizeReduced(shortcohortdata, pixelGroupMap, "speciesLeading", "pixelGroup")
  #   vegTypeMap <- setValues(vegTypeMap, as.integer(getValues(vegTypeMap)))
  #   levels(vegTypeMap) <- as.data.frame(attritable)
  #   projection(vegTypeMap) <- projection(pixelGroupMap)
  # }
  vegTypeMap
}

#' Load kNN species layers from online data repository
#'
#' TODO: description needed
#'
#' @param dPath path to the data directory
#'
#' @param rasterToMatch passed to \code{\link[reproducible]{prepInputs}}
#'
#' @param studyArea passed to \code{\link[reproducible]{prepInputs}}
#'
#' @param sppEquiv table with species name equivalencies between the
#'                           kNN format and the final naming format.
#'                           See \code{data("sppEquivalencies_CA", "pemisc")}.
#'
#' @param knnNamesCol character string indicating the column in \code{sppEquiv}
#'                    containing kNN species names.
#'                    Default \code{"KNN"} for when \code{sppEquivalencies_CA} is used.
#'
#' @param sppEquivCol character string indicating the column in \code{sppEquiv}
#'                    to use for final species names.
#'
#' @param thresh the minimum number of pixels where the species must have
#'               \code{biomass > 0} to be considered present in the study area.
#'               Defaults to 1.
#'
#' @param url the source url for the data, passed to \code{\link[reproducible]{prepInputs}}
#'
#' @param ... Additonal arguments passed to \code{\link[reproducible]{Cache}}
#'            and \code{\link{equivalentName}}.
#'
#' @return A raster stack of percent cover layers by species.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom raster ncell raster
#' @importFrom reproducible Cache .prefix preProcess
#' @importFrom utils untar
loadkNNSpeciesLayers <- function(dPath, rasterToMatch, studyArea, sppEquiv,
                                 knnNamesCol = "KNN", sppEquivCol, thresh = 1, url, ...) {
  dots <- list(...)

  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]),]
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  sppMerge <- unique(sppEquiv[[sppEquivCol]][duplicated(sppEquiv[[sppEquivCol]])])

  if ("cachePath" %in% names(dots)) {
    cachePath <- dots$cachePath
  } else {
    cachePath <- getOption("reproducible.cachePath")
  }

  ## get .tar file first - no extraction
  outPreProcess <- preProcess(targetFile = file.path(dPath, "kNN-Species.tar"),
                              archive = file.path(dPath, "kNN-Species.tar"),
                              url = url, destinationPath = dPath)

  ## get all kNN species - names only
  allSpp <- untar(tarfile = outPreProcess$targetFilePath, list = TRUE)
  allSpp <- allSpp %>%
    grep(".zip", ., value = TRUE) %>%
    sub("_v0.zip", "", .) %>%
    sub(".*Species_", "", .)

  ## get all species layers from .tar
  if (length(sppNameVector) == 1) ## avoids a warning in next if
    if (sppNameVector == "all")
      sppNameVector <- allSpp

  ## Make sure spp names are compatible with kNN names
  kNNnames <- as.character(equivalentName(sppNameVector, sppEquiv, column = knnNamesCol, multi = TRUE))

  ## if there are NA's, that means some species can't be found in kNN data base
  if (any(is.na(kNNnames))) {
    warning(paste0("Can't find ", sppNameVector[is.na(kNNnames)], " in `sppEquiv$",
                   knnNamesCol, ".\n Will use remaining matching species, but check if this is correct"))
    ## select only available species
    kNNnames <- kNNnames[!is.na(kNNnames)]
    sppNameVector <- sppNameVector[!is.na(kNNnames)]
  }

  ## same as above
  if (any(!kNNnames %in% allSpp)) {
    warning(paste0("Can't find ", sppNameVector[is.na(kNNnames)], " in kNN database.\n",
                   "Will use remaining matching species, but check if this is correct."))
    kNNnames <- kNNnames[kNNnames %in% allSpp]
    sppNameVector <- sppNameVector[kNNnames %in% allSpp]
  }

  ## define suffix to append to file names
  suffix <- if (basename(cachePath) == "cache") {
    paste0(as.character(ncell(rasterToMatch)), "px")
  } else {
    basename(cachePath)
  }
  suffix <- paste0("_", suffix)

  ## select which archives/targetFiles to extract
  targetFiles <- paste0("NFI_MODIS250m_kNN_Species_", kNNnames, "_v0.tif")
  # this sorting is necessary because split fn call below will sort
  #   regardless, so this will ensure archive and targetFiles are same
  targetFiles <- sort(targetFiles)
  names(targetFiles) <- targetFiles
  archives <- cbind(archive1 = file.path(dPath, "kNN-Species.tar"),
                    archive2 = paste0("NFI_MODIS250m_kNN_Species_", kNNnames, "_v0.zip"))
  archives <- split(archives, archives[, "archive2"])

  if (!all(file_path_sans_ext(names(archives))==file_path_sans_ext(targetFiles)))
    stop("Something is wrong. File a bug report re: loadkNNSpeciesLayers and archive ordering")

  postProcessedFilenames <- .suffix(targetFiles, suffix = suffix)

  message("Running prepInputs for ", paste(kNNnames, collapse = ", "))
  speciesLayers <- Cache(Map, targetFile = targetFiles, archive = archives,
                       filename2 = postProcessedFilenames,
                       MoreArgs = list(url = url,
                                       destinationPath = asPath(dPath),
                                       fun = "raster::raster",
                                       studyArea = studyArea,
                                       rasterToMatch = rasterToMatch,
                                       method = "bilinear",
                                       datatype = "INT2U",
                                       overwrite = TRUE,
                                       userTags = dots$userTags
                       ),
                       prepInputs)

  names(speciesLayers) <- kNNnames

  if (!is.null(sppMerge)) {
    speciesLayers <- mergeSppRaster(sppMerge = sppMerge, speciesLayers = speciesLayers,
                                    sppEquiv = sppEquiv, column = "KNN", suffix = suffix,
                                    dPath = dPath)
  }

  ## Rename species layers - There will be 2 groups -- one
  nameChanges <- equivalentName(names(speciesLayers), sppEquiv, column = sppEquivCol)
  nameChangeNA <- is.na(nameChanges)
  names(speciesLayers)[!nameChangeNA] <- nameChanges[!nameChangeNA]

  nameChangesNonMerged <- equivalentName(names(speciesLayers)[nameChangeNA], sppEquiv, column = sppEquivCol)
  names(speciesLayers)[nameChangeNA] <- nameChangesNonMerged

  ## remove layers that have less data than thresh (i.e. spp absent in study area)
  ## count no. of pixels that have biomass
  layerData <- Cache(sapply, X = speciesLayers, function(x) sum(x[] > 0, na.rm = TRUE))

  ## remove layers that had < thresh pixels with biomass
  belowThresh <- layerData < thresh
  if (any(belowThresh)) {
    message(names(belowThresh)[belowThresh], " was removed because it had no data")
    speciesLayers[belowThresh] <- NULL
  }

  ## return stack and updated species names vector
  stack(speciesLayers)
}

#' Function to sum rasters of species layers
#'
#' @param speciesLayers stack of species layers rasters
#' @param layersToSum names/indices of layers to be summed - optional
#' @param filenameToSave file path to save output raster
#' @param newLayerName name of the output raster layer
#'
#' @export
#' @importFrom raster calc stack writeRaster
sumRastersBySpecies <- function(speciesLayers, layersToSum, filenameToSave, newLayerName) {
  out <- raster::calc(raster::stack(speciesLayers[layersToSum]), sum)
  names(out) <- newLayerName
  writeRaster(out, filename = filenameToSave, datatype = "INT2U", overwrite = TRUE)
  out # Work around for Cache
}

#' Overlay layers within raster stacks
#'
#' Overlays rasters of different data resolution by filling in gaps in the highest
#' resolution raster with data available in lowest resolution one.
#' If only high or low resolution data are available, it will use it without
#' attempting to overlay.
#'
#' @param highQualityStack      high quality list/stack of rasters (will be used preferencially)
#' @param lowQualityStack       low quality list/stack of rasters (will be used to fill NAs in highQualityStack)
#' @param outputFilenameSuffix  file suffix to save raster if there was overlaying. Defautls to "overlay"
#' @param destinationPath       directory for saved rasters
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom quickPlot layerNames
#' @importFrom raster ncell res stack
overlayStacks <- function(highQualityStack, lowQualityStack, outputFilenameSuffix = "overlay",
                          destinationPath) {
  ## check if HQ resolution > LQ resolutions
  hqLarger <- ncell(lowQualityStack) * prod(res(lowQualityStack)) <
    ncell(highQualityStack) * prod(res(highQualityStack))

  ## make table of species layers in HQ and LQ
  dt1 <- data.table(SPP = layerNames(highQualityStack), HQ = layerNames(highQualityStack))
  dt2 <- data.table(SPP = layerNames(lowQualityStack), LQ = layerNames(lowQualityStack))
  setkey(dt1, SPP); setkey(dt2, SPP)
  dtj <- merge(dt1, dt2, all = TRUE)
  dtj[, c("HQ", "LQ") := list(!is.na(HQ), !is.na(LQ))]

  ## check which layers have species info in HQ and LQ
  #dtj[, HQ := any(!is.na(highQualityStack[[SPP]][])), by = 1:nrow(dtj)] #nolint
  #dtj[, LQ := any(!is.na(lowQualityStack[[SPP]][])), by = 1:nrow(dtj)] #nolint

  stackRas <- list()
  for (x in seq(nrow(dtj))) {
    stackRas[[x]] <- dtj[x, .overlay(SPP, HQ, LQ, hqLarger = hqLarger,
                                     highQualityStack = highQualityStack,
                                     lowQualityStack = lowQualityStack,
                                     outputFilenameSuffix = outputFilenameSuffix,
                                     destinationPath = destinationPath)]
  }
  names(stackRas) <- dtj$SPP

  stack(stackRas)
}

#' Overlaying function
#'
#' Used internally in \code{overlayStacks}. Function to be applied to each row
#' of a \code{data.table} containing information of whether the species layer
#' exists in the HQ and LQ data.
#' Only overlays if data exists in both layers, otherwise returns the layer with data.
#'
#' @inheritParams overlayStacks
#' @param SPP \code{data.table} column of species layer name
#' @param HQ \code{data.table} column of whether SPP is present in HQ layers
#' @param LQ \code{data.table} column of whether SPP is present in LQ layers
#'
#' @importFrom gdalUtils gdalwarp
#' @importFrom raster compareRaster crs extent filename res projectExtent raster
#' @importFrom raster writeRaster xmax xmin ymax ymin
#' @keywords internal
.overlay <- function(SPP, HQ, LQ, hqLarger, highQualityStack, lowQualityStack, #nolint
                     outputFilenameSuffix = "overlay", destinationPath) {
  ## if HQ & LQ have data, pool
  if (HQ & LQ) {
    ## check equality of raster attributes and correct if necessary
    if (!all(
      isTRUE(all.equal(extent(lowQualityStack), extent(highQualityStack))),
      isTRUE(all.equal(crs(lowQualityStack), crs(highQualityStack))),
      isTRUE(all.equal(res(lowQualityStack), res(highQualityStack))))) {
      message("  ", SPP, " extents, or resolution, or projection did not match; ",
              "using gdalwarp to make them overlap")
      LQRastName <- basename(tempfile(fileext = ".tif"))
      if (!nzchar(filename(lowQualityStack[[SPP]]))) {
        LQCurName <- basename(tempfile(fileext = ".tif"))
        lowQualityStack[[SPP]][] <- as.integer(lowQualityStack[[SPP]][])
        lowQualityStack[[SPP]] <- writeRaster(lowQualityStack[[SPP]],
                                              filename = LQCurName,
                                              datatype = "INT2U")
      }

      LQRastInHQcrs <- projectExtent(lowQualityStack, crs = crs(highQualityStack))
      # project LQ raster into HQ dimensions
      gdalwarp(overwrite = TRUE,
               dstalpha = TRUE,
               s_srs = as.character(crs(lowQualityStack[[SPP]])),
               t_srs = as.character(crs(highQualityStack[[SPP]])),
               multi = TRUE, of = "GTiff",
               tr = res(highQualityStack),
               te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                      xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
               filename(lowQualityStack[[SPP]]), ot = "Byte",
               LQRastName)

      LQRast <- raster(LQRastName)
      LQRast[] <- LQRast[]
      unlink(LQRastName)

      try(unlink(LQCurName), silent = TRUE)

      if (hqLarger) {
        tmpHQName <- basename(tempfile(fileext = ".tif"))

        gdalwarp(overwrite = TRUE,
                 dstalpha = TRUE,
                 s_srs = as.character(crs(highQualityStack[[SPP]])),
                 t_srs = as.character(crs(highQualityStack[[SPP]])),
                 multi = TRUE, of = "GTiff",
                 tr = res(highQualityStack),
                 te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                        xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
                 filename(highQualityStack[[SPP]]), ot = "Byte", tmpHQName)
        HQRast <- raster(tmpHQName)
        HQRast[] <- HQRast[]
        HQRast[HQRast[] == 255] <- NA_integer_
        unlink(tmpHQName)
      } else {
        HQRast <- highQualityStack[[SPP]]
      }
    } else {
      LQRast <- lowQualityStack[[SPP]]
      HQRast <- highQualityStack[[SPP]]
    }

    message("  Writing new, overlaid ", SPP, " raster to disk.")
    if (!compareRaster(LQRast, HQRast))
      stop("Stacks not identical, something is wrong with overlayStacks function.")

    NAs <- is.na(HQRast[])

    ## complete missing HQ data with LQ data
    HQRast[NAs] <- LQRast[][NAs]
    HQRast <- writeRaster(HQRast, datatype = "INT1U",
                          filename = file.path(destinationPath,
                                               paste0(SPP, "_", outputFilenameSuffix, ".tif")),
                          overwrite = TRUE)
    names(HQRast) <- SPP
    return(HQRast)
  } else {
    ## if only HQ/LQ exist return one of them
    ## if none have data return one of the empty to keep all layers
    if (HQ) {
      HQRast <- highQualityStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    } else if (LQ) {
      LQRast <- lowQualityStack[[SPP]]
      names(LQRast) <- SPP
      return(LQRast)
    } else {
      HQRast <- highQualityStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    }
  }
}

#' Merge species pecent-cover rasters
#'
#' Used internally in \code{overlayStacks}.
#'
#' @param sppMerge TODO
#' @param speciesLayers TODO
#' @param sppEquiv TODO
#' @param column TODO
#' @param dPath destination path TODO
#' @param suffix TODO
#' @param ... Additonal arguments TODO
#'
#' @importFrom raster calc stack writeRaster
#' @importFrom reproducible .suffix prepInputs
#' @keywords internal
mergeSppRaster <- function(sppMerge, speciesLayers, sppEquiv, column, suffix, dPath, ...) {
    ## make sure species names and list names are in the right formats
    names(sppMerge) <- sppMerge
    sppMerges <- lapply(sppMerge, FUN = function(x) {
      equivalentName(x, sppEquiv,  column = "KNN", multi = TRUE)
    })
    #names(sppMerges) <- equivalentName(names(sppMerges), sppEquiv,  column = sppEquivCol)

    ## keep species present in the data
    sppMerges <- lapply(sppMerges, FUN = function(x) x[x %in% names(speciesLayers)])

    for (i in seq(length(sppMerges))) {
      sumSpecies <- sppMerges[[i]]
      newLayerName <- names(sppMerges)[i]

      fname <- .suffix(file.path(dPath, paste0("kNN", newLayerName, ".tif")), suffix)
      a <- calc(stack(speciesLayers[sumSpecies]), sum, na.rm = TRUE)
      names(a) <- newLayerName
      a <- writeRaster(a, filename = fname, overwrite = TRUE, ...)
      ## replace spp rasters by the summed one
      speciesLayers[sumSpecies] <- NULL
      speciesLayers[[newLayerName]] <- a
      message("  Merging ", paste(sumSpecies, collapse = ", "), "; becoming: ", newLayerName)
    }
  speciesLayers
}
