utils::globalVariables(c(
  "Degree_Of_", "Permafrost", "thermokarst", "dists"
))

## ------------------------------------------------------
## FUNCTIONS TO AID IN PREPARATION OF PERMAFROST LAYERS
## ------------------------------------------------------
## LandR may not be the best package for these functions, but they will live
## here for now.

#' Source and prepare permafrost, wetland and land-cover
#'   data and create a permafrost presence-absence raster
#'
#' Obtains permafrost peatland complex % data from Gibson et al. (2021),
#'    land-cover data for 2017 from Hermosilla et al. (2022) and wetland
#'    presence/absence data from Wulder et al. (2018), which are crossed
#'    to create a permafrost raster layer.
#'    See `?makeSuitForPerm` and `?assignPermafrost` for details on how this
#'    is achieved.
#'
#' @param cores integer. Number of threads to use for parallelisation (during
#'   permafrost raster creation). Defaults to no parallelisation (1 thread)
#' @param outPath character. A folder path where temporary permafrost rasters
#'   will be saved. Note that these *will not* be deleted at the end of this function
#' @cacheTags character. Common tags to use across `reproducible::Cache` calls. These will
#'   be appended to other tags distinguishing the nested `reproducible::Cache` calls.
#' @param destinationPath character. Passed to `reproducible::prepInputs` and the
#'   location where the final permafrost layer will be saved.
#' @param studyArea character. Passed to `reproducible::postProcess`.
#'
#' @return a `SpatRaster` of permafrost presences (1) and absences (0) with the
#'   same properties as the land-cover and wetland rasters, cropped and masked
#'   to `studyArea`.
#'
#' @references Gibson, C., Morse, P. D., Kelly, J. M., Turetsky, M. R., Baltzer, J. L., Gingras-Hill, T., & Kokelj, S. V. (2020). Thermokarst Mapping Collective: Protocol for organic permafrost terrain and preliminary inventory from the Taiga Plains test area, Northwest Territories (NWT Open Report 2020-010, p. 29). Northwest Territories Geological Survey.
#' @references Hermosilla, T., Bastyr, A., Coops, N. C., White, J. C., & Wulder, M. A. (2022). Mapping the presence and distribution of tree species in Canada’s forested ecosystems. Remote Sensing of Environment, 282, 113276. https://doi.org/10.1016/j.rse.2022.113276
#' @references Wulder, M., Li, Z., Campbell, E., White, J., Hobart, G., Hermosilla, T., & Coops, N. (2018). A National Assessment of Wetland Status and Trends for Canada’s Forested Ecosystems Using 33 Years of Earth Observation Satellite Data. Remote Sensing, 10(10), 1623. https://doi.org/10.3390/rs10101623
#'
#' @importFrom reproducible Cache prepInputs postProcess studyAreaName CacheDigest
#' @importFrom terra writeRaster writeVector project subst mask
#' @export
makePermafrostRas <- function(cores = 1L, outPath = getOption("spades.outputPath"),
                              cacheTags = NULL, destinationPath = getOption("reproducible.destinationPath"),
                              studyArea = NULL) {
  permafrostSA <- Cache(prepInputs,
                        targetFile = "TP_DP_StudyBoundary.shp",
                        alsoExtract = "similar",
                        archive = "DOI_2020-010.zip",
                        overwrite = TRUE,
                        fun = "terra::vect",
                        # useCache = "overwrite",
                        userTags = c(cacheTags, "permafrostSA"),
                        omitArgs = c("userTags"))
  permafrostSAName <- studyAreaName(permafrostSA)

  permafrostPoly <- Cache(getPermafrostDataGibson,
                          destinationPath = destinationPath,
                          studyArea = permafrostSA,
                          cacheTags = cacheTags,
                          userTags = c(cacheTags, "default", "permafrost"),
                          omitArgs = c("userTags"))

  rstLCC <- Cache(prepInputs,
                  targetFile = "CA_forest_VLCE2_2017.tif",
                  alsoExtract = "similar",
                  archive = "CA_forest_VLCE2_2017.zip",
                  studyArea = permafrostSA,
                  useSAcrs = FALSE,
                  overwrite = TRUE,
                  filename2 = .suffix("rstLCCPermafrost.tif", paste0("_", permafrostSAName)),
                  url = "https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_2017.zip",
                  userTags = c(cacheTags, "rstLCC", "permafrost", permafrostSAName),
                  omitArgs = c("userTags", "overwrite"))


  wetlands <- Cache(prepInputs,
                    targetFile = "CA_wetlands_post2000.tif",
                    archive = "CA_wetlands_post2000.zip",
                    studyArea = permafrostSA,
                    useSAcrs = FALSE,
                    overwrite = TRUE,
                    filename2 = .suffix("wetlandsPermafrost.tif", paste0("_", permafrostSAName)),
                    url = "https://opendata.nfis.org/downloads/forest_change/CA_wetlands_post2000.zip",
                    userTags = c(cacheTags, "wetlands", "permafrost", permafrostSAName),
                    omitArgs = c("userTags"))

  cTags <- c(cacheTags, "suitForPerm", permafrostSAName)
  suitForPerm <- Cache(makeSuitForPerm,
                       rstLCC = rstLCC,
                       wetlands = wetlands,
                       studyArea = permafrostSA,
                       cacheTags = cTags,
                       userTags = cTags,
                       omitArgs = c("userTags"))

  ## project permafrost polygons if need be
  if (!.compareCRS(permafrostPoly, suitForPerm)) {
    permafrostPoly <- project(permafrostPoly, crs(suitForPerm))
  }

  ## avoid large cached file by digesting earlier.
  preDigest <- CacheDigest(suitForPerm, permafrostPoly)
  preDigest <- preDigest$outputHash

  ## pass spatial objs as file names for compatibility with parallelisation
  tempFileRas <- paste0(tempfile(), ".tif")
  tempFileVect <- paste0(tempfile(), ".gpkg")
  writeRaster(suitForPerm, filename = tempFileRas)
  writeVector(permafrostPoly, filename = tempFileVect)

  # test <- list.files(outPath, "tempFiles"))
  # test <- sub("permafrost_polyID", "", sub("\\.tif", "", test))
  # Couldn't assign permafrost to enough pixels. id: 207217
  permafrostRas <- Cache(assignPermafrostWrapper,
                         gridPoly = tempFileVect,
                         ras = tempFileRas,
                         id = permafrostPoly$OBJECTID,
                         # id = setdiff(permafrostPoly$OBJECTID, test),
                         # id = 206782, ## a tiny polygon in SA edge
                         cores = cores,
                         # cores = 1L, ## test
                         saveDir = file.path(outPath, "tempFiles"),
                         outFilename = .suffix(file.path(destinationPath, "permafrost.tif"), paste0("_", permafrostSAName)),
                         dPath = destinationPath,
                         .cacheExtra = list(preDigest),
                         userTags = c(cacheTags, "permafrost", permafrostSAName),
                         omitArgs = c("userTags", "gridPoly", "ras", "cores", "saveDir",
                                      "outFilename", "dPath"))

  file.remove(tempFileRas, tempFileVect)

  ## set NAs to 0s inside permafrost study area
  permafrostRas <- subst(permafrostRas, NA, 0)
  permafrostSA <- project(permafrostSA, crs(permafrostRas))
  permafrostRas <- mask(permafrostRas, permafrostSA)

  ## now crop/mask to SA
  permafrostRas <- postProcess(permafrostRas,
                               studyArea = studyArea,
                               useSAcrs = FALSE,
                               filename2 = .suffix(file.path(destinationPath, "permafrost.tif"), paste0("_", .studyAreaName)),
                               overwrite = TRUE,
                               userTags = c(cacheTags, "permafrost", .studyAreaName),
                               omitArgs = c("userTags", "filename2"))
  ## make sure layer name is "permafrost"
  names(permafrostRas) <- "permafrost"
  permafrostRas
}


#' Make a mask of areas that suitable/unsuitable for permafrost
#'
#' Uses a land-cover map and a wetlands map to create a mask
#'  of areas suitable for permafrost presence. By default, it
#'  assumes land-cover classification follows Hermosilla et al. (2022)
#'  and that the wetlands map is a wetland presence/absence map.
#'  See details below.
#'
#' @param rstLCC a land-cover raster. Should already have been cropped to
#'   `studyArea`.
#' @param wetlands a wetlands raster (with `1L` and `NA`s only). Should already
#'   have been cropped to `studyArea`.
#' @param suitableCls land-cover classes (in `rstLCC`) suitable
#'  for permafrost presence
#' @param unsuitableCls land-cover classes (in `rstLCC`) unsuitable
#'  for permafrost presence
#' @param studyArea a polygon of the study area to mask the
#'  output raster to.
#' @param cacheTags passed to `reproducible::Cache`
#'
#' @details Suitable areas are classified as `1L` and unsuitable
#'  areas as `0L` in the output raster.
#'  Unsuitable land-cover classes are:
#'    * non-forested ecozone, which was not mapped for LC (0)
#'    * water (20)
#'    * snow/ice (31)
#'    * rock/rubble (32)
#'    * exposed/barren land (33)
#'    * wetland (80)
#'    * wetland-treed (81)
#'  Wetland areas in `wetlands` (coded as 1L, i.e. wetland presence)
#'  are also coded as unsuitable.
#'
#' @return a raster layer with `1`s and `0`s for pixels with suitable
#'  and unsuitable land cover for permafrost, respectively.
#'
#' @references Hermosilla, T., Bastyr, A., Coops, N. C., White, J. C., & Wulder, M. A. (2022). Mapping the presence and distribution of tree species in Canada’s forested ecosystems. Remote Sensing of Environment, 282, 113276. https://doi.org/10.1016/j.rse.2022.113276
#'
#' @importFrom terra classify subst mask
#' @importFrom raster raster
#' @importFrom reproducible Cache postProcess
#'
#' @export
makeSuitForPerm <- function(rstLCC, wetlands, suitableCls = c(40, 50, 100, 210, 220, 230),
                            unsuitableCls = c(0, 20, 31, 32, 33, 80, 81),
                            studyArea, cacheTags) {
  unsuit <- cbind(unsuitableCls, NA_integer_)
  suit <- cbind(suitableCls, 1L)
  m <- rbind(unsuit, suit)

  rstLCC <- Cache(classify,
                  x = rstLCC,
                  rcl = m,
                  userTags = c(cacheTags, "rstLCC", "suitablePermafrost"),
                  omitArgs = c("userTags"))
  rstLCC <- Cache(mask,
                  x = rstLCC,
                  mask = wetlands,
                  inverse = TRUE,
                  userTags = c(cacheTags, "rstLCC", "suitablePermafrost"),
                  omitArgs = c("userTags"))

  if (getOption("LandR.assertions", TRUE)) {
    test <- unique(rstLCC[!is.na(as.vector(values(wetlands)))])
    if (any(!is.na(test))) {
      stop("Something went wrong. All wetland areas should have 'NA' in suitable areas for permafrost")
    }
  }

  ## covert NAs to zeros then mask again
  rstLCC <- subst(rstLCC, NA, 0L)
  rstLCC <- postProcess(x = rstLCC,
                        studyArea = studyArea,
                        useSAcrs = FALSE,
                        userTags = c(cacheTags, "rstLCC", "suitablePermafrost"),
                        omitArgs = c("userTags"))

  return(rstLCC)
}



#' Create permafrost P/A layer
#'
#' Creates a permafrost presence/absence raster layer based
#'  on information about land-cover suitability for permafrost
#'  (a raster layer at high resolution) and % of permafrost
#'  present (a polygon layer that groups several cells of the
#'  land-cover layer, i.e. at a larger spatial scale). The function
#'  has been designed to be looped over polygons, and so only processes
#'  only polygon at a time.
#'
#' @param gridPoly a `SpatVector`, `PackedSpatVector` or `character`.
#'  The full polygon layer containing information about % permafrost
#'  (as a `SpatVector` or `PackedSpatVector`) or the file name to that
#'  layer.
#' @param ras a `SpatRaster`, `PackedSpatRaster` or `character`.
#'  The land-cover suitability for permafrost raster layer (as a
#'  `SpatRaster` or `PackedSpatRaster`) or the file name to that layer.
#'  Suitable cells are coded as `rasClass`.
#' @param saveOut logical. Should the processed rasters for the focal
#'  polygon be saved (as the file name returned) or directly returned?
#'  If parallelising, choose to save, as `SpatRasters` cannot be serialized.
#' @param saveDir character. The directiory to save output rasters.
#' @param id character or numeric. The polygon ID in `gridPoly[IDcol]`
#'  to process (the focal polygon).
#' @param IDcol character. Column name in `gridPoly` containing polygon IDs.
#' @param rasClass Class
#'
#' @return a file name or the permafrost presence/absence raster for the focal
#'  polygon.
#'
#' @importFrom terra rast vect unwrap writeRaster as.polygons as.points as.lines
#' @importFrom terra mask crop cellFromXY crds disagg distance expanse nearby fillHoles
#' @importFrom data.table as.data.table
#' @importFrom crayon cyan
#' @export
assignPermafrost <- function(gridPoly, ras, saveOut = TRUE, saveDir = NULL,
                             id = NULL, IDcol = "OBJECTID", rasClass = 1L) {
  message(cyan("Preparing layers..."))

  if (is(gridPoly, "PackedSpatVector")) {
    gridPoly <- unwrap(gridPoly)
  }
  if (is(ras, "PackedSpatRaster")) {
    ras <- unwrap(ras)
  }

  if (is(gridPoly, "character")) {
    gridPoly <- vect(gridPoly)
  }

  if (is(ras, "character")) {
    ras <- rast(ras)
  }

  if (is.null(id)) {
    id <- gridPoly[[IDcol]][1,]
  }
  idd <- which(gridPoly[[IDcol]] == id)
  landscape <- gridPoly[idd]
  sub_ras <- crop(ras, landscape, mask = TRUE, touches = FALSE)
  ## certain polygons may be in ras boundaries and so small
  ## that they don't overlap any cell centroids. So try again with touches = TRUE
  if (isFALSE(any(!is.na(as.vector(sub_ras[]))))) {
    sub_ras <- crop(ras, landscape, mask = TRUE)
  }

  ## the solution may not have worked and there may simply be no
  ## raster cells intersected by this polygon.
  skipThisPoly <- isFALSE(any(!is.na(as.vector(sub_ras[]))))

  ## make storage raster
  sub_rasOut <- sub_ras
  sub_rasOut[] <- NA_integer_

  if (skipThisPoly) {
    warning(paste0("None of 'ras' touch polygon id ", id,
                   ".\n  You may want to consider extending 'ras'."))
  } else {
    permpercent <- as.numeric(landscape[["Permafrost"]])

    ## and in pixels
    suitablePixNo <- sum(as.vector(sub_ras[]) == rasClass, na.rm = TRUE)
    permpercentPix <- (permpercent/100)*ncell(sub_ras)  ## don't round here, otherwise values <0.5 become 0.

    ## use max(..., 1) to guarantee that values lower than 1, get one pixel.
    if (permpercentPix > 0)
      permpercentPix <- max(permpercentPix, 1)

    permpercentPix <- round(permpercentPix)

    ## if there's less permafrost than suitable areas
    ## find the largest patch and a point that is distant from its edge
    ## then assign permafrost starting from this point
    ## until the percentage is reached

    message(cyan("Assigning permafrost..."))
    ## if there are more suitable areas than permafrost start
    ## "filling in" from focal pixels using distance to edges
    ## as the probability of a pixel being selected (more distant = higher prob)
    if (suitablePixNo > permpercentPix & permpercentPix > 0) {
      ## make polygons, so that we can identify separate patches in raster
      sub_poly <- as.polygons(sub_ras) |> disagg()   ## 0s in sub_ras are ignored
      names(sub_poly) <- "patchType"
      ## subset to patches of interest
      sub_poly <- sub_poly[sub_poly$patchType == rasClass, ]
      sub_poly$ID <- 1:nrow(sub_poly)

      ## now make sub_ras with poly IDs
      sub_ras2 <- rasterize(sub_poly, sub_ras, field = "ID")
      # terra::plot(sub_ras2, col = viridis::viridis(length(sub_poly)))

      ## compute distances to edges.
      ## first make an "inverse" raster with NAs
      sub_rasDist <- sub_ras2
      sub_rasDist[] <- 1L
      sub_rasDist <- mask(sub_rasDist, sub_ras2, inverse = TRUE)   ## use sub_ras2, because 0s were NAed
      sub_rasDist <- distance(sub_rasDist)

      ## make "probabilities" by scaling 0-1
      spreadProb <- sub_rasDist/max(sub_rasDist[], na.rm = TRUE)
      # terra::plot(spreadProb, col = viridis::inferno(100))

      ## we may need to try several times until we get the number of pixels
      ## at each attempt increase spreadProb
      sub_rasOut <- assignPresences(assignProb = spreadProb,
                                    landscape = sub_ras2,
                                    pixToConvert = permpercentPix,
                                    probWeight = 0.5, numStartsDenom = 10)
      # terra::plot(sub_rasOut, col = viridis::inferno(100))
    }

    ## if there's more permafrost than suitable areas
    ## assign all suitable areas as permafrost, then increase
    ## with a buffer (starting with largest patch)
    if (suitablePixNo <= permpercentPix & permpercentPix > 0) {
      pixToConvert <- permpercentPix

      if (suitablePixNo > 0) {
        ## convert all available cells
        cellIDs <- which(as.vector(sub_ras[]) == rasClass)
        sub_rasOut[cellIDs] <- 1L

        ## if there are not enough points within the patches,
        ## try to fill neighbouring pixels (leave holes alone s we want to be
        ## able to start from a swiss cheese pattern)
        pixToConvert2 <- pixToConvert - length(cellIDs)

        while (pixToConvert2 > 0) {
          sub_poly <- as.polygons(sub_rasOut) |> disagg()

          sub_poly_filled <- fillHoles(sub_poly)
          sub_poly_filled <- rasterize(sub_poly_filled, sub_ras)

          ## buffer around patch
          sub_poly_filledBuffer <- buffer(sub_poly_filled, width = unique(res(sub_poly_filled)),
                                          background = NA)
          cellIDs <- which(as.vector(sub_poly_filledBuffer[]) == 1)

          ## remove cellIDs that have been converted already
          vals <- sub_rasOut[cellIDs][]
          cellIDs <- cellIDs[is.na(vals)]

          ## check that we don't have too many cells. if we do, keep
          ## pixels around larger patches
          if (length(cellIDs) > pixToConvert2) {
            sub_rasOut2 <- sub_rasOut
            sub_rasOut2[cellIDs] <- 1L

            sub_poly2 <- as.polygons(sub_rasOut2) |> disagg()
            sub_poly2$area <- expanse(sub_poly2)
            sub_rasOutArea <- rasterize(sub_poly2, sub_rasOut2, field = "area")

            DT <- data.table(cells = 1:ncell(sub_rasOutArea), area = as.vector(sub_rasOutArea[]))
            DT <- DT[complete.cases(DT)][cells %in% cellIDs]
            DT <- DT[order(area, decreasing = TRUE)]
            cellIDs <- DT[1:pixToConvert2, cells]
          }

          ## we may have exhausted areas to fill outside the holes
          ## so we fill holes (preferentially the smaller ones)
          if (length(cellIDs) < 1) {
            sub_poly <- as.polygons(sub_rasOut) |> disagg()
            sub_poly_holes <- fillHoles(sub_poly, inverse = TRUE) |> disagg()
            sub_poly_holes$area <- expanse(sub_poly_holes)

            sub_rasHolesArea <- rasterize(sub_poly_holes, sub_rasOut, field = "area")
            ## there may be patches inside holes (like islands) that need
            ## to be masked out, as they are ignored by fillHoles
            sub_rasHolesArea <- mask(sub_rasHolesArea, sub_rasOut, inverse = TRUE)

            DT <- data.table(cells = 1:ncell(sub_rasHolesArea), area = as.vector(sub_rasHolesArea[]))
            DT <- DT[complete.cases(DT)]
            setorder(DT, area, cells)

            cellIDs <- DT[1:pixToConvert2, cells]
          }

          sub_rasOut[cellIDs] <- 1L

          pixToConvert2 <- pixToConvert2 - length(cellIDs)
        }
      } else {
        ## there may be no available pixels, in which case permafrost can be assigned
        ## starting in a random polygon within unsuitable areas (hence equal prob below)
        spreadProb <- subst(sub_ras, c(NA, 0L), c(0L, 1L))

        sub_rasOut <- assignPresences(assignProb = spreadProb,
                                      landscape = sub_ras,
                                      pixToConvert = pixToConvert,
                                      probWeight = 1, numStartsDenom = 10)

      }
    }

    if (sum(!is.na(sub_rasOut[])) < permpercentPix) {
      warning(paste("Couldn't assign permafrost to enough pixels. id:", id))
    }
    message(cyan("Done!"))
  }

  if (saveOut) {
    tmpFile <- paste0("permafrost_polyID", id, ".tif")
    if (!is.null(saveDir)) {
      if (!dir.exists(saveDir)) dir.create(saveDir, showWarnings = FALSE)
      tmpFile <- file.path(saveDir, basename(tmpFile))
    }
    writeRaster(sub_rasOut, tmpFile, overwrite = TRUE)
    return(tmpFile)
  } else {
    return(sub_rasOut)
  }
}

#' Assign presences to patches based on a raster or presence
#'   probabilities
#'
#' @param assignProb a SpatRaster of presence probabilities
#' @param landscape a SpatRaster of the entire lanscape with NAs outside
#'   patches
#' @param pixToConvert numeric. Number of pixels to convert to presence
#'   across `landscape`. If `NULL`, `round(sum(!is.na(landscape[]))/2)`
#'   pixels will be converted to presences.
#' @param probWeight numeric. If `pixToConvert` cannot be reached
#'  `assignProb` will be weighted using `assignProb^probWeight`
#'  and the algorothm will try again.
#' @param numStartsDenom integer. Used to calculate the number of starting pixels
#'  to assign presences (at each try; as `pixToConvert/numStartsDenom`)
#'
#' @details This function attempts to iteratively assign
#'  presences starting in `pixToConvert/numStartsDenom` pixels
#'  sampled from areas with high probabilities in `assignProb` (weighted if
#'  `probWeight != 1`).
#'  For each starting pixel, it assigns a maximum of `numStartsDenom * probWeight`
#'  pixels. If not enough pixels are assigned presences (i.e. the number of
#'  presences does not reach `pixToConvert`), the function tries again
#'  after increasing `assignProb` by a factor of `1.5` (with values > 1 capped at 1).
#'  If too many pixels are assigned presences, pixels are sampled according to
#'  `assignProb ^ 10` (where `assignProb` corresponds to the original supplied values.)
#'
#' @return a SpatRaster
#' @importFrom SpaDES.tools spread
#' @importFrom raster raster
#' @importFrom terra rast
#' @importFrom data.table data.table
#'
#' @export
assignPresences <- function(assignProb, landscape, pixToConvert = NULL, probWeight = 1,
                            numStartsDenom = 10) {
  if (is.null(pixToConvert)) {
    pixToConvert <- round(sum(!is.na(landscape[]))/2)
  }

  if (probWeight < 0.5 || probWeight > 7)
    stop("probWeight must be between 0.5 and 7")

  ## save original probabilites for later
  assignProbOrig <- as.vector(assignProb[])

  ## if the mean is too high, then bring it down to 0.35 to avoid creating
  ## square patches
  meanSP <- mean(assignProbOrig, na.rm = TRUE)
  if (meanSP > 0.35)
    assignProb <- assignProb / meanSP * 0.35

  convertedPix <- 0

  while (convertedPix < pixToConvert) {
    ## exponentiate probabilities to provide more weight to pixels further from edges
    assignProbEx <- assignProb^probWeight
    # terra::plot(assignProbEx, col = viridis::inferno(100))

    ## try to spread in from many focal pixels
    numStarts <- ceiling(pixToConvert/numStartsDenom)
    startPoints <- sample(1:ncell(assignProbEx), size = numStarts, prob = assignProbEx[])

    # temp <- assignProbEx
    # temp[] <- NA_integer_
    # temp[startPoints] <- 1L
    # terra::plot(assignProbEx, col = viridis::inferno(100))
    # terra::plot(temp, add = TRUE, col = "blue")

    outRas <- SpaDES.tools::spread(landscape = raster(landscape),
                                   loci = startPoints,
                                   assignProb = raster(assignProbEx),
                                   maxSize = numStartsDenom * probWeight)
    outRas <- rast(outRas) |>
      mask(mask = landscape)
    convertedPix <- sum(!is.na(outRas[]))

    ## increase spread probabilities in case we need to try again
    assignProb <- min(assignProb * 1.5, 1)
  }

  ## if we spread too much remove pixels that are closest to edges
  if (convertedPix > pixToConvert) {
    ## "convert" to distances
    DT <- data.table(cells = 1:ncell(outRas), dists = assignProbOrig)
    DT <- DT[dists > 0 & !is.na(dists)][order(dists, decreasing = TRUE)]

    # pixToRm <- DT[(pixToConvert + 1):nrow(DT), cells] ## creates concave patches
    pixToKeep <- sample(DT$cells, pixToConvert, prob = DT$dists^10)  ## creates a little more noise

    outRas[] <- NA_integer_
    outRas[pixToKeep] <- 1L
  }
  convertedPix <- sum(!is.na(outRas[]))
  message(convertedPix, " were converted")
  return(outRas)
}
