utils::globalVariables(c(
  "Degree_Of_", "Permafrost", "thermokarst", "dists", "cells", "area"
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
#' @param cacheTags character. Common tags to use across `reproducible::Cache` calls. These will
#'   be appended to other tags distinguishing the nested `reproducible::Cache` calls.
#' @param destinationPath character. Passed to `reproducible::prepInputs` and the
#'   location where the final permafrost layer will be saved.
#' @param studyArea character. Passed to `reproducible::postProcess`.
#' @param .studyAreaName character. A string to append to `filename2` in `postProcess`.
#' @param doMissing logical. If TRUE, the permafrost layer will only be computed in
#'   polygon IDs missing from raster file names in `outPath`.
#' @param ... further arguments passed to `assignPermafrost`
#'
#'
#' @return a `SpatRaster` of permafrost presences (1) and absences (0) with the same
#'   properties as the land-cover and wetland rasters, cropped and masked to `studyArea`.
#'
#' @references Gibson, C., Morse, P. D., Kelly, J. M., Turetsky, M. R., Baltzer, J. L.,
#'             Gingras-Hill, T., & Kokelj, S. V. (2020). Thermokarst Mapping Collective:
#'             Protocol for organic permafrost terrain and preliminary inventory from the Taiga
#'             Plains test area, Northwest Territories (NWT Open Report 2020-010, p. 29).
#'             Northwest Territories Geological Survey.
#' @references Hermosilla, T., Bastyr, A., Coops, N. C., White, J. C., & Wulder, M. A. (2022).
#'             Mapping the presence and distribution of tree species in Canada's forested ecosystems.
#'             Remote Sensing of Environment, 282, 113276. https://doi.org/10.1016/j.rse.2022.113276
#' @references Wulder, M., Li, Z., Campbell, E., White, J., Hobart, G., Hermosilla, T., &
#'             Coops, N. (2018). A National Assessment of Wetland Status and Trends for Canada's
#'             Forested Ecosystems Using 33 Years of Earth Observation Satellite Data.
#'             Remote Sensing, 10(10), 1623. https://doi.org/10.3390/rs10101623
#'
#' @importFrom reproducible Cache prepInputs postProcess studyAreaName CacheDigest
#' @importFrom terra writeRaster writeVector project subst mask
#' @export
makePermafrostRas <- function(cores = 1L, outPath = getOption("spades.outputPath"),
                              cacheTags = NULL, destinationPath = getOption("reproducible.destinationPath"),
                              studyArea = NULL, .studyAreaName = NULL, doMissing = FALSE, ...) {
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

  toDoIDs <- permafrostPoly$OBJECTID
  tempSaveDir <- file.path(outPath, "tempFiles")
  if (doMissing) {
    doneIDs <- list.files(tempSaveDir)
    doneIDs <- sub("permafrost_polyID", "", sub("\\.tif", "", doneIDs))
    toDoIDs <- setdiff(toDoIDs, doneIDs)
  }

  argsList <- list(
    id = toDoIDs,
    # id = c(99, 401, 27557, 200841, 206782), ## several different cases
    cores = cores,
    # cores = 1L, ## test
    dPath = destinationPath,
    ## further arguments passed to assignPermafrost:
    gridPoly = tempFileVect,
    ras = tempFileRas,
    saveDir = tempSaveDir)

  dots <- list(...)
  notRecognised <- dots[!names(dots) %in% names(formals("assignPermafrost"))]
  if (length(notRecognised)) {
    stop("The following arguments do not match any 'assignPermafrost' arguments:",
         "\n", paste(names(notRecognised), collapse = ", "))
  }

  dots <- dots[names(dots) %in% names(formals("assignPermafrost"))]
  argsList <- append(argsList, dots)

  permafrostFiles <- do.call(.assignPermafrostWrapper, argsList) |>
    Cache(.cacheExtra = list(preDigest, toDoIDs),
          userTags = c(cacheTags, "permafrost", permafrostSAName),
          omitArgs = c("userTags", "gridPoly", "ras", "cores", "saveDir",
                       "dPath"))

  file.remove(tempFileRas, tempFileVect)

  ## read rasters and merge
  if (doMissing) {
    ## re-make list with all files
    permafrostRasFiles <- list.files(tempSaveDir, full.names = TRUE)
  }
  permafrostRasLs <- lapply(permafrostRasFiles, rast)
  permafrostRas <- sprc(permafrostRasLs)
  permafrostRas <- terra::merge(permafrostRas)

  message("Permafrost layer done!")

  outFilename <- .suffix(file.path(destinationPath, "permafrost.tif"), paste0("_", permafrostSAName))
  writeRaster(permafrostRas, filename = outFilename, overwrite = TRUE)

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
  return(permafrostRas)
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
#' @references Hermosilla, T., Bastyr, A., Coops, N. C., White, J. C., & Wulder, M. A. (2022).
#'             Mapping the presence and distribution of tree species in Canada's forested ecosystems.
#'             Remote Sensing of Environment, 282, 113276. https://doi.org/10.1016/j.rse.2022.113276
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

#' Wrapper for `assignPermafrost`
#'
#' @param cores how many threads to use for parallelisation. If
#'  `1`, no parallelisation occurs
#' @param id passed to `assignPermafrost`
#' @param dPath directory where output permafrost raster file will be written.
#' @param ... further arguments passed to `assignPermafrost`, like `gridPoly`
#'  `ras`, etc.
#'
#' @importFrom terra writeRaster sprc
.assignPermafrostWrapper <- function(id, cores, dPath, ...) {
  if (cores > 1) {
    message("Creating permafrost layer with parallelisation")
    if (.Platform$OS.type == "unix") {
      cl <- parallel::makeForkCluster(cores)
    } else {
      if (!requireNamespace("parallelly", quietly = TRUE)) {
        stop("install 'parallelly'")
      }
      cl <- parallelly::makeClusterPSOCK(cores, rscript_libs = .libPaths())
    }

    if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      stop("install 'RhpcBLASctl'")
    }

    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library("MASS")

      ## limit spawning of additional threads from workers
      data.table::setDTthreads(1)
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    })

    future::plan(future::cluster, workers = cl)
    applyFUN <- future.apply::future_Map
  } else {
    message("Creating permafrost layer sequentially")
    applyFUN <- Map
  }

  dots <- list(...)
  notRecognised <- dots[!names(dots) %in% names(formals("assignPermafrost"))]
  if (length(notRecognised)) {
    stop("The following arguments do not match any 'assignPermafrost' arguments:",
         "\n", paste(names(notRecognised), collapse = ", "))
  }

  assignPermafrostArgs <- dots[names(dots) %in% names(formals("assignPermafrost"))]

  ## use do.call so that future.seed can be added if running in parallel
  applyFUNArgs <- list(f = assignPermafrost,
                       id = id,
                       MoreArgs = assignPermafrostArgs)
  if (cores > 1) {
    applyFUNArgs <- append(applyFUNArgs, list(future.seed = TRUE))
  }
  permafrostRasFiles <- do.call(applyFUN, applyFUNArgs)
  return(permafrostRasFiles)
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
#' @param saveDir character. The directory to save output rasters.
#' @param id character or numeric. The polygon ID in `gridPoly[IDcol]`
#'  to process (the focal polygon).
#' @param IDcol character. Column name in `gridPoly` containing polygon IDs.
#' @param permafrostcol character. Column name in `gridPoly` containing polygon amount
#'  of permafrost (in %).
#' @param thermokarstcol character. Optional. Column name in `gridPoly` containing
#'  permafrost condition (thermokarst or degree of vulnerability). The function
#'  expects levels "low", "medium" or "high". See details below.
#' @param useMidpoint logical. If TRUE the midpoint of thermokarst percentage ranges
#'  is used to recalculate the amount of pixels to be converted to permafrost. Only
#'  applicable if `thermokarstcol` is not NULL.
#' @param rasClass Class
#'
#' @return a file name or the permafrost presence/absence raster for the focal
#'  polygon.
#'
#' @details If `thermokarstcol` is not NULL, the function expects one of "low",
#'  "medium" or "high" in `gridPoly[[thermokarstcol]]`, each corresponding to the
#'  following thermokarst percentages:
#'  * low: 0-33%
#'  * medium: 34-66%
#'  * high: 67-100%
#'  If `useMidpoint == TRUE` the midpoint value of each range is used as the percentage
#'  of pixels that have been thawed already and are not to be converted to permafrost
#'  (`thawedPercent` below). Otherwise, `thawedPercent` is randomly drawn from a uniform
#'  distribution bounded by the range minimum and maximum. The final number of pixels to convert
#'  to permafrost (`pixToConvert` below) is then calculated as:
#'  `pixToConvert = no. of suitable pixels \* (gridPoly[[Permafrost]] \* (100-thawedPercent))`.
#'
#' @importFrom crayon cyan
#' @importFrom data.table as.data.table setorder
#' @importFrom SpaDES.tools neutralLandscapeMap
#' @importFrom terra as.polygons buffer crop disagg distance expanse fillHoles
#' @importFrom terra mask minmax rast vect unwrap writeRaster
#' @export
assignPermafrost <- function(gridPoly, ras, saveOut = TRUE, saveDir = NULL,
                             id = NULL, permafrostcol = "Permafrost", thermokarstcol = NULL,
                             IDcol = "OBJECTID", rasClass = 1L, useMidpoint = TRUE) {
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
    ## extract % PPC and thermokarst level from polygon
    permpercent <- as.numeric(landscape[[permafrostcol]])

    if (!is.null(thermokarstcol) & permpercent > 0) {
      thermLevel <- tolower(landscape[[thermokarstcol]])
      ## there is at least one polygon with permafrost > 0 and NA thermokarst.
      ## We assume 0 thermokarst, which belongs to low thermLevel
      if (is.na(thermLevel)) {
        thermPercent <- 0
        thermLevel <- "low"
      } else {
        thermPercent <- .thermPercent(thermLevel, useMidpoint)
      }

      ## remove "thermokarsted" % amount
      ## because % thermokarst is relative to % PPC (gridPoly[[permafrostcol]])
      ## we cannot simply subtract
      permpercent <- permpercent * ((100 - thermPercent) / 100)
    }

    ## How many suitable pixels are there?
    suitablePixNo <- sum(as.vector(sub_ras[]) == rasClass, na.rm = TRUE)

    ## how many pixels do we need to convert (relative to all non-NA pixels -- note that ncell would count NAs)
    allPix <- sum(!is.na(sub_ras[]))
    permpercentPix <- (permpercent/100)*allPix  ## don't round here, otherwise values <0.5 become 0.

    ## use max(..., 1) to guarantee that values lower than 1, get one pixel.
    if (permpercentPix > 0)
      permpercentPix <- max(permpercentPix, 1)

    permpercentPix <- round(permpercentPix)

    ## if there's less permafrost than suitable areas
    ## find the largest patch and a point that is distant from its edge
    ## then assign permafrost starting from this point
    ## until the percentage is reached

    message(cyan("Assigning permafrost..."))
    ## if there are more suitable areas than permafrost:
    ## A) start "filling in" from focal pixels using distance to edges
    ## as the probability of a pixel being selected (more distant = higher prob)
    ## or B) fill largest patches first
    ## (Ceres tried option A alone, but assignPresences got stuck in fragmented landscape very often
    ## and it was hard to find a suitable weight/startPixel combination across varying degrees
    ## of fragmentation)
    if (suitablePixNo > permpercentPix & permpercentPix > 0) {
      pixToConvert <- permpercentPix

      ## make polygons, so that we can identify separate patches in raster
      sub_poly <- as.polygons(sub_ras) |> disagg()   ## 0s in sub_ras are ignored
      names(sub_poly) <- "patchType"
      ## subset to patches of interest
      sub_poly <- sub_poly[sub_poly$patchType == rasClass, ]
      sub_poly$ID <- 1:nrow(sub_poly)

      ## now make sub_ras with poly IDs
      sub_ras2 <- rasterize(sub_poly, sub_ras, field = "ID")
      # terra::plot(sub_ras2, col = viridis::viridis(length(sub_poly)))

      patchSizes <- as.data.table(table(sub_ras2[]))
      patchSizes <- patchSizes[order(N, decreasing = TRUE)]
      patchFillFirst <- which(cumsum(patchSizes$N) <= pixToConvert)
      patchFillFirst <- as.integer(as.character(patchSizes[patchFillFirst, V1]))

      sub_rasOut[as.vector(sub_ras2[]) %in% patchFillFirst] <- 1L
      # terra::plot(sub_ras)
      # terra::plot(sub_rasOut, col = "lightblue", add = TRUE)

      pixToConvert <- pixToConvert - sum(as.vector(sub_rasOut[]) == 1, na.rm = TRUE)

      ## the code bellow was often getting stuck in fragmented landscapes
      ## so we use it to fill in the missing pixels only.
      if (pixToConvert > 0) {
        ## compute distances to edges of patches that need to be filled
        sub_ras2 <- subst(sub_ras, c(0), NA_integer_)
        sub_ras2 <- mask(sub_ras2, sub_rasOut, inverse = TRUE)

        ## first make an "inverse" raster with NAs where patches are
        sub_rasDist <- sub_ras
        sub_rasDist[] <- 1L
        sub_rasDist <- mask(sub_rasDist, sub_ras2, inverse = TRUE)   ## use sub_ras2, because 0s were NAed
        sub_rasDist <- distance(sub_rasDist)   ## compute distance from pathes (NAs) to edges

        ## make "probabilities" by scaling 0-1
        spreadProb <- sub_rasDist / minmax(sub_rasDist)["max",]
        spreadProb <- mask(spreadProb, sub_ras2)   ## only NAs here are respected by spread
        # terra::plot(spreadProb, col = viridis::inferno(100))

        ## weight probabilities by patch size
        patchSizes[, V1 := as.integer(V1)]
        sub_rasPS <- rasterize(sub_poly, sub_ras, field = "ID")
        sub_rasPS <- classify(sub_rasPS, patchSizes)
        sub_rasPS <- sub_rasPS / minmax(sub_rasPS)["max",]   ## normalize
        spreadProb <- spreadProb * sub_rasPS

        suitablePixNo2 <- sum(sub_ras2[] == 1, na.rm = TRUE)
        availRatio <- suitablePixNo2/allPix


        ## Thermokarst levels (also) affect degree of fragmentation
        ## April 14th 2023: Ceres could not find a robust way of varying weights with
        ## thermokarst while ensuring feasible computation times and sensible patterns
        ## (small weights create very random patches when availRatio is high, large
        ## weights massively increase computation times when availRatio is <0.30 )
        ## April 18th 2023: thermokarst effects on weights were reintroduced and seem to
        ## be working as well as possible, when the midpoint is used in .thermPercent.
        ## Computation times are long (>10h for 12000 polygons)
        # if (FALSE) {
        if (!is.null(thermokarstcol)) {
          weight <- .thermWeight(thermLevel)
        } else {
          weight <- 3 ## a good average level of aggregation
        }

        ## in landscapes with very high availRatio, we can end up with very
        ## disaggregated/fragmented patterns if the weight < 7 and the number of
        ## starting pixels is high.
        ## So we bump the weights and keep a low no. starting pixels to increase
        ## aggregation. When availRatio is very low, increase number of starting pixels
        if (availRatio > 0.5) {
          weight <- weight*(1+availRatio*2)
        }
        # } else {
        #   weight <- 1  ## use probs as they are
        # }

        ## the number of starting pixels varies with availRatio
        ## with a negative exponential decay. If weight is high, we'll
        ## have a high degree of clumping which can lead to it being harder to
        ## generate presences, so we multiply the no. pixels by weight.
        # plot(((exp(seq(0, 1, length.out = 100)))^-4)*200 ~ seq(0,1,  length.out = 100), xlab = "availRatio", ylab = "noStartPix")
        noStartPix <- round((exp(availRatio)^-5)*200*weight)
        noStartPix <- pmax(noStartPix, 7)   ## don't go lower than 7
        noStartPix <- pmin(suitablePixNo2, noStartPix)   ## can't have more than the number of suitable pixels.

        sub_rasOut2 <- assignPresences(assignProb = spreadProb,
                                       landscape = sub_ras2,
                                       pixToConvert = pixToConvert,
                                       probWeight = weight,
                                       numStartsDenom = pixToConvert/noStartPix)
        # terra::plot(sub_ras)
        # terra::plot(sub_rasOut, col = "lightblue", add = TRUE)
        # terra::plot(sub_rasOut2, col = "darkblue", add = TRUE)
        sub_rasOut <- cover(sub_rasOut, sub_rasOut2)
      }
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
        ## try to fill neighbouring pixels, closer to edges
        pixToConvert2 <- pixToConvert - length(cellIDs)

        if (pixToConvert2 > 0) {
          sub_ras_filled <- sub_rasOut

          ## no longer necessary with assignPresences, bcs it creates a fragmented
          ## pattern by itself
          if (FALSE) {
            ## leave holes alone at first as we want to start from a
            ## swiss-cheese pattern
            # while (pixToConvert2 > 0) {   ## only needed when using buffers
            sub_poly <- as.polygons(sub_rasOut) |> disagg()   ## when using buffers
            ## don't fill holes first to keep a swiss-cheese pattern
            sub_poly_filled <- fillHoles(sub_poly)
            sub_ras_filled <- rasterize(sub_poly_filled, sub_ras)
          }

          ## calculate distance to patch from NAs
          sub_ras_filledDist <- distance(sub_ras_filled)
          ## to reverse to higher spreadProb closest to patches
          ## first rescale to 0-1 by dividing by max, then subtract 1
          ## and take abs (all -1s are converted to 0s) This is better than taking the inverse and recalling
          spreadProb <- (sub_ras_filledDist/minmax(sub_ras_filledDist)["max",])
          spreadProb <- spreadProb - 1
          spreadProb <- subst(spreadProb, -1, 0)
          spreadProb <- abs(spreadProb)
          spreadProb <- mask(spreadProb, sub_ras)   ## only NAs here are respected by spread

          weight <- 5    ## should ensure a good degree of clumping

          suitablePixNo2 <- sum(spreadProb[] > 0, na.rm = TRUE)
          availRatio <- suitablePixNo2/allPix
          noStartPix <- round((exp(availRatio)^-5)*200*weight)
          noStartPix <- pmax(noStartPix, 7)   ## don't go lower than 7
          noStartPix <- pmin(suitablePixNo2, noStartPix)   ## can't have more than the number of suitable pixels.

          sub_ras_filledBuffer <- assignPresences(assignProb = spreadProb,
                                                  landscape = sub_ras,
                                                  pixToConvert = pixToConvert2,
                                                  probWeight = weight,
                                                  numStartsDenom = pixToConvert2/10)
          # terra::plot(sub_ras)
          # terra::plot(sub_rasOut, col = "lightblue", add = TRUE)
          # terra::plot(sub_ras_filledBuffer, col = "darkblue", add = TRUE)

          ## buffer around patch -- this is creating square patterns when starting from a single pixel
          # sub_ras_filledBuffer <- buffer(sub_ras_filled, width = unique(res(sub_poly_filled)),
          # background = NA)
          cellIDs <- which(as.vector(sub_ras_filledBuffer[]) == 1)

          ## remove cellIDs that have been converted already -- shouldn't be necessary with assignPresences
          ## but kept here in case we switch to buffer again.
          vals <- sub_rasOut[cellIDs][]
          cellIDs <- cellIDs[is.na(vals)]

          ## this is no longer necessary with assignPresences
          if (FALSE) {
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

            ## This will no longer happen with assignPresences
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
          }

          sub_rasOut[cellIDs] <- 1L
          # pixToConvert2 <- pixToConvert2 - length(cellIDs)  ## only needed for the while, when using buffers

          # terra::plot(sub_ras)
          # terra::plot(sub_rasOut, col = "lightblue", add = TRUE)
        }
      } else {
        ## there may be no available pixels, in which case permafrost can be assigned
        ## starting in random pixels within unsuitable areas, using a neutral landscape
        ## of probabilities (to ensure some degree of agglomeration)
        spreadProb <- neutralLandscapeMap(sub_ras, type = "nlm_gaussianfield",
                                          autocorr_range = 100, mag_var = 50)
        weight <- 100   ## weights have to be "stupidly" high to ensure clumping

        sub_rasOut <- assignPresences(assignProb = spreadProb,
                                      landscape = sub_ras,
                                      pixToConvert = pixToConvert,
                                      probWeight = weight,
                                      numStartsDenom = pixToConvert/10)
        # terra::plot(sub_ras)
        # terra::plot(sub_rasOut, col = "lightblue", add = TRUE)

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

#' Assign presences to patches based on a raster or presence probabilities
#'
#' @param assignProb a `SpatRaster` of presence probabilities (`NA` allowed)
#' @param landscape a `SpatRaster` of the entire landscape with `NA` outside patches
#' @param pixToConvert numeric. Number of pixels to convert to presence across `landscape`.
#'                     If `NULL`, it will convert `round(sum(!is.na(landscape[]))/2)` pixels.
#' @param probWeight numeric. A weight for `assignProb` (`assignProb^probWeight`)
#'  which affects the degree of clumping (higher weights result in clumpier patterns)
#' @param numStartsDenom integer. Used to calculate the number of starting pixels
#'  (as `pixToConvert/numStartsDenom`) at each try of assigning presences.
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
#' @examples
#' library(terra)
#' library(raster)
#' library(SpaDES.tools)
#' set.seed(123)
#' ras <- raster(extent(0, 100, 0, 100), res = 1, vals = 0)
#' ras[] <- 1L
#' ras_nas <- randomPolygons(raster(extent(0, 100, 0, 100), res = 1, vals = 0), numTypes = 200)
#' ras_nas[ras_nas[] %in% sample(ras_nas[], 120)] <- NA
#' ras <- mask(ras, ras_nas)
#' ras <- rast(ras)

#' ## NA patches to calculate distance
#' ras_nas <- ras
#' ras_nas[] <- 1L
#' ras_nas <- mask(ras_nas, ras, inverse = TRUE)
#' plot(ras, col = "lightgreen")
#' plot(ras_nas, add = TRUE, col = "coral")
#'
#' ## distance from NAs
#' ras_distance <- distance(ras_nas)
#' plot(ras_distance)
#'
#' spreadProbRas <- ras_distance/minmax(ras_distance)["max",]
#' # spreadProbRas <- mask(spreadProbRas, ras)
#' # spreadProbRas[is.na(spreadProbRas[])] <- 0
#' spreadProbRas[as.vector(spreadProbRas[]) == 0] <- NA
#' plot(spreadProbRas)
#'
#' ## varying pixels and weights
#' weights <- seq(0.5, 7, length.out = 5)
#' ## we are in a case where we want fewer pixels than available
#' pixToConvert <- round(sum(!is.na(ras[]))/4 * seq(0.5, 2, length.out = length(weights)))
#' tests <- Map(probWeight = weights,   ## not much difference in clumping levels
#'              pixToConvert = pixToConvert,
#'              numStartsDenom = pixToConvert/10,
#'              f = assignPresences,
#'              MoreArgs = list(assignProb = spreadProbRas,
#'                              landscape = ras))
#'
#' names(tests) <- paste(paste("Prob. weight", weights),
#'                       paste("no pix", c(pixToConvert)),
#'                       sep = ";")
#'
#' plot(rast(tests), col = "lightblue")
#'
#' ## Function to iterate over several landscape configurations,
#' ## weights and no. pixels desired.
#' testFun <- function(NAgroups = 5, varyWeights = TRUE, varyPixels = TRUE) {
#'   set.seed(123)
#'   ras <- raster(extent(0, 100, 0, 100), res = 1, vals = 0)
#'   ras[] <- 1L
#'   ras_nas <- randomPolygons(raster(extent(0, 100, 0, 100), res = 1, vals = 0), numTypes = 200)
#'   ras_nas[ras_nas[] %in% sample(unique(ras_nas[]), NAgroups)] <- NA
#'   ras <- mask(ras, ras_nas)
#'   ras <- rast(ras)
#'   ## NA patches to calculate distance
#'   ras_nas <- ras
#'   ras_nas[] <- 1L
#'   ras_nas <- mask(ras_nas, ras, inverse = TRUE)
#'
#'   ## distance from NAs
#'   ras_distance <- distance(ras_nas)
#'
#'   spreadProbRas <- ras_distance/minmax(ras_distance)["max",]
#'   spreadProbRas[as.vector(spreadProbRas[]) == 0] <- NA
#'
#'   ## varying weights
#'   if (varyWeights) {
#'     weights <- c(2, 3, 7)   ## higher values increase aggregation (decrease edge:area)
#'   } else {
#'     weights <- rep(2, 3)
#'   }
#'
#'   ## we are in a case where we want fewer pixels than available
#'   totalPix <- ncell(ras)
#'   totalAvailPix <- sum(!is.na(ras[]))
#'   pixToConvert <- round(totalAvailPix/4)
#'   availToNARatio <- totalAvailPix/totalPix
#'
#'   ## varying no. pixels to convert
#'   if (varyPixels) {
#'     pixToConvert <- round(pixToConvert * c(0.5, 1, 1.5))
#'   }
#'
#'   ## for landscapes with very high availToNARatio we
#'   ## bump the weights to increase aggregation
#'   weights <- weights*(1+availToNARatio*2)
#'
#'   if (!varyWeights & !varyPixels) {
#'     weights <- unique(weights)
#'     pixToConvert <- unique(pixToConvert)
#'   }
#'
#'   tests <- Map(probWeight = weights,   ## not much difference in clumping levels
#'                pixToConvert = pixToConvert,
#'                numStartsDenom = pixToConvert/10,
#'                f = assignPresences,
#'                MoreArgs = list(assignProb = spreadProbRas,
#'                                landscape = ras))
#'   if (length(tests) == 1) {
#'     tests <- tests[[1]]
#'   } else {
#'     tests <- rast(tests)
#'   }
#'   names(tests) <- paste(paste("Prob. weight", weights),
#'                         paste("no pix", pixToConvert),
#'                         sep = ";")
#'   return(tests)
#' }
#'
#' @return a `SpatRaster` object
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

  if (probWeight < 0.5)
    warning("probWeight < 0.5 may create very fragmented patterns")
  if (probWeight > 7)
    warning("probWeight > 7 may create very clumped patterns")

  ## save original probabilities for later
  assignProbOrig <- as.vector(assignProb[])

  convertedPix <- 0

  ## exponentiate probabilities to provide more weight to pixels further from edges
  assignProbEx <- assignProb^probWeight
  # terra::plot(assignProbEx, col = viridis::inferno(100))

  numStarts <- round(pixToConvert/numStartsDenom)

  while (convertedPix < pixToConvert) {
    ## if the mean is too high, then bring it down to 0.35 to avoid creating
    ## square patches
    meanSP <- mean(as.vector(assignProbEx[]), na.rm = TRUE)
    if (meanSP > 0.35)
      assignProbEx <- assignProbEx / meanSP * 0.35

    ## try to spread in from many focal pixels
    ## draw starting points from areas furthest from edge
    probs <- as.vector(assignProbEx[])
    probs[is.na(probs)] <- 0
    startPoints <- sample(1:ncell(assignProbEx), size = numStarts, prob = probs)

    maxSize <- numStarts ^ probWeight

    # temp <- assignProbEx
    # temp[] <- NA_integer_
    # temp[startPoints] <- 1L
    # terra::plot(assignProbEx, col = viridis::inferno(100))
    # terra::plot(temp, add = TRUE, col = "blue")

    assignProbEx[as.vector(assignProbEx[]) > 1] <- 1  ## otherwise spread throws an error (Error in sample.int(length(x), ...) : invalid 'size' argument)

    outRas <- SpaDES.tools::spread(landscape = raster(landscape),
                                   loci = startPoints,
                                   spreadProb = raster(assignProbEx),
                                   maxSize = maxSize)
    outRas[outRas == 0] <- NA   ## 0s mean no successful spread
    outRas <- rast(outRas) |>
      mask(mask = landscape)

    outRas[!is.na(as.vector(outRas[]))] <- 1L
    # terra::plot(outRas, col = viridis::inferno(100))
    convertedPix <- sum(as.vector(outRas[]), na.rm = TRUE)

    ## increase spread probabilities in case we need to try again
    assignProbEx <- min(assignProbEx * 1.5, 1)
    ## also increase no starting points
    numStarts <- numStarts + 10
  }

  ## if we spread too much remove pixels that are closest to edges
  if (convertedPix > pixToConvert) {
    ## "convert" to distances
    DT <- data.table(cells = 1:ncell(outRas), dists = assignProbOrig)
    DT <- DT[dists > 0 & !is.na(dists)][order(dists, decreasing = TRUE)]

    # pixToKeep <- DT[1:pixToConvert, cells] ## creates concave patches
    pixToKeep <- sample(DT$cells, pixToConvert, prob = DT$dists^probWeight)  ## creates a little more noise

    outRas[] <- NA_integer_
    outRas[pixToKeep] <- 1L
  }
  convertedPix <- sum(!is.na(outRas[]))
  message(convertedPix, " were converted")
  return(outRas)
}


.thermPercent <- function(thermLevel = c("low", "medium", "high"), useMidpoint) {
  thermLevel <- match.arg(thermLevel)
  switch(thermLevel,
         "low" = ifelse(useMidpoint, mean(c(0, 33)), round(runif(1, 0, 33))),
         "medium" = ifelse(useMidpoint, mean(c(34, 66)), round(runif(1, 34, 66))),
         "high" = ifelse(useMidpoint, mean(c(67, 100)), round(runif(1, 67, 100))))
}

.thermWeight <- function(thermLevel = c("low", "medium", "high")) {
  ## higher weights increase the level of aggregation, and thus are more suited
  ## for lower thermokarst levels
  thermLevel <- match.arg(thermLevel)
  switch(thermLevel,
         "low" = 4,
         "medium" = 3,
         "high" = 2)
}

