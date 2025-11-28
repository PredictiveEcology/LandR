#' `LandR` package
#'
#' Utilities for 'LandR' suite of landscape simulation models.
#' These models simulate forest vegetation dynamics based on LANDIS-II, and incorporate
#' fire and insect disturbance, as well as other important ecological processes.
#' Models are implemented as 'SpaDES' modules.
#'
#' @section Package options:
#'
#' `LandR` packages use the following [options()] to configure behaviour:
#'
#' - `LandR.assertions`: If `TRUE`, additional code checks are run during function calls.
#'   Default `FALSE`.
#'
#' @import ggplot2
#' @import patchwork
#' @import methods
#' @name LandR-package
#' @rdname LandR-package
"_PACKAGE"

## usethis namespace: start
#' @importFrom cli col_blue col_cyan col_green col_magenta col_red
#' @importFrom data.table as.data.table copy data.table dcast fifelse fread
#' @importFrom data.table is.data.table last melt rbindlist
#' @importFrom data.table set setattr setcolorder setDT setDTthreads
#' @importFrom data.table setkey setkeyv setnames setorderv
#' @importFrom Formula Formula
#' @importFrom fpCompare %==% %>>% %<<% %<=%
#' @importFrom ggspatial annotation_north_arrow layer_spatial north_arrow_minimal
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom httr2 request
#' @importFrom lme4 glmer lmer
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom parallel mclapply
#' @importFrom pemisc factorValues2 termsInData
#' @importFrom quickPlot layerNames numLayers Plot setColors setColors<-
#' @importFrom raster calc deratify dropLayer extension levels NAvalue<- projectExtent
#' @importFrom raster raster rasterOptions ratify reclassify stack unstack
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom reproducible .prefix .requireNamespace .sortDotsUnderscoreFirst .suffix
#' @importFrom reproducible asPath basename2 Cache CacheDigest cropInputs Filenames fixErrors
#' @importFrom reproducible maxFn messageDF minFn paddedFloatToChar
#' @importFrom reproducible postProcess postProcessTo postProcessTerra
#' @importFrom reproducible prepInputs preProcess projectInputs
#' @importFrom reproducible rasterRead
#' @importFrom reproducible writeOutputs writeTo
#' @importFrom sf as_Spatial st_as_sf st_cast st_coordinates st_intersects st_crs
#' @importFrom sf st_read st_transform st_union st_zm
#' @importFrom sp CRS proj4string SpatialPoints
#' @importFrom SpaDES.tools inRange neutralLandscapeMap randomPolygons rasterizeReduced runifC
#' @importFrom SpaDES.tools spread2
#' @importFrom stats approx as.formula complete.cases fitted glm na.omit
#' @importFrom stats predict quantile runif setNames terms update vcov
#' @importFrom terra app as.factor as.int cellFromRowCol cellFromXY classify coltab<- compareGeom
#' @importFrom terra crop crosstab crs crs<- deepcopy ext extract focalMat
#' @importFrom terra intersect is.factor is.int is.points is.valid
#' @importFrom terra levels mask minmax NAflag<- ncell nlyr project
#' @importFrom terra rast rasterize res rowColFromCell set.names set.values terraOptions
#' @importFrom terra values vect writeRaster xmax xmin ymax ymin xyFromCell
#' @importFrom tidyterra geom_spatraster
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils capture.output combn count.fields data getFromNamespace head install.packages
#' @importFrom utils str tail untar
#' @importFrom viridis scale_fill_viridis
## usethis namespace: end
NULL
