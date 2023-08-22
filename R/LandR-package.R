#' `LandR` package
#'
#' Utilities for 'LandR' suite of landscape simulation models.
#' These models simulate forest vegetation dynamics based on LANDIS-II, and
#' incorporate fire and insect disturbance, as well as other important ecological
#' processes. Models are implemented as 'SpaDES' modules.
#'
#' @section Package options:
#'
#' `LandR` packages use the following [options()] to configure behaviour:
#'
#' - `LandR.assertions`: If `TRUE`, additional code checks are run during function calls.
#'   Default `FALSE`.
#'
#' @import methods
#' @name LandR-package
#' @rdname LandR-package
"_PACKAGE"

## usethis namespace: start
#' @importFrom crayon blue cyan green magenta red
#' @importFrom data.table as.data.table copy data.table fread is.data.table last melt rbindlist
#' @importFrom data.table copy
#' @importFrom data.table set setattr setcolorder setDT setDTthreads setkey setkeyv setnames setorderv
#' @importFrom fpCompare %==% %>>% %<<%
#' @importFrom ggplot2 aes coord_equal coord_sf element_blank element_text facet_wrap
#' @importFrom ggplot2 geom_bar geom_hline geom_line geom_point labs geom_raster geom_ribbon geom_sf
#' @importFrom ggplot2 ggplot guide_legend guides labs
#' @importFrom ggplot2 scale_color_distiller
#' @importFrom ggplot2 scale_fill_distiller scale_fill_manual scale_fill_viridis_c scale_fill_viridis_d
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 stat stat_summary sym theme theme_classic unit
#' @importFrom ggpubr theme_pubr
#' @importFrom ggspatial annotation_north_arrow layer_spatial north_arrow_minimal
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom lme4 glmer lmer
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom parallel mclapply
#' @importFrom pemisc factorValues2 termsInData
#' @importFrom quickPlot layerNames numLayers Plot setColors setColors<-
#' @importFrom raster calc deratify dropLayer extension levels NAvalue<- projectExtent
#' @importFrom raster raster rasterOptions ratify reclassify stack unstack
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom reproducible .prefix .requireNamespace .sortDotsUnderscoreFirst .suffix
#' @importFrom reproducible asPath basename2 Cache cropInputs Filenames fixErrors
#' @importFrom reproducible maxFn messageDF minFn paddedFloatToChar
#' @importFrom reproducible postProcess postProcessTerra prepInputs preProcess projectInputs
#' @importFrom reproducible writeOutputs
#' @importFrom sf as_Spatial st_as_sf st_cast st_coordinates  st_intersects st_crs
#' @importFrom sf st_read st_transform st_zm
#' @importFrom sp CRS proj4string SpatialPoints
#' @importFrom SpaDES.tools inRange neutralLandscapeMap randomPolygons rasterizeReduced runifC
#' @importFrom SpaDES.tools spread2
#' @importFrom stats as.formula complete.cases fitted glm na.omit predict quantile runif terms update vcov
#' @importFrom terra app as.int cellFromRowCol cellFromXY classify coltab<- compareGeom crop crs crs<- ext extract focalMat
#' @importFrom terra intersect is.factor is.int is.points levels mask minmax NAflag<- ncell nlyr project
#' @importFrom terra rast rasterize res rowColFromCell terraOptions values vect writeRaster xmax xmin ymax ymin
#' @importFrom terra xyFromCell
#' @importFrom tidyterra geom_spatraster
#' @importFrom tools file_path_sans_ext
#' @importFrom utils capture.output combn count.fields data getFromNamespace install.packages
#' @importFrom utils str tail untar
## usethis namespace: end
NULL
