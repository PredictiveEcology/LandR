## TODO: Ceres: these functions can be improved by using options(reproducible...),
## but for now I didn't want to assume those options were directly affecting
## these objects

#' Set NA value in Raster/SpatRaster
#'
#' @param ras a Raster* , or SpatRaster object
#' @param NAval the value to use as `NA`
#'
#' @importFrom raster NAvalue<-
#' @importFrom terra NAflag<-
#'
#' @return a raster with attributed NA values
.NAvalueFlag <- function(ras, NAval) {
  if (is(ras, "Raster")) {
    NAvalue(ras) <- NAval
  } else {
    if (is(ras, "SpatRaster")) {
      NAflag(ras) <- NAval
    } else {
      stop("ras should be a Raster* or SpatRaster")
    }
  }
  ras
}


#' Make stacked raster
#'
#' @param rasList a list of Raster* or SpatRaster layers
#'
#' @importFrom raster stack
#' @importFrom terra rast
#'
#' @return a stacked raster
.stack <- function(rasList) {
  isRaster <- sapply(rasList, function(ras) is(ras, "Raster"))
  isSpatRas <- sapply(rasList, function(ras) is(ras, "SpatRaster"))
  if (all(isRaster)) {
    return(raster::stack(rasList))
  } else {
    if (all(isSpatRas)) {
      return(rast(rasList))
    } else {
      stop("List entries should all be RasterLayer or SpatRaster")
    }
  }
}

#' Project raster extent
#'
#' @param ras a Raster* or SpatRaster object
#' @param crs passed to `raster::projectRaster(..., crs = crs)`
#'   and `terra::project(..., y = crs)`
#'
#' @importFrom raster projectExtent
#' @importFrom terra ext project
#'
#' @return the projected extent
.projectExtent <- function(ras, crs) {
  if (is(ras, "Raster")) {
    return(projectExtent(ras, crs = crs))
  } else {
    if (is(ras, "SpatRaster")) {
      rasProj <- project(ras, y = crs)
      return(ext(rasProj))
    } else {
      stop("ras should be a Raster* or SpatRaster")
    }
  }
}

#' Compare two raster's properties
#'
#' Note: this function internally converts Raster* to
#' `SpatRaster` to allow using `compareGeom` and benefit
#' from its complexity
#' @param ras a Raster* or SpatRaster object
#' @param ras2 a Raster* or SpatRaster object
#' @param ... passed to `terra::compareGeom`
#'
#' @importFrom terra compareGeom rast
#' @rdname compare
#' @export
#' @return the projected extent
#'
## TODO: should be extended to many rasters
.compareRas <- function(ras1, ras2, ...) {
  if (is(ras1, "Raster")) {
    ras1 <- rast(ras1)
  }
  if (is(ras2, "Raster")) {
    ras2 <- rast(ras2)
  }

  .compareCRS(ras1, ras2)

  dots <- list(...) # need to pass the ...
  dots$crs <- FALSE
  do.call(compareGeom, append(list(ras1, ras2), dots))
}

#' @export
#' @rdname compare
#' @importFrom sf st_crs
.compareCRS <- function(ras1, ras2, ...) {
  st_crs(ras1) != st_crs(ras2)
}

#' Method to read raster
#'
#' TODO: Move to `reproducible`
#' @param ... passed to `terra::rast` or
#'   `raster::raster`
#'
#' @return the function
#' @export
rasterRead <- function(...)
  eval(parse(text = getOption("reproducible.rasterRead")))(...)
