## TODO: Ceres: these functions can be improved by using options(reproducible...),
## but for now I didn't want to assume those options were directly affecting
## these objects

#' Set NA value in Raster/SpatRaster
#'
#' @param ras a Raster* , or SpatRaster object
#'
#' @importFrom raster NAvalue
#' @importFrom terra NAflag
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


#' Make stacked raster
#'
#' @param ras a Raster* or SpatRaster object
#'
#' @importFrom raster filename
#' @importFrom terra sources
#'
#' @return the filename or source of raster. See `raster::filename` and `terra::sources`
.filename <- function(ras) {
  if (is(ras, "Raster")) {
    return(filename(ras))
  } else {
    if (is(ras, "SpatRaster")) {
      return(sources(ras))
    } else {
      stop("ras should be a Raster* or SpatRaster")
    }
  }
}


#' Project raster extent
#'
#' @param ras a Raster* or SpatRaster object
#'
#' @importFrom raster filename
#' @importFrom terra sources
#'
#' @return the projected extent
.projectExtent <- function(ras, crs) {
  if (is(ras, "Raster")) {
    return(projectExtent(ras, crs = crs))
  } else {
    if (is(ras, "SpatRaster")) {
      rasProjExt <- rast(extent = ext(ras), crs = crs(ras))
      rasProjExt <- project(rasProjExt, y = crs)
      return(ext(rasProjExt))
    } else {
      stop("ras should be a Raster* or SpatRaster")
    }
  }
}
