## TODO: Ceres: these functions can be improved by using options(reproducible...),
## but for now I didn't want to assume those options were directly affecting
## these objects

#' Set NA value in Raster/SpatRaster
#'
#' @param ras a Raster*, or SpatRaster object
#' @param NAval the value to use as `NA`
#' @return a raster with attributed NA values
#'
#' @rdname rasterTerraHelpers
#' @importFrom raster NAvalue<-
#' @importFrom terra NAflag<-
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
#' @return a stacked raster
#'
#' @rdname rasterTerraHelpers
#' @importFrom raster stack
#' @importFrom terra rast
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
#' @return the projected extent
#'
#' @rdname rasterTerraHelpers
#' @importFrom raster projectExtent
#' @importFrom terra ext project
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

#' Compare raster properties
#'
#' TODO: Move to `reproducible`
#'
#' Note: this function internally converts Raster* to
#' `SpatRaster` to allow using `compareGeom` and benefit
#' from its complexity
#' @param ras1 a Raster* or SpatRaster object
#' @param ... additional Raster* or SpatRaster objects, and arguments
#'   passed to `terra::compareGeom`
#'
#' @importFrom terra compareGeom rast
#' @rdname compare
#' @export
#' @return the projected extent
#'
.compareRas <- function(ras1, ...) {
  mc <- match.call()

  dots <- list(...) # need to pass the ...
  whRast <- vapply(dots, inherits, c("Raster", "SpatRaster"), FUN.VALUE = logical(1))
  rasts <- append(list(ras1), dots[whRast])

  rasts <- Map(function(ras) {
    if (is(ras, "Raster")) {
      rast(ras)
    } else ras
  }, ras = rasts)

  dotsNotRasters <- dots[!whRast]
  dotsNotRasters$crs <- FALSE

  ras1 <- rasts[[1]]

  for (i in 2:length(rasts)) {
    out <- .compareCRS(ras1, rasts[[i]])
    if (!isTRUE(out)) {
      message(".compareCRS fail: ", format(mc[[i + 1]]), " is not same as ", format(mc[["ras1"]]))
      break
    } else{
      out <- do.call(compareGeom, append(list(ras1, rasts[[i]]), dotsNotRasters))
      if (!isTRUE(out)) {
        message(".compareRas fail: ", format(mc[[i + 1]]), " is not same as ", format(mc[["ras1"]]))
        break
      }
    }
  }
  return(out)
}

#' Compare two raster's projections
#'
#' TODO: Move to `reproducible`
#' TODO: expand to multiple objects
#'
#' @param ras1 a Raster* or SpatRaster object
#' @param ras2 a Raster* or SpatRaster object
#' @export
#' @rdname compare
#' @importFrom sf st_crs
.compareCRS <- function(ras1, ras2) {
  st_crs(ras1) == st_crs(ras2)
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



#' Helpers for transition to `terra`
#'
#' These all create a single function that can be used for either `Raster` or `SpatRaster`
#' objects.
#'
#' @export
#' @rdname rasterTerraHelpers
#' @return
#' `asInt` returns a `*Raster` with values converted to `integer`, if they weren't already.
#' @importFrom terra as.int
asInt <- function(ras) {
  if (!isInt(ras)) {
    if (inherits(ras, "SpatRaster"))
      ras <- as.int(ras)
    else
      ras[] <- as.integer(ras)

  }
  ras
}

#' @export
#' @rdname rasterTerraHelpers
#' @return
#' `isInt` returns a logical as per `is.integer`.
#' @importFrom terra is.int
isInt <- function(ras) {
  if (inherits(ras, "SpatRaster"))
    is.int(ras)
  else
    is.integer(values(ras))
}


#' @export
#' @rdname rasterTerraHelpers
#' @return
#' `reclass` returns a `*Raster` with values reclassified as per `terra::classify`
#'   and `raster::reclassify`.
#' @importFrom terra classify
#' @importFrom raster reclassify
reclass <- function(ras, tab) {
  if (is(ras, "SpatRaster")) {
    classify(ras, tab)
  } else {
    reclassify(ras, tab)
  }
}


