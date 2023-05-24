utils::globalVariables(c(
  "variable", "value"
))

## TODO: Ceres: these functions can be improved by using options(reproducible...),
## but for now I didn't want to assume those options were directly affecting
## these objects

#' Set NA value in Raster/SpatRaster
#'
#' @param ras a Raster*, or SpatRaster object
#' @param NAval the value to use as `NA`
#'
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
    } else {
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


#' @param ras a Raster*, or SpatRaster object
#' @param tab matrix of values to reclassify. See `terra::classify`
#'   and `raster::reclassify`.
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


#' GENERIC EXTRACT POINTS
#'
#' Extracts points from raster layers using the
#'   original raster layer projection.
#'
#' @param x a raster or polygon object (`sp`, `raster` or `terra`)
#' @param y a points or polygons spatial object (`sp`, `sf`, or `terra`)
#' @param field character. The field(s) to extract when x is a polygon.
#'   If `NULL`, all fields are extracted and returned. IDs of y are always
#'   returned (`ID` column).
#' @param ... passed to `terra::extract`
#'
#' @return a data.table with extracted values and an `ID` column of
#'   y point IDs
#'
#' @details
#'  If `x` and `y` are both polygons, `extract` often outputs `NA`
#'  due to misalignments (this can happen even when x == y), even after
#'  `snap(y, x)`. To circumvent this problem, `intersect` is used internally
#'  and, if the `extract` argument `fun` is passed, it is applied to values
#'  of `y` per polygon ID of `x`
#'
#'
#' @importFrom terra rast vect project extract is.points intersect
#' @importFrom data.table as.data.table
#' @export
genericExtract <- function(x, y, field = NULL, ...) {
  if (inherits(x, c("SpatVector", "Spatial", "sf", "sfc"))) {
    if (!is(y, c("SpatVector"))) {
      x <- vect(x)
    }
  } else {
    if (!inherits(x, c("SpatRaster", "Raster"))) {
      stop("x must be a `terra`/`raster`/`sp` raster or polygon")
    } else {
      if (!is(x, "SpatRaster")) {
        x <- rast(x)
      }
    }
  }

  if (!inherits(y, c("SpatVector", "Spatial", "sf", "sfc"))) {
    stop("y must be `terra`/`raster`/ `sf` points/polygons")
  }
  if (!is(y, c("SpatVector"))) {
    y <- vect(y)
  }

  if (!.compareCRS(y, x)) {
    y <- project(y, y = crs(x, proj = TRUE))
  }

  if (is(x, "SpatVector")) {
    if (is.points(y)) {
      out <- terra::extract(x = x, y = y)
    } else {
      ## Even when polygons should be the same, there may be small
      ## misalignments that generate NAs (even after snapping)
      names(y) <- toupper(names(y))
      if (!"ID" %in% names(y)) {
        y$ID <- 1:length(y)
      }
      out <- terra::intersect(x = x, y = y)
      out <- as.data.table(out)

      dots <- list(...)
      if ("fun" %in% names(dots)) {
        if (is.null(field)) {
          field <- names(x)
        }
        out <- out[, lapply(.SD, get(dots$fun), ...), .SDcols = field, by = "ID"]
      }
    }
  } else {
    out <- terra::extract(x = x, y = y, ...)
  }
  out <- as.data.table(out)

  if (!is.null(field)) {
    if (all(field %in% names(x))) {
      IDcol <- names(out)[1]
      out <- out[, .SD, .SDcols = c(IDcol, field)]
      setnames(out, IDcol, "ID")
    } else {
      stop(paste(field, "not found in x"))
    }
  }

  out
}


#' PLOT RASTER FUNCTION
#'
#' A function that creates a `ggplot` of a raster
#'   layer. Can be used with `SpaDES.core::Plots`
#'
#' @param x `SpatRaster`, `RasterLayer`, `SpatVector` or `sf` object
#'
#' @param plotTitle character. A title for the plot passed to
#'   `ggplot::labs(title = plotTitle)`.
#'
#' @param limits TODO
#'
#' @param field character. If x is `sf` or `SpatVector`, a field
#'   to plot.
#'
#' @importFrom ggplot2 geom_raster scale_fill_viridis_c coord_equal theme_classic
#' @importFrom ggplot2 geom_sf coord_sf scale_fill_viridis_d facet_wrap labs aes
#' @importFrom terra nlyr rast vect ncell is.factor
#' @export
plotSpatial <- function(x, plotTitle, limits = NULL, field = NULL) {

  if (inherits(x, c("Raster", "SpatRaster"))) {
    plotRaster <- TRUE
  } else if (inherits(x, c("SpatVector", "sf"))) {
    plotRaster <- FALSE
  } else {
    stop("x must be a SpatRaster, RasterLayer, SpatVector or sf object")
  }

  if (inherits(x, c("Raster"))) {
    x <- rast(x)
    if (is.null(limits)) {
      limits <- range(as.vector(x[]))
    }
  }

  if (inherits(x, "SpatVector")) {
    x <- st_as_sf(x)
    if (is.null(limits) & !is.null(field)) {
      limits <- range(as.vector(x[, field]))
    }
  }

  if (plotRaster) {
    if (!requireNamespace("rasterVis")) {
      stop("Please install 'rasterVis'")
    }
    maxpixels <- ncell(x)
    if (maxpixels > 1E6) {
      message("Raster to plot has >1E6 pixels -- 'maxpixels' set to 1E6 to avoid overly long plotting times")
      maxpixels <- 1E6
    }
    plot1 <- rasterVis::gplot(x, maxpixels = maxpixels)

    if (terra::is.factor(x)) {
      ## need to do this manually for NAs
      cols <- rasVals <- unique(as.vector(x[]))
      cols[!is.na(cols)] <- viridis::cividis(sum(!is.na(cols)))
      names(cols) <- rasVals
      cols[is.na(rasVals)] <- "grey"

      plot1 <- plot1 +
        geom_raster(aes(fill = as.character(value))) +
        scale_fill_manual(values = cols, guide = "legend", limits = limits)
    } else {
      plot1 <- plot1 +
        geom_raster(aes(fill = value)) +
        scale_fill_viridis_c(option = "cividis", na.value = "grey", limits = limits)
    }
    plot1 <- plot1 + coord_equal() +
      theme_classic()

  } else {
    if (!requireNamespace("rlang")) {
      stop("Please install 'rlang'")
    }
    plot1 <- ggplot() + theme_classic()  ## theme needs to come before fill scale because it overrides it in geom_sf

    if (is.null(field)) {
      plot1 <- plot1 +
        geom_sf(data = x)
    } else {
      plot1 <- plot1 +
        geom_sf(data = x, aes(fill = !!sym(field)))
       if (inherits(x[[field]], c("character", "factor"))) {
         plot1 <- plot1 +
           scale_fill_viridis_d(option = "cividis", na.value = "grey", limits = limits)
       }  else {
         plot1 <- plot1 +
           scale_fill_viridis_c(option = "cividis", na.value = "grey", limits = limits)
       }
    }
    plot1 <- plot1 +
      coord_sf()
  }

  if (plotRaster) {
    if (nlyr(x) > 1) {
      plot1 <- plot1 + facet_wrap(~ variable)
      plotTitle <- paste0(plotTitle, " -- ", paste(names(x), collapse = " and "))
    }
  }

  plot1 <- plot1 +
    labs(title = plotTitle, x = "Longitude", y = "Latitude",
         fill = "")
  plot1
}
