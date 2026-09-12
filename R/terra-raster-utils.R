utils::globalVariables(c(
  "value", "variable"
))

## TODO: Ceres: these functions can be improved by using options(reproducible...),
## but for now I didn't want to assume those options were directly affecting
## these objects

#' Set NA values in `Raster` or `SpatRaster`
#'
#' @param ras a `Raster`, or `SpatRaster` object
#'
#' @param NAval the value to use as `NA`
#'
#' @return a raster with attributed NA values
#'
#' @rdname rasterTerraHelpers
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
#' @param rasList a list of `Raster` or `SpatRaster` objects
#'
#' @return a stacked raster
#'
#' @export
#' @rdname rasterTerraHelpers
.stack <- function(rasList) {
  isRaster <- sapply(rasList, function(ras) is(ras, "Raster"))
  isSpatRas <- sapply(rasList, function(ras) is(ras, "SpatRaster"))
  if (all(isRaster)) {
    out <- raster::stack(rasList)
    names(out) <- names(rasList)
  } else if (all(isSpatRas)) {
    out <- rast(rasList)
    names(out) <- names(rasList)
  } else {
    stop("List entries should all be RasterLayer or SpatRaster")
  }

  return(out)
}

#' Project raster extent
#'
#' @param ras a `Raster` or `SpatRaster` object
#'
#' @param crs passed to `raster::projectRaster(..., crs = crs)` and `terra::project(..., y = crs)`
#'
#' @return the projected extent
#'
#' @rdname rasterTerraHelpers
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
#' @note this function internally converts `Raster` to `SpatRaster` to allow using `compareGeom()`,
#' and benefit from its complexity.
#'
#' @param ... additional `Raster` or `SpatRaster` objects, and arguments
#'            passed to [terra::compareGeom()].
#'
#' @return the projected extent
#'
#' @export
#' @rdname compare
.compareRas <- function(x, ...) {
  mc <- match.call()
  dots <- list(...) # need to pass the ...

  ## subset spatial and non-spatial arguments
  spatialDots <- vapply(dots,
    FUN = inherits,
    what = c("Raster", "SpatRaster", "Spatial", "sf", "SpatVector"),
    FUN.VALUE = logical(1)
  )

  objs <- append(list(x), dots[spatialDots])

  otherArgs <- dots[!spatialDots]

  if (is.null(otherArgs$crs)) {
    otherArgs$crs <- TRUE
  }

  if (is.null(otherArgs$ext)) {
    otherArgs$ext <- TRUE
  }

  if (is.null(otherArgs$stopOnError)) {
    otherArgs$stopOnError <- TRUE
  }


  whRast <- vapply(objs, inherits, c("Raster", "SpatRaster"), FUN.VALUE = logical(1))

  if (any(whRast)) {
    objs <- Map(function(ras) {
      if (is(ras, "Raster")) {
        rast(ras)
      } else {
        ras
      }
    }, ras = objs)
  }

  whSpatialSf <- vapply(objs, inherits, c("Spatial", "sf"), FUN.VALUE = logical(1))
  if (any(whSpatialSf)) {
    objs <- Map(function(vec) {
      if (inherits(vec, c("Spatial", "sf"))) {
        vect(vec)
      } else {
        vec
      }
    }, vec = objs)
  }

  x <- objs[[1]]

  for (i in 2:length(objs)) {
    out <- TRUE
    if (isTRUE(otherArgs$crs)) {
      out <- .compareCRS(x, objs[[i]])
      otherArgs$crs <- FALSE ## prevent re-checking below.
    }

    if (!isTRUE(out)) {
      if (isTRUE(otherArgs$stopOnError)) {
        stop(".compareCRS fail: ", format(mc[[i + 1]]), " is not same as ", format(mc[["x"]]))
      }
      return(out)
    } else {
      if (is(x, "SpatRaster") && is(objs[[i]], "SpatRaster")) {
        ## note: compareGeom is failing even with 2 spatvect
        out <- do.call(compareGeom, append(list(x, objs[[i]]), otherArgs))
      } else {
        if (isTRUE(otherArgs$ext)) {
          out <- ext(x) == ext(objs[[i]])

          if (isFALSE(out) && isTRUE(otherArgs$stopOnError)) {
            stop(
              ".compareRas fail: ", format(mc[[i + 1]]), " and ",
              format(mc[["x"]]), " have different extents."
            )
          }
        }
      }
      if (!isTRUE(out)) {
        message(
          "compareGeom/.compareRas fail: ", format(mc[[i + 1]]),
          " is not same as ", format(mc[["x"]])
        )
        return(out)
      }
    }
  }
  return(out)
}

#' Compare two rasters' projections
#'
#' TODO: Move to `reproducible`
#' TODO: expand to multiple objects
#'
#' @param x,y a `Raster`, `SpatRaster`, `sf`, `SpatVector`, or `Spatial` object
#'
#' @export
#' @rdname compare
.compareCRS <- function(x, y) {
  st_crs(x) == st_crs(y)
}

#' Helpers for transition to `terra`
#'
#' These all create a single function that can be used for either `Raster` or `SpatRaster` objects.
#'
#' @return
#' `asInt` returns a `*Raster` with values converted to `integer`, if they weren't already.
#'
#' @export
#' @rdname rasterTerraHelpers
asInt <- function(ras) {
  if (!isInt(ras)) {
    if (inherits(ras, "SpatRaster")) {
      ras <- as.int(ras)
    } else {
      ras[] <- as.integer(ras)
    }
  }
  ras
}

#' @export
#' @rdname rasterTerraHelpers
#' @return
#' `isInt` returns a logical as per `is.integer`.
isInt <- function(ras) {
  if (inherits(ras, "SpatRaster")) {
    is.int(ras)
  } else {
    is.integer(values(ras))
  }
}

#' @param ras a `Raster`, or `SpatRaster` object
#'
#' @param tab matrix of values to reclassify. See `terra::classify` and `raster::reclassify`.
#'
#' @return
#' `reclass` returns a `*Raster` with values reclassified as per `terra::classify`
#' and `raster::reclassify`.
#'
#' @export
#' @rdname rasterTerraHelpers
reclass <- function(ras, tab) {
  if (is(ras, "SpatRaster")) {
    classify(ras, tab)
  } else {
    reclassify(ras, tab)
  }
}

#' Generic extract points
#'
#' Extracts points from raster layers using the original raster layer projection.
#'
#' @details
#' If `x` and `y` are both polygons, `extract` often outputs `NA` due to misalignments
#' (this can happen even when `x == y`), even after `snap(y, x)`.
#' To circumvent this problem, `intersect` is used internally and,
#' if the `extract` argument `fun` is passed, it is applied to values of `y` per polygon ID of `x`.
#'
#' @param x a raster or polygon object (`sp`, `raster` or `terra`)
#'
#' @param y a points or polygons spatial object (`sp`, `sf`, or `terra`)
#'
#' @param field character. The field(s) to extract when x is a polygon.
#' If `NULL`, all fields are extracted and returned. IDs of y are always returned (`ID` column).
#'
#' @param ... passed to [terra::extract()]
#'
#' @return a `data.table` with extracted values and an `ID` column of y point IDs
#'
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
        y$ID <- seq_along(y)
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

#' Create a `ggplot` of a raster or `sf` object.
#'
#' Can be used with `SpaDES.core::Plots`.
#'
#' @param x `SpatRaster`, `RasterLayer`, `SpatVector` or `sf` object
#'
#' @param plotTitle character. A title for the plot passed to `ggplot::labs(title = plotTitle)`.
#'
#' @param limits TODO
#'
#' @param field character. If `x` is `sf` or `SpatVector`, a field to plot.
#'
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
    if (is.null(limits) && !is.null(field)) {
      limits <- range(as.vector(x[, field]))
    }
  }

  if (plotRaster) {
    if (!requireNamespace("rasterVis", quietly = TRUE)) {
      stop("Please install 'rasterVis'")
    }
    maxpixels <- ncell(x)
    if (maxpixels > 1e6) {
      message(
        "Raster to plot has >1e6 pixels. ",
        "'maxpixels' set to 1e6 to avoid overly long plotting times."
      )
      maxpixels <- 1e6
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
    if (!requireNamespace("rlang", quietly = TRUE)) {
      stop("Please install 'rlang'")
    }
    ## theme needs to come before fill scale because it overrides it in geom_sf
    plot1 <- ggplot() +
      theme_classic()

    if (is.null(field)) {
      plot1 <- plot1 +
        geom_sf(data = x)
    } else {
      plot1 <- plot1 +
        geom_sf(data = x, aes(fill = !!sym(field)))
      if (inherits(x[[field]], c("character", "factor"))) {
        plot1 <- plot1 +
          scale_fill_viridis_d(option = "cividis", na.value = "grey", limits = limits)
      } else {
        plot1 <- plot1 +
          scale_fill_viridis_c(option = "cividis", na.value = "grey", limits = limits)
      }
    }
    plot1 <- plot1 +
      coord_sf()
  }

  if (plotRaster) {
    if (nlyr(x) > 1) {
      plot1 <- plot1 + facet_wrap(~variable)
      plotTitle <- paste0(plotTitle, " -- ", paste(names(x), collapse = " and "))
    }
  }

  plot1 <- plot1 +
    labs(
      title = plotTitle, x = "Longitude", y = "Latitude",
      fill = ""
    )
  plot1
}
