utils::globalVariables(c(
  "band1"
))

#' Summary plots of leading vegetation types
#'
#' Create raster of leading vegetation types and `Plot` a bar chart summary
#' and a vegetation type map. NOTE: plot order will follow `colors` order.
#'
#' @param speciesStack A `SpatRaster`, `RasterStack` or `RasterBrick`
#'   of percent-cover-by-species layers.
#'
#' @param vtm An optional vegetation type map (`RasterLayer` or `SpatRaster`).
#'            If not supplied, will be produced internally by `makeVegTypeMap`.
#'
#' @template vegLeadingProportion
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param colors Named vector of colour codes, named using species names. NOTE:
#'               plot order will follow this order.
#'
#' @param title The title to use for the generated plots.
#'
#' @author Eliot McIntire
#' @export
plotVTM <- function(speciesStack = NULL, vtm = NULL, vegLeadingProportion = 0.8,
                    sppEquiv, sppEquivCol, colors, title = "Leading vegetation types") {

  if (is(speciesStack, "RasterBrick")) {
    speciesStack <- raster::stack(speciesStack)
  }

  colorsEN <- equivalentName(names(colors), sppEquiv, "EN_generic_short")
  colDT <- data.table(cols = colors, species = colorsEN,
                      speciesOrig = names(colors),
                      speciesOrigOrder = seq(colors))
  mixedString <- "Mixed"
  hasMixed <- isTRUE(mixedString %in% names(colors))
  if (hasMixed) {
    whMixedColors <- which(names(colors) == mixedString)
    colDT[whMixedColors, species := mixedString]
  }

  setkeyv(colDT, "speciesOrigOrder")

  newStackOrder <- na.omit(match(colDT$speciesOrig, names(speciesStack)))
  speciesStack <- speciesStack[[newStackOrder]]

  if (is.null(vtm)) {
    if (!is.null(speciesStack)) {
      vtm <- Cache(vegTypeMapGenerator,
                   x = speciesStack,
                   vegLeadingProportion = vegLeadingProportion,
                   mixedType = 2,
                   sppEquiv = sppEquiv,
                   sppEquivCol = sppEquivCol,
                   colors = colors,
                   doAssertion = getOption("LandR.assertions", TRUE))
    } else
      stop("plotVTM requires either a speciesStack of percent cover or a",
           " vegetation type map (vtm).")
  }

  ## the ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]
  if (is(vtm, "RasterLayer")) {
    facLevels <- raster::levels(vtm)[[1]]
  } else {
    facLevels <- levels(vtm)[[1]]
  }
  vtmTypes <- as.character(factorValues2(vtm, na.omit(as.vector(unique(vtm[]))), att = 2)) ## 'species', 'Species', 'VALUE'
  vtmCols <- colors[match(vtmTypes, names(colors))]
  whMixed <- which(vtmTypes == "Mixed")

  vtmTypes <- equivalentName(vtmTypes, sppEquiv, "EN_generic_short")
  vtmTypes[whMixed] <- "Mixed"
  names(vtmCols) <- vtmTypes
  facLevels$Species <- vtmTypes #nolint

  ## plot initial types bar chart
  facVals <- factorValues2(vtm, as.vector(vtm[]), att = 2, na.rm = TRUE) ## 'species', 'Species', 'VALUE'
  df <- data.table(species = as.character(facVals), stringsAsFactors = FALSE)
  df <- df[!is.na(df$species)]

  speciesEN <- equivalentName(df$species, sppEquiv, "EN_generic_short")
  if (all(na.omit(speciesEN) %in% colorsEN)) {
    whMixed <- which(df$species == mixedString)

    df$species <- speciesEN

    if (hasMixed)
      df[whMixed, species := mixedString]

    df <- colDT[df, on = "species"] # merge color and species
  } else {
    stop("Species names of 'colors' must match those in 'speciesStack'.")
  }

  # Needs to be factor so ggplot2 knows that there may be missing levels
  df$species <- factor(df$species, levels = unique(colDT$species), ordered = FALSE)

  cols2 <- colDT$cols
  names(cols2) <- colDT$species

  initialLeadingPlot <- ggplot(data = df, aes(species, fill = species)) +
    scale_x_discrete(drop = FALSE) +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = cols2, drop = FALSE) +
    geom_bar(position = "stack") +
    theme(legend.text = element_text(size = 6), legend.title = element_blank(),
          axis.text = element_text(size = 6))

  Plot(initialLeadingPlot, title = title)

  ## plot inital types raster
  levels(vtm) <- facLevels
  if (is(vtm, "RasterLayer")) {
    setColors(vtm, length(vtmTypes)) <- vtmCols ## setColors for factors must have an
    ## entry for each row in raster::levels
  } else {
    ## TODO: setColors needs to be adapted to SpatRaster...
  }

  cols2 <- colDT$cols
  names(cols2) <- colDT$speciesOrig

  labs <- colDT$species
  names(labs) <- colDT$speciesOrig

  vtmPlot <- if (is(vtm, "RasterLayer")) {
     ggplot() + geom_raster(data = vtm)
  } else {
    ggplot() + geom_spatraster(data = vtm)
  }
  vtmPlot <- vtmPlot +
    scale_fill_manual(values = cols2,  labels = labs,
                        na.value = "grey80") +
    theme(legend.text = element_text(size = 6), legend.title = element_blank(),
          axis.text = element_text(size = 6))

  Plot(vtmPlot, title = title)
}

#' Helper for setting Raster or `SpatRaster` colors
#'
#' This is a wrapper to help with migration to \pkg{terra}.
#' Currently can only be used for a single layer `SpatRaster` or a `RasterLayer`.
#'
#' @param ras A `Raster*` or `SpatRaster` class object.
#'
#' @param cols a character vector of colours. See examples. Can also be a `data.frame`,
#'   see `terra::coltab`
#'
#' @param n A numeric scalar giving the number of colours to create. Passed to
#'   `quickPlot::setColors(ras, n = n) <- `. If missing, then `n` will be `length(cols)`
#'
#' @examples
#' \donttest{
#' cols <- colorRampPalette(c("blue", "red"))(12)
#' ras <- terra::rast(matrix(1:100, 10, 10))
#' ras <- Colors(ras, cols)
#' terra::plot(ras)
#'
#' ras <- raster::raster(matrix(1:100, 10, 10))
#' ras <- Colors(ras, cols)
#' raster::plot(ras)
#' }
#'
#' @export
Colors <- function(ras, cols, n = NULL) {
  if (is(ras, "SpatRaster")) {
    theSeq <- round(minFn(ras)):round(maxFn(ras))
    if (!is(cols, "data.frame")) {
      if (length(cols) < length(theSeq)) {
        message("not enough colours, interpolating")
        cols <- colorRampPalette(cols)(length(theSeq) + 1)
      }
    }
    coltab(ras, layer = 1) <- cols
  } else {
    if (is.null(n)) {
      n <- length(cols)
    }
    setColors(ras, n = n) <- cols
  }

  ras
}

#' Create species colour vector from a `sppEquiv` table
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#' @param newVals An optional character vector of extra names to use, e.g., "Mixed".
#' @param palette An \pkg{RColorBrewer} palette, e.g., "Accent".
#'     Can get \pkg{RColorBrewer} palette names from `rownames(RColorBrewer::brewer.pal.info)`.
#'
#' @return A named vector of colour codes, where the names are the species names
#' plus any extra names passed with `newVals`.
#'
#' @export
sppColors <- function(sppEquiv, sppEquivCol, newVals = NULL, palette) {
  sppColorNames <- c(na.omit(unique(sppEquiv[[sppEquivCol]])), newVals)

  sppColors <- NULL
  sppColors <- if (is.character(palette))
    if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      colorPalette <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 7, name = palette))
      colorPalette(length(sppColorNames))
    }

  if (is.null(sppColors))
    stop("Currently palette must be one of the RColorBrewer::brewer.pal names")

  names(sppColors) <- sppColorNames
  sppColors
}

plotFunction <- function(ras, studyArea, limits = NULL) {
  if (is.null(limits))
    limits <- range(as.vector(ras[]), na.rm = TRUE)
  ggplot() +
    layer_spatial(ras, aes(fill = stat(band1))) +
    layer_spatial(data = studyArea, fill = "transparent", colour = "black") +
    annotation_north_arrow(style = north_arrow_minimal,
                           height = unit(1, "cm"), width = unit(1, "cm"),
                           location = "tr", which_north = "true") +
    theme_pubr(legend = "bottom") +
    theme(plot.margin = unit(c(0,0,0,0), units = "mm")) +
    scale_fill_distiller(palette = "Greys", na.value = "transparent",
                         direction = 1,
                         breaks = seq(limits[1], limits[2], length.out = 6),
                         limits = limits) +
    labs(x = "longitude", y = "latitude", fill = "Cover",
         title = sub("\\.|_", " ", names(ras)))
}
