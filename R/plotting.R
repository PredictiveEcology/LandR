utils::globalVariables(c(
  "band1"
))

#' Summary plots of leading vegetation types
#'
#' Create raster of leading vegetation types and \code{Plot} a bar chart summary
#' and a vegetation type map. NOTE: plot order will follow \code{colors} order.
#'
#' @param speciesStack A \code{RasterStack} of percent-cover-by-species layers.
#'
#' @param vtm An optional vegetation type map (\code{RasterLayer}).
#'            If not supplied, will be produced internally by \code{makeVegTypeMap}.
#'
#' @param vegLeadingProportion The minimum proportion cover required to consider
#'                             a species to be the "leading" one. Default 0.8.
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
#' @importFrom data.table data.table setkeyv
#' @importFrom ggplot2 aes element_blank element_text geom_bar ggplot guide_legend guides
#' @importFrom ggplot2 scale_fill_manual scale_x_discrete theme
#' @importFrom pemisc factorValues2
#' @importFrom quickPlot Plot setColors<-
#' @importFrom raster factorValues maxValue minValue
#' @importFrom reproducible Cache
#' @importFrom stats na.omit
plotVTM <- function(speciesStack = NULL, vtm = NULL, vegLeadingProportion = 0.8,
                    sppEquiv, sppEquivCol, colors, title = "Leading vegetation types") {
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
    if (!is.null(speciesStack))
      vtm <- Cache(vegTypeMapGenerator,
                   x = speciesStack,
                   vegLeadingProportion = vegLeadingProportion,
                   mixedType = 2,
                   sppEquiv = sppEquiv,
                   sppEquivCol = sppEquivCol,
                   colors = colors,
                   doAssertion = getOption("LandR.assertions", TRUE))
    else
      stop("plotVTM requires either a speciesStack of percent cover or a",
           " vegetation type map (vtm).")
  }

  ## the ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]
  facLevels <- raster::levels(vtm)[[1]]
  vtmTypes <- as.character(factorValues2(vtm, facLevels$ID, att = 2)) ## 'species', 'Species', 'VALUE'
  vtmCols <- colors[match(vtmTypes, names(colors))]
  whMixed <- which(vtmTypes == "Mixed")

  vtmTypes <- equivalentName(vtmTypes, sppEquiv, "EN_generic_short")
  vtmTypes[whMixed] <- "Mixed"
  names(vtmCols) <- vtmTypes
  facLevels$Species <- vtmTypes #nolint

  ## plot initial types bar chart
  facVals <- factorValues2(vtm, vtm[], att = 2, na.rm = TRUE) ## 'species', 'Species', 'VALUE'
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
  df$species <- factor(df$species, levels = colDT$species, ordered = FALSE)

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
  setColors(vtm, length(vtmTypes)) <- vtmCols ## setColors for factors must have an
                                              ## entry for each row in raster::levels

  Plot(vtm, title = title)
}

#' Create species colour vector from a \code{sppEquiv} table
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#' @param newVals An optional character vector of extra names to use, e.g., "Mixed".
#' @param palette An \pkg{RColorBrewer} palette, e.g., "Accent".
#'     Can get \pkg{RColorBrewer} palette names from \code{rownames(RColorBrewer::brewer.pal.info)}.
#'
#' @return A named vector of colour codes, where the names are the species names
#' plus any extra names passed with \code{newVals}.
#'
#' @export
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom stats na.omit
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

#' @importFrom ggplot2 ggplot scale_fill_distiller theme labs stat unit
#' @importFrom ggspatial layer_spatial annotation_north_arrow north_arrow_minimal
#' @importFrom ggpubr theme_pubr
plotFunction <- function(ras, studyArea, limits = NULL) {
  if (is.null(limits))
    limits <- range(getValues(ras), na.rm = TRUE)
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
