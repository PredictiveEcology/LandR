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
#' @param sppEquiv a species equivalency table TODO: description needed
#'
#' @param sppEquivCol the column name to use from \code{sppEquiv}.
#'
#' @param colors Named vector of colour codes, named using species names. NOTE:
#'               plot order will follow this order.
#'
#' @param title The title to use for the generated plots.
#'
#' @author Eilot McIntire
#' @export
#' @importFrom data.table data.table setkeyv
#' @importFrom ggplot2 aes element_blank element_text geom_bar ggplot scale_fill_manual theme
#' @importFrom ggplot2 guides guide_legend guide_legend scale_x_discrete
#' @importFrom quickPlot Plot setColors<-
#' @importFrom raster factorValues maxValue minValue
#' @importFrom reproducible Cache
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
      vtm <- Cache(makeVegTypeMap, speciesStack, vegLeadingProportion, mixed = TRUE)
    else
      stop("plotVTM requires either a speciesStack of percent cover or a",
           " vegetation type map (vtm).")
  }

  ## the ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]),]
  facLevels <- raster::levels(vtm)[[1]]
  vtmTypes <- as.character(factorValues2(vtm, facLevels$ID, att = "Species"))
  vtmCols <- colors[match(vtmTypes, names(colors))]
  whMixed <- which(vtmTypes == "Mixed")

  vtmTypes <- equivalentName(vtmTypes, sppEquiv, "EN_generic_short")
  vtmTypes[whMixed] <- "Mixed"
  names(vtmCols) <- vtmTypes
  facLevels$Species <- vtmTypes

  ## plot initial types bar chart
  facVals <- factorValues2(vtm, vtm[], att = "Species", na.rm = TRUE)
  df <- data.table(species = as.character(facVals), stringsAsFactors = FALSE)
  df <- df[!is.na(df$species)]

  speciesEN <- equivalentName(df$species, sppEquiv, "EN_generic_short")
  if (all(na.omit(speciesEN) %in% colorsEN) ){
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
  setColors(vtm, length(vtmTypes)) <- vtmCols # setColors for factors must have an entry for each row in raster::levels

  Plot(vtm, title = title)
}

#' Create species color vector from a sppEquiv table
#'
#' Create species color vector from a sppEquiv table
#'
#' @param sppEquiv A species equivalency table, e.g., \code{data("sppEquivalencies_CA")}.
#' @param sppEquivCol The name of the column to get names from.
#' @param newVals An optional character vector of extra names to use, e.g., "Mixed".
#' @param palette An RColorBrewer palette, e.g., "Accent".
#'                Can get RColorBrewer palette names from
#'                \code{rownames(RColorBrewer::brewer.pal.info)}.
#'
#' @return A named vector of color codes, where the names are the species names
#' plus any extra names passed with \code{newVals}.
#'
#' @export
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
sppColors <- function(sppEquiv, sppEquivCol, newVals, palette) {
  sppColorNames <- c(na.omit(unique(sppEquiv[[sppEquivCol]])), newVals)

  sppColors <- NULL
  sppColors <- if (is.character(palette))
    if (palette %in% rownames(RColorBrewer::brewer.pal.info))
      RColorBrewer::brewer.pal(length(sppColorNames), palette)

  if (is.null(sppColors))
    stop("Currently palette must be one of the RColorBrewer::brewer.pal names")
  names(sppColors) <- sppColorNames
  sppColors
}
