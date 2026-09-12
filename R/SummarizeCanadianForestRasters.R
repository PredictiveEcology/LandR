utils::globalVariables(c("bin", "count", "geometry", "ID", "id_col"))

prep_polygons <- function(raster, polygons, polygon_id = NULL, filter_ids = NULL) {
  if (inherits(polygons, "SpatVector")) {
    polygons <- sf::st_as_sf(polygons)
  }

  ## Don't name the geometry column: it is "geometry" for a shapefile but "geom"
  ## for a GeoPackage, and sf stores the actual name in the `sf_column`
  ## attribute. `summarise()` on a grouped sf dissolves the geometries itself
  ## (`do_union = TRUE`), so there is no need to reference the column at all.
  polygons <- sf::st_transform(polygons, sf::st_crs(raster)) |>
    dplyr::rename(ID = !!polygon_id) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(.groups = "drop")

  if (!is.null(filter_ids)) {
    polygons <- dplyr::filter(polygons, ID %in% filter_ids)
  }

  terra::vect(polygons)
}

#' Raster summary statistics and maps by polygon
#'
#' `calc_raster_stats` iteratively calculates the frequency of raster values within each polygon.
#' Supports user-provided polygons or automatic retrieval of ecoregions/ecozones.
#'
#' `plot_raster_stats` calculates and returns statistics (min, max, median, mean,
#' 25th percentile, 75th percentile, and proportion of zeroes) for a raster,
#' and creates and saves figures for each polygon in a vector layer, including
#' a histogram of raster values, a map of raster values within the polygon, and
#' a summary of statistics . Optionally, zero values can be removed from the
#' analysis, and the raster can be aggregated to a lower resolution to speed up plotting.
#'
#' @note Large rasters should be processed on disk wherever possible (see
#' [terra::terraOptions()]). Cap the RAM terra uses per operation with `memmax`
#' (in GB) and force results to disk with `todisk = TRUE`, e.g.
#' `terraOptions(memmax = 16, todisk = TRUE, tempdir = <fast local dir>)`.
#' Prefer `memmax` over `memfrac`: `memfrac` is a *fraction of free RAM* (not a hard
#' cap), and `memfrac = 0.0` crashes (terra's block size collapses to 0). Point
#' `tempdir` at fast *local* storage (NVMe) -- never NFS, and never a tmpfs such
#' as `/tmp` on some systems, which is RAM-backed (so a raster "copied to disk"
#' there really sits in RAM, e.g. via `withr::local_tempdir()`).
#'
#' For an integer raster, `calc_raster_stats()` can be derived from the
#' value-frequency table produced by `calc_raster_counts()` (every statistic is
#' an exact function of it), which avoids reading every pixel into memory.
#'
#' @param raster A `SpatRaster` object.
#'
#' @param polygons A `SpatVector` or `sf` object of polygons.
#'
#' @param polygon_id Name of the column in `polygons` to use as polygon ID.
#'
#' @param filter_ids Optional. Polygon IDs to filter.
#'
#' @returns A `data.frame` with columns: `ID`, `value`, `count`.
#'
#' @examples
#' if (interactive()) {
#'   if (requireNamespace("dplyr", quietly = TRUE) &&
#'       requireNamespace("geodata", quietly = TRUE) &&
#'       requireNamespace("withr", quietly = TRUE)) {
#'     tmp_pth <- withr::local_tempdir()
#'
#'     ## process the raster on disk, with bounded RAM (see @note)
#'     terraOptions(memmax = 16, todisk = TRUE, tempdir = tmp_pth)
#'
#'     age <- prepInputsStandAgeMap(
#'       dataSource = "KNN",
#'       dataYear = 2011,
#'       destinationPath = tmp_pth
#'     )
#'
#'     ecozones <- reproducible::prepInputs(
#'       url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
#'       destinationPath = tmp_pth
#'     ) |>
#'      sf::st_make_valid()
#'
#'     id_col <- "ZONE_NAME"
#'
#'     ## WARNING: don't use Cache() below -- excessive RAM usage
#'
#'     ecozone_counts <- calc_raster_counts(
#'       raster = age,
#'       polygons = ecozones,
#'       polygon_id = id_col
#'     )
#'
#'     ## reuse the counts to avoid recomputing the frequency table
#'     ecozone_stats <- calc_raster_stats(
#'       raster = age,
#'       polygons = ecozones,
#'       polygon_id = id_col,
#'       counts_df = ecozone_counts
#'     )
#'
#'     gc()
#'
#'     ## below will call calc_raster_counts and calc_raster_stats internally,
#'     ## so we shortcut this by passing the previously-calculated tables
#'     stats_ecozone <- plot_raster_stats(
#'       raster = age,
#'       polygons = ecozones,
#'       polygon_id = id_col,
#'       filter_ids = NULL,
#'       counts_df = ecozone_counts,
#'       stats_df = ecozone_stats,
#'       aggregate_factor = 1,
#'       remove_zeros = FALSE,
#'       raster_label = "age",
#'       inset_canada = TRUE,
#'       bin_width = 10,
#'       output_dir = tmp_pth,
#'       csv_file = NULL
#'     )
#'
#'     gc()
#'
#'     ## review the plots produced in tmp_pth
#'     list.files(tmp_pth, pattern = "[.]png", full.names = TRUE)
#'
#'     identical(ecozone_stats, stats_ecozone) ## TRUE
#'   }
#' }
#'
#' @export
#' @rdname raster_stats
calc_raster_counts <- function(raster, polygons = NULL, polygon_id = NULL, filter_ids = NULL) {
  stopifnot(
    requireNamespace("dplyr", quietly = TRUE),
    inherits(raster, "SpatRaster"),
    inherits(polygons, c("sf", "SpatVector"))
  )

  polygons <- prep_polygons(raster, polygons, polygon_id, filter_ids)
  poly_names_df <- data.frame(zone = seq_along(polygons$ID), ID = polygons$ID)

  ## terra::freq outputs data.frame with colnames: layer, value, count, zone (zone is numeric)
  ## we need data.frame with colnames: ID, value, count (with ID being the poly name)
  terra::freq(raster, zones = polygons) |>
    dplyr::left_join(poly_names_df, by = "zone") |>
    dplyr::relocate(ID, .before = layer) |>
    dplyr::mutate(layer = NULL, zone = NULL) |>
    dplyr::arrange(ID, value)
}

## Weighted quantile from a value-frequency table, using the inverse-ECDF
## definition (smallest value whose cumulative proportion reaches `p`). This is
## well-behaved for integer rasters with large point masses (e.g. many zeros),
## where an interpolated definition would report a non-zero median even when the
## majority of cells are zero.
weighted_quantile <- function(value, count, probs) {
  o <- order(value)
  value <- value[o]
  ## as.numeric() before accumulating: `count` comes from a pixel tally, and
  ## cumsum() on an integer vector overflows silently past 2^31 - 1. A single
  ## ecozone of SCANFI at 30 m is already ~1.5e9 pixels.
  count <- as.numeric(count[o])
  cprop <- cumsum(count) / sum(count)
  vapply(probs, function(p) value[which(cprop >= p)[1]], numeric(1))
}

## Derive per-polygon summary statistics from a value-frequency table (as
## returned by `calc_raster_counts()`). For an integer raster every statistic is
## an exact function of that table, so only a few hundred rows per zone are held
## in memory rather than every pixel.
stats_from_counts <- function(counts_df) {
  stopifnot(all(c("ID", "value", "count") %in% names(counts_df)))

  counts_df <- counts_df[!is.na(counts_df$value), , drop = FALSE]

  split(counts_df, counts_df$ID) |>
    lapply(function(d) {
      ## Work in doubles. `value` and `count` are both integer as read, and at
      ## national 30 m scale the products and totals blow past the 2^31 - 1
      ## integer limit: one ecozone holds ~1.5e9 pixels, so `value * count`
      ## overflows to NA (silently, apart from a warning) and even `sum(count)`
      ## is within a factor of two of overflowing.
      value <- as.numeric(d$value)
      count <- as.numeric(d$count)
      tot <- sum(count)
      q <- weighted_quantile(value, count, c(0.25, 0.50, 0.75))
      data.frame(
        ID = d$ID[1],
        min = min(d$value),
        mean = sum(value * count) / tot,
        max = max(d$value),
        q25 = q[1],
        q50 = q[2], ## median
        q75 = q[3],
        prop_zero = sum(count[value == 0]) / tot
      )
    }) |>
    do.call(what = rbind) |>
    `rownames<-`(NULL)
}

#' @export
#' @rdname raster_stats
calc_raster_stats <- function(raster, polygons = NULL, polygon_id = NULL,
                              filter_ids = NULL, counts_df = NULL) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))

  ## Statistics are derived from the on-disk value-frequency table rather than by
  ## reading every pixel, which keeps memory use bounded on very large rasters.
  ## (For a *continuous* / float raster the frequency table is no longer small;
  ## summarize out-of-core instead -- e.g. with exactextractr/zonal, or by
  ## extracting to a parquet dataset partitioned by zone and using `arrow`.)
  ## Pass a precomputed `counts_df` (from calc_raster_counts()) to avoid
  ## recomputing the frequency table.
  if (is.null(counts_df)) {
    stopifnot(!is.null(polygons), !is.null(polygon_id))
    counts_df <- calc_raster_counts(raster, polygons, polygon_id, filter_ids)
  }

  stats_from_counts(counts_df)
}

## ---------------------------------------------------------------------------
## Panel builders for `plot_raster_stats()`
##
## These are factored out of the per-polygon loop so that the vignette artifact
## manifest (see `.artifact_hash()`) can pin the committed figures to exactly the
## code that draws them, and so the loop body stays readable.
## ---------------------------------------------------------------------------

.hist_panel <- function(poly_counts, bin_width, raster_label, region_val) {
  poly_counts_binned <- poly_counts |>
    dplyr::mutate(bin = floor(value / bin_width) * bin_width) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(count = sum(count), .groups = "drop")

  ggplot2::ggplot(poly_counts_binned, ggplot2::aes(x = bin, y = count)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Histogram:", region_val), x = raster_label, y = "Count")
}

## raster + polygon outline -- use tidyterra for SpatRaster/SpatVector (mixing
## terra objects with ggplot2::geom_sf() does not work). geom_spatraster()
## auto-maps the (single) layer's values to `fill`; a continuous scale is
## required for the continuous raster (scale_fill_viridis_d() collapsed every
## cell to one colour -> solid rectangle).
.map_panel <- function(r_mask, poly, raster_label, region_val) {
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = r_mask) +
    tidyterra::geom_spatvector(data = poly, color = "black", fill = NA) +
    ggplot2::scale_fill_viridis_c(name = raster_label, na.value = "transparent") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Map:", region_val), x = "Longitude", y = "Latitude")
}

.inset_panel <- function(canada, poly) {
  ggplot2::ggplot() +
    tidyterra::geom_spatvector(data = canada, fill = "grey85", color = NA) +
    tidyterra::geom_spatvector(data = poly, fill = "darkred", color = NA) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
}

.stats_panel <- function(poly_stats, region_val) {
  ## Label and value are drawn as two separate right/left-aligned columns about
  ## a shared centre. Pasting them into one centred string (as this used to do)
  ## runs them together as soon as a value is wide -- "Min : 0Mean : 60.6Max".
  labels <- c("Min", "Mean", "Max", "25%", "Median", "75%", "% Zero")
  values <- c(
    round(poly_stats$min, 2),
    round(poly_stats$mean, 2),
    round(poly_stats$max, 2),
    round(poly_stats$q25, 2),
    round(poly_stats$q50, 2),
    round(poly_stats$q75, 2),
    paste0(round(poly_stats$prop_zero * 100, 2), "%")
  )

  ## three columns, three rows; the last row holds a single centred entry
  centres <- c(1.9, 5.0, 8.1)
  x <- c(centres, centres, 5.0)
  y <- c(6.5, 6.5, 6.5, 5.0, 5.0, 5.0, 3.5)
  gap <- 0.25

  ggplot2::ggplot(data.frame(x = 0:10, y = 0:10)) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y), alpha = 0.0) +
    ggplot2::annotate(
      "text",
      x = x - gap,
      y = y,
      label = paste0(labels, ":"),
      hjust = 1,
      vjust = 0.5
    ) +
    ggplot2::annotate(
      "text",
      x = x + gap,
      y = y,
      label = as.character(values),
      hjust = 0,
      vjust = 0.5
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste0("Statistics: ", region_val))
}

## Assemble and write the combined figure for a single polygon. Everything that
## determines the contents of the saved .png lives here (panels, layout, device
## size/resolution, filename), so hashing this function plus the panel builders
## is sufficient to detect that committed figures have gone stale.
.write_polygon_figure <- function(poly, r_mask, poly_counts, poly_stats, region_val,
                                  raster_label, bin_width, canada, output_dir,
                                  fig_width = 12, fig_height = 7.5, fig_dpi = 300) {
  hist_plot <- .hist_panel(poly_counts, bin_width, raster_label, region_val)
  map_plot <- .map_panel(r_mask, poly, raster_label, region_val)
  stats_plot <- .stats_panel(poly_stats, region_val)

  if (!is.null(canada)) {
    map_plot <- map_plot +
      patchwork::inset_element(
        .inset_panel(canada, poly),
        left = 0.75,
        right = 0.95,
        bottom = 0.75,
        top = 0.95,
        align_to = "plot"
      )
  }

  final_plot <- hist_plot |
    (map_plot / stats_plot + patchwork::plot_layout(heights = c(3, 1))) +
      patchwork::plot_annotation(title = region_val)

  fig_file <- file.path(
    output_dir,
    paste0("region_", gsub("[^A-Za-z0-9]", "_", region_val), ".png")
  )

  ggplot2::ggsave(
    filename = fig_file,
    plot = final_plot,
    width = fig_width,
    height = fig_height,
    dpi = fig_dpi
  )

  fig_file
}

#' @param csv_file Optional. Base filename to save the output CSV file.
#' If provided, a suffix will be added to this filename to denote the 'counts' and 'stats' tables
#' (e.g., if `csv_file = "ecozones.csv"`, the output files saved to `output_dir` will be
#' `ecozones_counts.csv` and `ecozones_stats.csv`).
#'
#' @param aggregate_factor Optional. Integer factor to aggregate raster for plotting. Default is `1` (no aggregation).
#'
#' @param remove_zeros Logical. Removes zeroes from raster and counts for plotting.
#' Useful to eliminate non-forested areas from forest rasters.
#'
#' @param counts_df optional `data.frame` of frequency counts output by `calc_raster_counts`.
#' If not provided, that function will be called internally.
#'
#' @param stats_df optional `data.frame` of raster statistics output by `calc_raster_stats`.
#' If not provided, that function will be called internally.
#'
#' @param raster_label Optional. Character label for raster used in plots.
#' Defaults to raster layer name.
#'
#' @param inset_canada Logical. Generates a small inset map of plotted polygons within Canada.
#'
#' @param inset_polygons Optional `sf` or `SpatVector` outline to draw as the
#' backdrop of the inset map. When `NULL` (the default) the Canada boundary is
#' downloaded via [gadm_canada()]. Supply it to keep the call offline -- e.g.
#' `sf::st_union()` of the polygons being summarized. Ignored when
#' `inset_canada = FALSE`.
#'
#' @param bin_width Numeric. Bin size for histogram plots.
#'
#' @param output_dir Directory to save figures and optional csv file.
#'
#' @param fig_width,fig_height Numeric. Size of each saved figure, in inches.
#'
#' @param fig_dpi Numeric. Resolution of each saved figure. The defaults suit
#' print; drop `fig_dpi` (and the dimensions) well below them when the figures
#' are destined for a web page or a vignette, where a 300 dpi 12 x 7.5 in PNG is
#' about 20x larger than it needs to be.
#'
#' @return A `data.frame` of computed statistics (via `calc_raster_stats()`), with side effects
#' of saving summary plots (and optionally a \file{.csv}) to disk.
#'
#' @export
#' @rdname raster_stats
plot_raster_stats <- function(
  raster,
  polygons = NULL,
  polygon_id = NULL,
  filter_ids = NULL,
  counts_df = NULL,
  stats_df = NULL,
  aggregate_factor = 1,
  remove_zeros = FALSE,
  raster_label = NULL,
  inset_canada = TRUE,
  inset_polygons = NULL,
  bin_width = 10,
  output_dir = ".",
  csv_file = NULL,
  fig_width = 12,
  fig_height = 7.5,
  fig_dpi = 300
) {
  stopifnot(
    requireNamespace("dplyr", quietly = TRUE),
    inherits(raster, "SpatRaster"),
    inherits(polygons, c("sf", "SpatVector")),
    is.logical(inset_canada) && !is.na(inset_canada)
  )

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (is.null(raster_label)) {
    raster_label <- terra::names(raster)[1]
  }

  ## Aggregate raster for plotting
  raster_plot <- if (aggregate_factor > 1) {
    terra::aggregate(raster, fact = aggregate_factor)
  } else {
    terra::deepcopy(raster)
  }

  ## Remove zeros if requested
  if (remove_zeros) {
    raster_plot <- terra::clamp(raster_plot, lower = 1, value = FALSE)
    raster <- terra::clamp(raster, lower = 1, value = FALSE) # also clamp original raster for consistency
  }

  ## Calculate raster stats only when not supplied. These can be very expensive
  ## on large rasters, so callers may pass precomputed `counts_df`/`stats_df`
  ## (e.g. from calc_raster_counts()/calc_raster_stats()) to avoid recomputation.
  if (is.null(counts_df)) {
    counts_df <- calc_raster_counts(raster, polygons, polygon_id, filter_ids)
  }
  if (is.null(stats_df)) {
    stats_df <- calc_raster_stats(raster, polygons, polygon_id, filter_ids)
  }

  polygons <- prep_polygons(raster, polygons, polygon_id, filter_ids) ## do after calculating stats

  if (!is.null(csv_file)) {
    utils::write.csv(
      counts_df,
      file.path(output_dir, .suffix(csv_file, "_counts")),
      row.names = FALSE
    )
    utils::write.csv(
      stats_df,
      file.path(output_dir, .suffix(csv_file, "_stats")),
      row.names = FALSE
    )
  }

  ## Outline for the inset map: resolve ONCE (not once per polygon, which is slow
  ## and fragile). A caller that already has a suitable outline can pass it as
  ## `inset_polygons`, which avoids the download entirely -- useful offline, and
  ## in vignettes and tests. Otherwise fetch it, and degrade gracefully to no
  ## inset if the data source is temporarily unavailable rather than aborting the
  ## whole run. A NULL `canada` is the signal to `.write_polygon_figure()` to
  ## skip the inset.
  canada <- NULL
  if (inset_canada) {
    canada <- if (!is.null(inset_polygons)) {
      terra::project(terra::vect(inset_polygons), raster)
    } else {
      tryCatch(
        gadm_canada(src = "geodata", dst_path = tempdir()) |> terra::project(raster),
        error = function(e) NULL
      )
    }
    if (is.null(canada)) {
      warning(
        "Could not retrieve the Canada boundary for the inset map; ",
        "plotting without inset.",
        call. = FALSE
      )
    }
  }

  ## Build plots for each polygon
  for (i in seq_len(nrow(polygons))) {
    region_val <- polygons[["ID"]][i, ]

    message("Plotting polygon: ", region_val)

    poly <- polygons[i, ]

    ## Skip if polygon doesn't intersect raster extent
    check_intersect <- suppressWarningsSpecific(
      terra::relate(poly, raster, relation = "intersects"),
      "partial argument match of 'ext' to 'extent'"
    )
    if (!check_intersect) {
      message("Skipping ", region_val, " (no spatial intersection)")
      next
    }

    poly_counts <- subset(counts_df, ID == region_val)
    poly_stats <- subset(stats_df, ID == region_val)

    ## Crop aggregated raster for plotting
    r_mask <- terra::mask(terra::crop(raster_plot, poly), poly)

    ## Skip if raster is empty (all NA)
    if (all(is.na(terra::values(r_mask)))) {
      message("Skipping ", region_val, " (raster empty after crop)")
      next
    }

    .write_polygon_figure(
      poly = poly,
      r_mask = r_mask,
      poly_counts = poly_counts,
      poly_stats = poly_stats,
      region_val = region_val,
      raster_label = raster_label,
      bin_width = bin_width,
      canada = canada,
      output_dir = output_dir,
      fig_width = fig_width,
      fig_height = fig_height,
      fig_dpi = fig_dpi
    )

    gc()
  }

  return(stats_df)
}
