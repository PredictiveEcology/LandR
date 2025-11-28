utils::globalVariables(c(".data", "bin", "count", "geometry", "ID", "id_col"))

prep_polygons <- function(raster, polygons, polygon_id = NULL, filter_ids = NULL) {
  polygons <- sf::st_transform(polygons, sf::st_crs(raster)) |>
    dplyr::rename(ID = .data[[polygon_id]]) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")

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
#' @note Large rasters should be processed on disk, where possible, e.g. by setting
#' `terraOptions(memfrac = 0.0)`.
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
#'       requireNamespace("purrr", quietly = TRUE) &&
#'       requireNamespace("withr", quietly = TRUE) &&
#'       requireNamespace("zonal", quietly = TRUE)) {
#'     tmp_pth <- withr::local_tempdir()
#'
#'     ## work with raster on disk instead of in memory
#'     terraOptions(memfrac = 0.0)
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
#'     ecozone_stats <- calc_raster_stats(
#'      raster = age,
#'      polygons = ecozones,
#'       polygon_id = id_col
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
calc_raster_counts <- function(
    raster,
    polygons = NULL,
    polygon_id = NULL,
    filter_ids = NULL
) {
  stopifnot(
    inherits(raster, "SpatRaster"),
    inherits(polygons, c("sf", "SpatVector"))
  )

  polygons <- prep_polygons(raster, polygons, polygon_id, filter_ids)
  poly_names_df <- data.frame(zone = seq_along(polygons$ID), ID = polygons$ID)

  ## terra::freq outputs data.frame with colnames: layer, value, count, zone (zone is numeric)
  ## we need data.frame with colnames: ID, value, count (with ID being the poly name)
  out <- terra::freq(raster, zones = polygons) |>
    dplyr::left_join(poly_names_df, by = "zone") |>
    dplyr::relocate(ID, .before = layer) |>
    dplyr::mutate(layer = NULL, zone = NULL) |>
    dplyr::arrange(ID, value)

  return(out_df)
}

prop_zero <- function(df, ...) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))

  na.omit(df) |> dplyr::summarise(prop_zero = sum(value == 0) / length(value))
}

#' @export
#' @rdname raster_stats
calc_raster_stats <- function(raster, polygons = NULL, polygon_id = NULL, filter_ids = NULL) {
  stopifnot(
    requireNamespace("dplyr", quietly = TRUE),
    requireNamespace("purrr", quietly = TRUE),
    requireNamespace("zonal", quietly = TRUE),
    !is.null(polygons),
    !is.null(polygon_id)
  )

  polygons <- prep_polygons(raster, polygons, polygon_id, filter_ids) |> sf::st_as_sf()

  ## functions known to exactextractr::extract_extract()
  ee_funs <- list("min", "mean", "max", "quantile")

  ## custom funs for use with
  my_funs <- list("prop_zero")

  all_funs <- append(ee_funs, my_funs)
  ll <- lapply(all_funs, function(fun) {
    gc()

    if (fun == "quantile") {
      zonal::execute_zonal(
        data = raster,
        geom = polygons,
        ID = "ID",
        fun = fun,
        join = FALSE,
        quantiles = c(0.25, 0.50, 0.75)
      ) |>
        setNames(c("ID", "q25", "q50", "q75")) ## q50 is the median
    } else if (fun == "prop_zero") {
      zonal::execute_zonal(
        data = raster,
        geom = polygons,
        ID = "ID",
        fun = prop_zero, ## pass function, not character
        join = FALSE,
        summarize_df = TRUE
      ) |>
        setNames(c("ID", "prop_zero"))
    } else {
      zonal::execute_zonal(data = raster, geom = polygons, ID = "ID", fun = fun, join = FALSE) |>
        setNames(c("ID", fun))
    }
  }) |>
    purrr::reduce(dplyr::inner_join, by = "ID")
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
#' @param bin_width Numeric. Bin size for histogram plots.
#'
#' @param output_dir Directory to save figures and optional csv file.
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
    bin_width = 10,
    output_dir = ".",
    csv_file = NULL
) {
  stopifnot(
    requireNamespace("dplyr", quietly = TRUE),
    requireNamespace("purrr", quietly = TRUE),
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

  ## Calculated raster stats
  counts_df <- calc_raster_counts(raster, polygons, polygon_id, filter_ids)
  stats_df <- calc_raster_stats(raster, polygons, polygon_id, filter_ids)

  polygons <- prep_polygons(raster, polygons, polygon_id, filter_ids) ## do after calculating stats

  if (!is.null(csv_file)) {
    utils::write.csv(counts_df, file.path(output_dir, .suffix(csv_file, "_counts")), row.names = FALSE)
    utils::write.csv(stats_df, file.path(output_dir, .suffix(csv_file, "_stats")), row.names = FALSE)
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

    value_col <- terra::names(r_mask)

    ## Histogram from counts_df
    poly_counts_binned <- poly_counts |>
      dplyr::mutate(bin = floor(value / bin_width) * bin_width) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(count = sum(count), .groups = "drop")

    hist_plot <- ggplot2::ggplot(poly_counts_binned, ggplot2::aes(x = bin, y = count)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Histogram:", region_val), x = raster_label, y = "Count")

    ## Main map
    map_plot_base <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = r_mask, aes(fill = .data[[value_col]])) +
      ggplot2::geom_sf(data = poly, color = "black", fill = NA) +
      ggplot2::scale_fill_viridis_c(name = raster_label, na.value = "transparent") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Map:", region_val), x = "Longitude", y = "Latitude")

    if (inset_canada) {
      ## Canada outline for inset
      canada <- gadm_canada(src = "geodata", dst_path = tempdir()) |> terra::project(raster)
      canada_sf <- sf::st_as_sf(canada)

      ## Inset map
      inset_plot <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = canada, fill = "grey85", color = NA) +
        ggplot2::geom_sf(data = poly, fill = "darkred", color = NA) +
        ggplot2::theme_void() +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
        )

      map_with_inset <- map_plot_base +
        patchwork::inset_element(
          inset_plot,
          left = 0.75,
          right = 0.95,
          bottom = 0.75,
          top = 0.95,
          align_to = "plot"
        )
    }

    ## Stats plot (centered text)
    stats_text <- c(
      paste0("Min    : ", round(poly_stats$min, 2)),
      paste0("Mean   : ", round(poly_stats$mean, 2)),
      paste0("Max    : ", round(poly_stats$max, 2)),
      paste0("25%    : ", round(poly_stats$q25, 2)),
      paste0("Median : ", round(poly_stats$q50, 2)),
      paste0("75%    : ", round(poly_stats$q75, 2)),
      paste0("% Zero : ", round(poly_stats$prop_zero * 100, 2), "%")
    )
    stats_plot <- ggplot2::ggplot(data.frame(x = 0:10, y = 0:10)) +
      ggplot2::geom_point(aes(x = x, y = y), alpha = 0.0) +
      ggplot2::annotate(
        "text",
        x = c(2.5, 5.0, 7.5, 2.5, 5.0, 7.5, 5.0),
        y = c(6.5, 6.5, 6.5, 5.0, 5.0, 5.0, 3.5),
        label = stats_text,
        hjust = 0.5,
        vjust = 0.5
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0("Statistics: ", region_val))

    if (inset_canada) {
      final_plot <- hist_plot |
        (map_with_inset / stats_plot + patchwork::plot_layout(heights = c(3, 1))) +
        patchwork::plot_annotation(title = region_val)
    } else {
      final_plot <- hist_plot |
        (map_plot_base / stats_plot + patchwork::plot_layout(heights = c(3, 1))) +
        patchwork::plot_annotation(title = region_val)
    }
    fig_name <- paste0("region_", gsub("[^A-Za-z0-9]", "_", region_val), ".png")

    ggplot2::ggsave(
      filename = file.path(output_dir, fig_name),
      plot = final_plot,
      width = 12,
      height = 7.5,
      dpi = 300
    )

    gc()
  }

  return(stats_df)
}
