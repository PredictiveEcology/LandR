utils::globalVariables(c(
  "time", "vegType"
))

#' Create a summaries of vegetation type transitions
#'
#' @param vtm character vector of file paths to vegetation type maps (see [vegTypeMapGenerator()].
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param zones `sf` polygons object delineating the reporting area (`zone`) boundaries.
#'
#' @return
#' - `vtm2conifdecid()` returns a character vector of file paths to the conifer-deciduous maps.
#'
#' @export
#' @rdname vegetation-transitions
vtm2conifdecid <- function(vtm, sppEquiv = NULL, sppEquivCol = "LandR", zones) {
  if (is.null(sppEquiv)) {
    sppEquiv <- get(
      data("sppEquivalencies_CA", package = "LandR", envir = environment()),
      inherits = FALSE
    )
  }

  vapply(
    seq_along(vtm),
    function(i) {
      r <- terra::rast(vtm[i]) |> terra::crop(zones, mask = TRUE)
      lvls_vt <- levels(r)[[1]]
      names(lvls_vt) <- tolower(names(lvls_vt))

      vegType <- lvls_vt[["values"]][match(values(r, mat = FALSE), lvls_vt[["id"]])]

      conifdecid <- sppEquiv[["Type"]][match(vegType, sppEquiv[[sppEquivCol]])] |> as.factor()

      fout <- .suffix(vtm[i], "_conifdecid")
      rout <- terra::rast(r)
      terra::values(rout) <- conifdecid
      terra::writeRaster(rout, fout, overwrite = TRUE)

      return(fout)
    },
    character(1)
  )
}

#' @param field character string of the column name in `zones` to use for grouping.
#'
#' @param times numeric vector of years corresponding to `vtm`.
#'
#' @param na.rm logical. If `TRUE`, remove rows with `NA` values in `vegType`
#'        (they won't appear as a stratum in the alluvial diagram).
#'        If `FALSE`, these `NA` values will be replaced with `"_NA_"` so transitions
#'        between vegetated and non-vegetated pixels can be visualized.
#'
#' @param dest character, specifying a destination directory
#'
#' @return
#' - `vegTransitions()` returns a `arrow` dataset with columns `pixelID`, `zone`, `vegType`,
#' and `time`.
#'
#' @export
#' @rdname vegetation-transitions
vegTransitions <- function(vtm, zones, field, times, na.rm = FALSE, dest = ".") {
  stopifnot(requireNamespace("arrow", quietly = TRUE), requireNamespace("dplyr", quietly = TRUE))

  dest <- file.path(dest, "vegetation-transitions")

  rtm <- terra::rast(vtm[1]) |> terra::rast() ## remove values
  rstZones <- terra::rasterize(zones, rtm, field = field) |> terra::crop(zones, mask = TRUE)
  levels(rstZones) <- data.frame(ID = seq_len(nrow(zones))) |>
    dplyr::mutate({{ field }} := zones[[field]])

  purrr::walk(.x = seq_along(times), .f = function(yr) {
    r <- terra::rast(vtm[yr]) |> terra::crop(zones, mask = TRUE)
    lvls_vt <- terra::levels(r)[[1]]
    names(lvls_vt) <- tolower(names(lvls_vt))
    idcol_vt <- grep("^(id|value)$", names(lvls_vt), ignore.case = TRUE, value = TRUE)

    lvls_zn <- terra::levels(rstZones)[[1]]
    names(lvls_zn) <- tolower(names(lvls_zn))
    idcol_zn <- grep("^(id|value)$", names(lvls_zn), ignore.case = TRUE, value = TRUE)
    field <- tolower(field)

    tdf <- data.table(
      pixelID = seq_len(terra::ncell(r)),
      zone = lvls_zn[[field]][match(terra::values(rstZones, mat = FALSE), lvls_zn[[idcol_zn]])],
      vegType = lvls_vt[["values"]][match(terra::values(r, mat = FALSE), lvls_vt[[idcol_vt]])],
      time = times[yr]
    ) |>
      na.omit("zone")

    if (isTRUE(na.rm)) {
      tdf <- na.omit(tdf, "vegType")
    } else {
      tdf <- tdf[is.na(vegType), vegType := "_NA_"]
    }

    as.data.frame(tdf) |>
      dplyr::mutate(time = factor(time, levels = as.character(times))) |>
      dplyr::group_by(zone, time) |>
      arrow::write_dataset(dest, existing_data_behavior = "overwrite")
  })

  transitions_df <- arrow::open_dataset(dest)

  return(transitions_df)
}

#' Plot vegetation type transitions
#'
#' @note creating these plots for large landscapes can be computationally intensive
#' (time and memory use).
#'
#' @param transitions_df A data frame with columns `pixelID`, `zone`, `vegType`, and `time`.
#'                       (i.e., output of `vegTransitions()`).
#'
#' @return
#'  - `plotVegTransitions()` returns a list of `ggplot` objects, one for each zone.
#'
#' @export
#' @rdname vegetation-transitions
plotVegTransitions <- function(transitions_df) {
  stopifnot(
    requireNamespace("arrow", quietly = TRUE),
    requireNamespace("dplyr", quietly = TRUE),
    requireNamespace("ggalluvial", quietly = TRUE),
    requireNamespace("ggrepel", quietly = TRUE),
    all(c("pixelID", "zone", "vegType", "time") %in% colnames(transitions_df))
  )

  zoneNames <- transitions_df |> dplyr::distinct(zone) |> dplyr::collect() |> dplyr::pull(zone)
  times <- transitions_df |>
    dplyr::distinct(time) |>
    dplyr::collect() |>
    dplyr::pull(time) |>
    sort()

  transition_ggs <- lapply(zoneNames, function(er) {
    gg <- dplyr::filter(transitions_df, zone == er) |>
      dplyr::collect() |>
      dplyr::mutate(time = factor(time, levels = as.character(times))) |>
      ggplot(aes(
        x = time,
        stratum = vegType,
        alluvium = pixelID,
        fill = vegType,
        label = vegType
      )) +
      scale_x_discrete(expand = c(0.1, 0)) +
      ggalluvial::geom_flow(color = "darkgray") +
      ggalluvial::geom_stratum(width = 1 / 8) +
      scale_linetype_manual(values = c("blank", "solid")) +
      ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(as.character(time)) == min(times), vegType, NA)),
        stat = ggalluvial::StatStratum,
        size = 3,
        direction = "y",
        nudge_x = -0.5
      ) +
      ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(as.character(time)) == max(times), vegType, NA)),
        stat = ggalluvial::StatStratum,
        size = 3,
        direction = "y",
        nudge_x = +0.5
      ) +
      theme(legend.position = "none") +
      ggtitle(paste("Vegetation type transitions in", er))
  })
  names(transition_ggs) <- zoneNames

  transition_ggs
}
