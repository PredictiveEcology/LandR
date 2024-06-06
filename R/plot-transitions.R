#' Create a vegetation transitions data frame
#'
#' @param fvtm character vector of filepaths to vegetation type maps.
#' @param ecoregion `SpatRaster` of ecoregion (or other) codes by which to group produce plots.
#' @param studyArea `sf` polygons object delineating the area to use for cropping and masking
#'                  of `ecoregion` (e.g., `studyAreaReporting`).
#' @param times numeric vector of years corresponding to the `fvtm` files.
#'
#' @return `data.frame` with columns `pixelID`, `ecoregion`, `vegType`, and `time`.
#'
#' @export
#' @rdname vegetation-transitions
vegTransitions <- function(fvtm, ecoregion, studyArea, times) {
  transitions_df <- lapply(seq_along(times), function(yr) {
    vtm <- terra::rast(fvtm[yr]) |>
      terra::crop(studyArea) |>
      terra::mask(studyArea)
    lvls_vt <- terra::levels(vtm)[[1]]

    tdf <- data.frame(
      pixelID = seq_len(terra::ncell(vtm)),
      ecoregion = lvls_er[["NDTBEC"]][match(values(ecoregion, mat = FALSE), lvls_er[["ID"]])],
      vegType = lvls_vt[["values"]][match(values(vtm, mat = FALSE), lvls_vt[["id"]])],
      time = times[yr]
    ) |> na.omit("vegType")

    return(tdf)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(time = factor(time, levels = as.character(times)))

  return(transitions_df)
}

#' Plot vegetation type transitions
#'
#' @note creating these plots for large landscapes can be computationally intensive
#' (time and memory use).
#'
#' @param transition_df A data frame with columns `pixelID`, `ecoregion`, `vegType`, and `time`.
#'                      (i.e., output of `vegTransitions()`).
#'
#' @return A list of ggplot objects, one for each ecoregion.
#' @export
#' @rdname vegetation-transitions
plotVegTransitions <- function(transitions_df) {
  stopifnot(
    requireNamespace("dplyr", quietly = TRUE),
    requireNamespace("ggalluvial", quietly = TRUE),
    requireNamespace("ggrepel", quietly = TRUE),

    all(c("pixelID", "ecoregion", "vegType", "time") %in% colnames(transitions_df))
  )

  erNames <- unique(transitions_df$ecoregion)
  transition_ggs <- lapply(erNames, function(er) {
    gg <- dplyr::filter(transitions_df, ecoregion == er) |>
      ggplot(aes(x = time, stratum = vegType, alluvium = pixelID, fill = vegType, label = vegType)) +
      scale_x_discrete(expand = c(0.1, 0)) +
      ggalluvial::geom_flow(color = "darkgray") +
      ggalluvial::geom_stratum(width = 1/8) +
      scale_linetype_manual(values = c("blank", "solid")) +
      ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(as.character(time)) == head(as.numeric(as.character(times)), 1), vegType, NA)),
        stat = "stratum", size = 4, direction = "y", nudge_x = -0.5
      ) +
      ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(as.character(time)) == tail(as.numeric(as.character(times)), 1), vegType, NA)),
        stat = "stratum", size = 4, direction = "y", nudge_x = +0.5
      ) +
      theme(legend.position = "none") +
      ggtitle(paste("Vegetation type transitions in", er))
  })
  names(transition_ggs) <- erNames

  transition_ggs
}
