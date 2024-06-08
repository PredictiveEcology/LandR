test_that("leading species transitions plots look good", {
  skip("needs at least 20GB RAM")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("dplyr")
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("map")
  skip_if_not_installed("withr")

  withr::local_package("data.table")
  withr::local_package("ggplot2")
  withr::local_package("ggalluvial")
  withr::local_package("terra")

  run <- 1L
  outputDir <- file.path("~/GitHub/BC_HRV/outputs",
                         "NRD_Quesnel_scfm_hrv_FRT_res125",
                         sprintf("rep%02d", run))

  ml <- readRDS(file.path(outputDir, "ml_preamble.rds"))
  rTM <- terra::rast(file.path(outputDir, "pixelGroupMap_year0000.tif")) |> terra::rast()
  studyArea2 <- map::studyArea(ml, 2) ## studyAreaReporting
  NDTBEC <- suppressWarnings({
    sf::st_crop(ml$`ecoregionLayer (NDTxBEC)`, studyArea2)
  })
  rstNDTBEC <- terra::rasterize(NDTBEC, rTM, field = "NDTBEC") |>
    terra::crop(studyArea2) |>
    terra::mask(studyArea2)
  rm(ml)

  years <- c(800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200)

  fvtm <- vapply(outputDir, function(d) {
    file.path(d, sprintf("vegTypeMap_year%04d.tif", years))
  }, character(length(years))) |> as.vector()

  ## using VTM as-is ---------------------------------------------------------------------------------

  transitions_df <- vegTransitions(
    fvtm = fvtm,
    ecoregion = NDTBEC,
    field = "NDTBEC",
    studyArea = studyArea2,
    times = years
  )

  if (interactive()) {
    transition_ggs <- plotVegTransitions(transitions_df)

    lapply(seq_along(transition_ggs), function(i) {
      ggsave(file.path(outputDir, "figures", paste0("transition_vegTypeMap_", i, ".png")),
             transition_ggs[[i]], width = 12, height = 6)
    })
  }

  ## using VTM to get conifer/deciduous/mixed --------------------------------------------------------

  # tmp <- SpaDES.core::loadSimList(file.path(outputDir, "simOutSpeciesLayers_NRD_Quesnel.rds"))
  # fwrite(tmp$sppEquiv, file.path(outputDir, "sppEquiv.csv"))
  # rm(tmp)

  sppEquiv <- fread(file.path(outputDir, "sppEquiv.csv"))
  sppEquiv <- sppEquiv[, c("BC_HRV", "Type")] |>
    rbind(data.table(BC_HRV = "Mixed", Type = "Mixed"))

  transitions_df2 <- lapply(seq_along(years), function(yr) {
    vtm <- terra::rast(fvtm[yr]) |>
      terra::crop(studyArea2) |>
      terra::mask(studyArea2)
    lvls_vt <- levels(vtm)[[1]]
    names(lvls_vt) <- tolower(names(lvls_vt))

    vegType <- lvls_vt[["values"]][match(values(vtm, mat = FALSE), lvls_vt[["id"]])]
    conifdecid <- sppEquiv$Type[match(vegType, sppEquiv$BC_HRV)]

    tdf <- data.table(
      pixelID = seq_len(terra::ncell(vtm)),
      ecoregion = lvls_er[["ndtbec"]][match(values(rstNDTBEC, mat = FALSE), lvls_er[["id"]])],
      vegType = conifdecid,
      year = years[yr]
    ) |>
      na.omit("ecoregion")

    tdf <- tdf[is.na(vegType), vegType := "none"]

    return(as.data.frame(tdf))
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(year = factor(year, levels = as.character(years)))

  if (interactive()) {
    transition_ggs2 <- plotVegTransitions(transitions_df2)

    lapply(seq_along(transition_ggs2), function(i) {
      ggsave(file.path(outputDir, "figures", paste0("transition_conifdecid_", er, ".png")), gg,
             width = 12, height = 6)
    })
  }

  withr::deferred_run()
})
