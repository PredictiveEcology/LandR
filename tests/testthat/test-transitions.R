test_that("leading species transitions plots look good", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("dplyr")
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("map")
  skip_if_not_installed("memuse")
  skip_if_not_installed("withr")

  withr::local_package("data.table")
  withr::local_package("ggplot2")
  withr::local_package("ggalluvial")
  withr::local_package("memuse")
  withr::local_package("terra")

  ## need ~20GB RAM to construct summary data.frames and plots
  skipifnot(isTRUE(Sys.meminfo()$freeram >= as.memuse(20 * 1024^3)))

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

    rm(transition_ggs)
  }

  rm(transition_df)

  ## using VTM to get conifer/deciduous/mixed --------------------------------------------------------

  # tmp <- SpaDES.core::loadSimList(file.path(outputDir, "simOutSpeciesLayers_NRD_Quesnel.rds"))
  # fwrite(tmp$sppEquiv, file.path(dirname(outputDir), "sppEquiv.csv"))
  # rm(tmp)

  sppEquiv <- fread(file.path(dirname(outputDir), "sppEquiv.csv"))
  sppEquiv <- sppEquiv[, c("BC_HRV", "Type")] |>
    rbind(data.table(BC_HRV = "Mixed", Type = "Mixed"))

  fvtm2 <- vtm2conifdecid(
    fvtm = fvtm,
    sppEquiv = sppEquiv,
    sppEquivCol = "BC_HRV",
    studyArea = studyArea2
  )

  transitions_df <- vegTransitions(
    vtm = fvtm2,
    ecoregion = NDTBEC,
    field = "NDTBEC",
    studyArea = studyArea2,
    times = years
  )

  if (interactive()) {
    transition_ggs2 <- plotVegTransitions(transitions_df2)

    lapply(seq_along(transition_ggs2), function(i) {
      ggsave(file.path(outputDir, "figures", paste0("transition_conifdecid_", er, ".png")), gg,
             width = 12, height = 6)
    })

    rm(transition_ggs2)
  }

  rm(transition_df)

  withr::deferred_run()
})
