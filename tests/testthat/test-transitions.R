testthat::test_that("leading species transitions plots look good", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("ggalluvial")
  testthat::skip_if_not_installed("ggrepel")
  testthat::skip_if_not_installed("map")
  testthat::skip_if_not_installed("memuse")
  testthat::skip_if_not_installed("SpaDES.core")
  testthat::skip_if_not_installed("withr")

  withr::local_package("arrow")
  withr::local_package("data.table")
  withr::local_package("dplyr")
  withr::local_package("ggplot2")
  withr::local_package("ggalluvial")
  withr::local_package("memuse")
  withr::local_package("terra")

  ## need ~10GB RAM to construct summary data.frames and plots
  testthat::skip_if_not(isTRUE(Sys.meminfo()$freeram >= as.memuse(10 * 1024^3)))

  run <- 1L
  outputDir <- file.path(
    "~/GitHub/BC_HRV/outputs",
    "NRD_Quesnel_scfm_LH_hrv_NDTBEC_FRT_res125",
    sprintf("rep%02d", run)
  )

  testthat::skip_if_not(dir.exists(outputDir))

  ml <- readRDS(file.path(dirname(outputDir), "ml_preamble.rds"))
  rTM <- terra::rast(file.path(outputDir, "pixelGroupMap_year0000.tif")) |> terra::rast()
  studyArea2 <- map::studyArea(ml, 2) ## studyAreaReporting
  NDTBEC <- suppressWarnings(sf::st_crop(ml[["BEC zones"]], studyArea2)) |>
    dplyr::mutate(NDTBEC = paste0(NATURAL_DISTURBANCE, "_", ZONE))
  rm(ml)

  years <- seq(800, 1200, 50)

  fvtm <- file.path(outputDir, sprintf("vegTypeMap_year%04d.tif", years))

  stopifnot(all(file.exists(fvtm)))

  ## using VTM as-is ---------------------------------------------------------------------------------

  transitions_df <- vegTransitions(
    vtm = fvtm,
    zones = NDTBEC,
    field = "NDTBEC",
    times = years,
    na.rm = TRUE,
    dest = outputDir
  )

  if (interactive()) {
    transition_ggs <- plotVegTransitions(transitions_df)

    plot_files <- purrr::map_chr(.x = names(transition_ggs), .f = function(i) {
      ggsave(
        file.path(outputDir, "figures", paste0("transition_vegTypeMap_", i, ".png")),
        transition_ggs[[i]],
        width = 12,
        height = 6
      )
    })

    expect_all_true(file.exists(plot_files))

    rm(transition_ggs)
  }

  rm(transition_df)

  ## using VTM to get conifer/deciduous/mixed --------------------------------------------------------

  sppEquiv_file <- file.path(dirname(outputDir), "sppEquiv.csv")

  # tmp <- file.path(dirname(outputDir), "simOutDataPrep_NRD_Quesnel.rds") |>
  #   SpaDES.core::loadSimList()
  # fwrite(tmp$sppEquiv, sppEquiv_file)
  # rm(tmp)

  skip_if_not(file.exists(sppEquiv_file))

  sppEquiv <- fread(sppEquiv_file)
  sppEquiv <- sppEquiv[, c("BC_HRV", "Type")] |> rbind(data.table(BC_HRV = "Mixed", Type = "Mixed"))

  fvtm2 <- vtm2conifdecid(
    vtm = fvtm,
    sppEquiv = sppEquiv,
    sppEquivCol = "BC_HRV",
    studyArea = studyArea2
  )

  transitions_df <- vegTransitions(vtm = fvtm2, zones = NDTBEC, field = "NDTBEC", times = years)

  if (interactive()) {
    transition_ggs2 <- plotVegTransitions(transitions_df)

    lapply(names(transition_ggs2), function(i) {
      ggsave(
        file.path(outputDir, "figures", paste0("transition_conifdecid_", er, ".png")),
        gg,
        width = 12,
        height = 6
      )
    })

    rm(transition_ggs2)
  }

  rm(transition_df)

  withr::deferred_run()
})
