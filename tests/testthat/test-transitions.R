# test_that("leading species transitions plots look good", {
#   skip_on_cran()
#   skip_on_ci()
#   skip_if_not_installed("dplyr")
#   skip_if_not_installed("ggalluvial")
#   skip_if_not_installed("ggrepel")
#   skip_if_not_installed("map")
#   skip_if_not_installed("withr")
#
#   withr::local_package("ggplot2")
#   withr::local_package("ggalluvial")
#
#   run <- 1L
#   outputDir <- file.path("~/GitHub/BC_HRV/outputs",
#                          "NRD_Quesnel_scfm_hrv_FRT_res125",
#                          sprintf("rep%02d", run))
#
#   skip_if_not(all(dir.exists(outputDir)))
#
#   ml <- readRDS(file.path(outputDir, "ml_preamble.rds"))
#   rTM <- terra::rast(file.path(outputDir, "pixelGroupMap_year0000.tif")) |> terra::rast()
#   studyArea3 <- map::studyArea(ml, 3)
#   NDTBEC <- suppressWarnings({
#     sf::st_crop(ml$`ecoregionLayer (NDTxBEC)`, studyArea3)
#   })
#   rstNDTBEC <- terra::rasterize(NDTBEC, rTM, field = "NDTBEC")
#   rm(ml)
#
#   lvls_er <- levels(rstNDTBEC)[[1]]
#   years <- c(800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200)
#
#   fvtm <- vapply(outputDir, function(d) {
#     file.path(d, sprintf("vegTypeMap_year%04d.tif", years))
#   }, character(length(years))) |> as.vector()
#
#
#   if (interactive()) {
#     transition_ggs <- plotVegTransitions()
#
#     lapply(erNames, function(er) {
#       ggsave(file.path(outputDir, "figures", paste0("transition_vegType_", er, ".png")),
#              transition_ggs[[er]], width = 12, height = 6)
#     })
#   }
#
#   withr::deferred_run()
# })
