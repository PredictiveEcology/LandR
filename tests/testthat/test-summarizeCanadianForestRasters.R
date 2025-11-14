test_that("SummarizeCanadianForestRasters works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  testthat::skip_if_not_installed("withr")

  td <- withr::local_tempdir("dest_")

  withr::local_package("terra")

  crs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  ras <- terra::rast(terra::ext(c(-105, -104, 55, 56)), crs = crs,
                      resolution = c(.1, .1), vals = 1:100)

  ## warn that rasterToMatchLarge and rasterToMatch are both missing
  testthat::expect_message(
    ras_stats <- plot_raster_statistics(
      ras,
      output_dir = td
    )
  )
})
