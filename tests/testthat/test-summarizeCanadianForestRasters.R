test_that("calc_raster_stats and plot_raster_stats work", {
  testthat::skip_on_cran()
  # testthat::skip_on_ci()

  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("geodata")
  testthat::skip_if_not_installed("purrr")
  testthat::skip_if_not_installed("withr")

  td <- withr::local_tempdir("dest_")

  crs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  ras <- terra::rast(
    terra::ext(c(-105, -104, 55, 56)),
    crs = crs,
    resolution = c(.1, .1),
    vals = 1:100
  )
  ecodistricts <- reproducible::prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    destinationPath = tempdir()
  ) |>
    sf::st_make_valid()

  ras_stats <- plot_raster_stats(
    raster = ras,
    polygons = ecodistricts,
    polygon_id = "ECOREGION",
    filter_ids = c("88", "147"), ## 88 intesects; 147 does not
    output_dir = td,
    csv_file = "ecoregion.csv"
  )

  expected_files <- c("region_88.png", "ecoregion_counts.csv", "ecoregion_stats.csv")
  expect_true(all(file.exists(file.path(td, expected_files))))
})
