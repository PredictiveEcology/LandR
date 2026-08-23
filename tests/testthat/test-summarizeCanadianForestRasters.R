test_that("stats_from_counts derives statistics from a value-frequency table", {
  ## Hand-computed from a table with a large point mass at zero: the
  ## inverse-ECDF definition must keep the median at 0 when most cells are 0
  ## (an interpolated definition would report a non-zero median here).
  counts <- data.frame(
    ID = "A",
    value = c(0L, 10L, 20L),
    count = c(5L, 3L, 2L)
  )

  res <- LandR:::stats_from_counts(counts)

  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 1L)
  expect_identical(res$min, 0L)
  expect_identical(res$max, 20L)
  expect_equal(res$mean, 7) ## (0*5 + 10*3 + 20*2) / 10
  expect_equal(res$q25, 0) ## cumulative proportion reaches 0.25 within the zeros
  expect_equal(res$q50, 0)
  expect_equal(res$q75, 10) ## 0.8 >= 0.75, at value 10
  expect_equal(res$prop_zero, 0.5)
})

test_that("stats_from_counts drops NA values and handles multiple polygons", {
  counts <- data.frame(
    ID = c("A", "A", "A", "B", "B"),
    value = c(0L, 10L, NA_integer_, 4L, 6L),
    count = c(5L, 5L, 99L, 1L, 1L)
  )

  res <- LandR:::stats_from_counts(counts)

  expect_identical(res$ID, c("A", "B"))
  ## the NA row must not contribute to the totals
  expect_equal(res$mean, c(5, 5))
  expect_equal(res$prop_zero, c(0.5, 0))
})

test_that("calc_raster_stats reuses a supplied counts table", {
  counts <- data.frame(
    ID = c("A", "A", "B", "B"),
    value = c(0L, 8L, 2L, 4L),
    count = c(3L, 1L, 2L, 2L)
  )

  ## no raster or polygons needed when counts are supplied
  expect_identical(
    calc_raster_stats(raster = NULL, counts_df = counts),
    LandR:::stats_from_counts(counts)
  )
})

test_that("prep_polygons accepts sf, GeoPackage-backed sf, and SpatVector", {
  testthat::skip_if_not_installed("dplyr")

  ras <- terra::rast(
    terra::ext(c(-105, -104, 55, 56)),
    crs = "EPSG:4326",
    resolution = c(.1, .1),
    vals = 1:100
  )
  shp <- readRDS(test_path("fixtures", "ecodistricts_88_115.rds"))

  ## a shapefile-derived sf names its geometry column "geometry" ...
  expect_identical(attr(shp, "sf_column"), "geometry")
  from_sf <- LandR:::prep_polygons(ras, shp, "ECOREGION")

  ## ... but a GeoPackage names it "geom", which must work just the same
  gpkg <- withr::local_tempfile(fileext = ".gpkg")
  sf::st_write(shp, gpkg, quiet = TRUE)
  round_tripped <- sf::st_read(gpkg, quiet = TRUE)
  expect_identical(attr(round_tripped, "sf_column"), "geom")
  from_gpkg <- LandR:::prep_polygons(ras, round_tripped, "ECOREGION")

  ## SpatVector input is documented, so it must work too
  from_vect <- LandR:::prep_polygons(ras, terra::vect(shp), "ECOREGION")

  for (res in list(from_sf, from_gpkg, from_vect)) {
    expect_s4_class(res, "SpatVector")
    expect_setequal(as.character(res$ID), c("88", "115"))
  }

  ## the dissolve must give the same geometry regardless of input flavour
  expect_equal(terra::geom(from_sf), terra::geom(from_gpkg))
  expect_equal(terra::geom(from_sf), terra::geom(from_vect))
})

test_that("plot_raster_stats produces figures and CSVs for each polygon", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("withr")

  td <- withr::local_tempdir("dest_")

  crs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  ras <- terra::rast(
    terra::ext(c(-105, -104, 55, 56)),
    crs = crs,
    resolution = c(.1, .1),
    vals = 1:100
  )

  ## Committed fixture rather than a live download: sis.agr.gc.ca throttles CI
  ## runners. Regenerate with data-raw/make-test-fixtures.R.
  ecodistricts <- readRDS(test_path("fixtures", "ecodistricts_88_115.rds"))

  ## inset_canada = FALSE keeps this test offline (the inset fetches GADM).
  ras_stats <- plot_raster_stats(
    raster = ras,
    polygons = ecodistricts,
    polygon_id = "ECOREGION",
    filter_ids = c("88", "115"), ## 88 intersects the raster; 115 does not
    inset_canada = FALSE,
    output_dir = td,
    csv_file = "ecoregion.csv"
  )

  expected_files <- c("region_88.png", "ecoregion_counts.csv", "ecoregion_stats.csv")
  expect_true(all(file.exists(file.path(td, expected_files))))

  ## 115 does not intersect the raster, so no figure is written for it
  expect_false(file.exists(file.path(td, "region_115.png")))

  ## only the intersecting polygon contributes cells to the frequency table
  expect_setequal(as.character(ras_stats$ID), "88")
  counts <- utils::read.csv(file.path(td, "ecoregion_counts.csv"))
  expect_equal(
    ras_stats[order(ras_stats$ID), ],
    LandR:::stats_from_counts(counts)[order(LandR:::stats_from_counts(counts)$ID), ],
    ignore_attr = TRUE
  )
})

test_that("the ecodistrict polygons are still downloadable", {
  ## Fidelity check on the upstream source backing the committed fixture.
  ## Skipped on CI, where sis.agr.gc.ca throttles the connection.
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("withr")

  td <- withr::local_tempdir("ecod_")

  ecodistricts <- reproducible::prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    destinationPath = td
  ) |>
    sf::st_make_valid()

  expect_s3_class(ecodistricts, "sf")
  expect_true("ECOREGION" %in% names(ecodistricts))
  expect_true(all(c(88L, 115L) %in% ecodistricts$ECOREGION))
})
