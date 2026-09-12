## a template raster standing in for `rasterToMatch`
rtmFixture <- function() {
  ras <- terra::rast(
    nrows = 6, ncols = 6, xmin = 0, xmax = 600, ymin = 0, ymax = 600, crs = "EPSG:3978"
  )
  terra::values(ras) <- 1
  ras
}

test_that("assertSpeciesLayers passes with zero species layers", {
  empty <- LandR:::.emptySpatRaster(rtmFixture())
  expect_equal(terra::nlyr(empty), 0L)
  expect_no_error(assertSpeciesLayers(empty, thresh = 10, doAssertion = TRUE))
})

test_that("assertSpeciesLayers still stops when all pixels are NA", {
  allNA <- rtmFixture()
  terra::values(allNA) <- NA_real_
  expect_error(
    assertSpeciesLayers(allNA, thresh = 10, doAssertion = TRUE),
    "no pixels found"
  )
})

test_that("the empty speciesLayers preserves the geometry of rasterToMatch", {
  rtm <- rtmFixture()
  empty <- LandR:::.emptySpatRaster(rtm)
  expect_equal(terra::nlyr(empty), 0L)
  expect_equal(names(empty), character(0))
  expect_true(LandR:::.compareRas(empty, rtm, stopOnError = FALSE))
})
