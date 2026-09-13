## a template raster standing in for `rasterToMatch`
rtmFixture <- function() {
  ras <- terra::rast(
    nrows = 6, ncols = 6, xmin = 0, xmax = 600, ymin = 0, ymax = 600, crs = "EPSG:3978"
  )
  terra::values(ras) <- 1
  ras
}

## "No tree species" is carried as speciesLayers = NULL. A zero-layer SpatRaster was tried
## first and rejected: terra can build one only through its private constructor, and cannot
## wrap(), unwrap() or write it, so it did not survive Cache (2026-09-12).
test_that("assertSpeciesLayers passes with no species layers (NULL)", {
  expect_no_error(assertSpeciesLayers(NULL, thresh = 10, doAssertion = TRUE))
  expect_silent(assertSpeciesLayers(NULL, thresh = 10, doAssertion = TRUE))
})

test_that("assertSpeciesLayers still stops when all pixels are NA", {
  allNA <- rtmFixture()
  terra::values(allNA) <- NA_real_
  expect_error(
    assertSpeciesLayers(allNA, thresh = 10, doAssertion = TRUE),
    "no pixels found"
  )
})
