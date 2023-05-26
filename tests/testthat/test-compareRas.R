test_that("test .compareRas, .compareCRS", {
  require(reproducible)
  require(raster)
  require(terra)
  opts <- options("reproducible.inputPaths" = NULL,
                  "reproducible.overwrite" = TRUE,
                  "reproducible.useTerra" = TRUE,
                  "reproducible.rasterRead" = "terra::rast")

  on.exit({
    options(opts)
  }, add = TRUE)

  targetFile <- "rasterTest.tif"
  url <- "https://github.com/tati-micheletti/host/raw/master/data/rasterTest.tif"

  ras <- prepInputs(url = url,
                    targetFile = targetFile)
  ras2 <- project(ras, "EPSG:2169")
  expect_true(.compareRas(ras, ras))
  expect_true(.compareRas(ras, ras, ras))
  expect_false(.compareRas(ras, ras, ras2))

  ras3 <- terra::extend(ras, 10)
  expect_error(.compareRas(ras, ras3))
  expect_false(.compareRas(ras, ras3, stopOnError = FALSE))

  ## and with RasterLayer
  ras <- prepInputs(url = url,
                    fun = "raster::raster",
                    targetFile = targetFile)
  ras2 <- projectRaster(ras, crs = crs("EPSG:2169", proj = TRUE))
  expect_true(.compareRas(ras, ras))
  expect_true(.compareRas(ras, ras, ras))
  expect_false(.compareRas(ras, ras, ras2))

  ras3 <- raster::extend(ras, 10)
  expect_true(.compareRas(ras, ras3, ext = FALSE, rowcol = FALSE))
  expect_error(.compareRas(ras, ras3))
  expect_false(.compareRas(ras, ras3, stopOnError = FALSE))
})
