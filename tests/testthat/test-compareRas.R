test_that("test .compareRas, .compareCRS -- rasters only", {
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
  expect_error(.compareRas(ras, ras, ras2))
  expect_true(.compareRas(ras, ras, ras2, crs = FALSE, ext = FALSE))
  expect_false(.compareRas(ras, ras, ras2, stopOnError = FALSE))

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
  expect_error(.compareRas(ras, ras, ras2))
  expect_false(.compareRas(ras, ras, ras2, stopOnError = FALSE))

  ras3 <- raster::extend(ras, 10)
  expect_true(.compareRas(ras, ras3, ext = FALSE, rowcol = FALSE))
  expect_error(.compareRas(ras, ras3))
  expect_false(.compareRas(ras, ras3, stopOnError = FALSE))
})

test_that("test .compareRas, .compareCRS -- vectors only", {
  require(reproducible)
  require(raster)
  require(terra)
  require(sf)
  opts <- options("reproducible.inputPaths" = NULL,
                  "reproducible.overwrite" = TRUE,
                  "reproducible.useTerra" = TRUE,
                  "reproducible.rasterRead" = "terra::rast")

  on.exit({
    options(opts)
  }, add = TRUE)

  f <- system.file("ex/lux.shp", package="terra")
  v <- vect(f)
  v2 <- project(v, crs(v))

  expect_true(.compareRas(v, v))
  expect_true(.compareRas(v, v, v))
  expect_true(.compareRas(v, v, v2))

  v2 <- project(v2, "EPSG:2169")
  expect_error(.compareRas(v, v2))
  expect_false(.compareRas(v, v2, stopOnError = FALSE))

  v3 <- terra::buffer(v, 10)
  expect_error(.compareRas(v, v3))
  expect_false(.compareRas(v, v3, stopOnError = FALSE))
  expect_true(.compareRas(v, v3, ext = FALSE))

  ## with SPDF
  v4 <- shapefile(f)
  v5 <- raster::buffer(v4, 10)
  expect_true(.compareRas(v, v4))
  expect_true(.compareRas(v, v3, v4, ext = FALSE))
  expect_error(.compareRas(v, v3, v4))
  expect_error(.compareRas(v, v5, v4))
  expect_false(.compareRas(v, v5, v4, stopOnError = FALSE))

  expect_true(.compareRas(v3, v5))

  ## with sf
  v6 <- st_read(f)
  v7 <- st_as_sf(v5)
  expect_true(.compareRas(v, v4, v6))
  expect_error(.compareRas(v6, v7))
  expect_false(.compareRas(v6, v7, stopOnError = FALSE))

  expect_error(.compareRas(v, v4, v6, v7))
  expect_false(.compareRas(v, v4, v6, v7, stopOnError = FALSE))
  expect_true(.compareRas(v, v4, v6, v7, ext = FALSE))

  expect_true(.compareRas(v3, v5, v7))
})

test_that("test .compareRas, .compareCRS -- vectors and rasters", {
  require(reproducible)
  require(raster)
  require(terra)
  require(sf)
  opts <- options("reproducible.inputPaths" = NULL,
                  "reproducible.overwrite" = TRUE,
                  "reproducible.useTerra" = TRUE,
                  "reproducible.rasterRead" = "terra::rast")

  on.exit({
    options(opts)
  }, add = TRUE)

  f <- system.file("ex/lux.shp", package="terra")
  v <- vect(f)

  ras <- rast(v, res = 0.1)

  expect_true(.compareRas(v, ras, ext = FALSE))
  expect_error(.compareRas(v, ras))
  expect_false(.compareRas(v, ras, stopOnError = FALSE))

  v2 <- project(v, "EPSG:2169")
  expect_error(.compareRas(ras, v2, ext = FALSE))
  expect_false(.compareRas(ras, v2, ext = FALSE, stopOnError = FALSE))
  expect_true(.compareRas(ras, v2, crs = FALSE, ext = FALSE))

  ## with raster and SPDF
  ras2 <- raster(ras)
  v3 <- shapefile(f)
  expect_true(.compareRas(v, v3, ras, ras2, ext = FALSE))
  expect_error(.compareRas(v, v3, ras, ras2))
  expect_false(.compareRas(v, v3, ras, ras2, stopOnError = FALSE))
})
