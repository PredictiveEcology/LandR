testthat::test_that("test .compareRas, .compareCRS -- rasters only", {
  testthat::skip_if_offline()
  testthat::skip_if_not_installed("withr")

  td <- withr::local_tempdir("dest_")

  withr::local_package("reproducible")

  withr::local_options(list(
    reproducible.inputPaths = NULL,
    reproducible.overwrite = TRUE,
    reproducible.useTerra = TRUE,
    reproducible.rasterRead = "terra::rast"
  ))

  url <- "https://github.com/tati-micheletti/host/raw/master/data/rasterTest.tif"
  targetFile <- basename(url)

  ras <- prepInputs(url = url, destinationPath = td, targetFile = targetFile)
  ras2 <- terra::project(ras, "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +datum=NAD83 +units=m +no_defs")
  testthat::expect_true(.compareRas(ras, ras))
  testthat::expect_true(.compareRas(ras, ras, ras))
  testthat::expect_error(.compareRas(ras, ras, ras2))
  testthat::expect_true(.compareRas(ras, ras, ras2, crs = FALSE, ext = FALSE))
  testthat::expect_false(.compareRas(ras, ras, ras2, stopOnError = FALSE))

  ras3 <- terra::extend(ras, 10)
  testthat::expect_error(.compareRas(ras, ras3))
  testthat::expect_false(.compareRas(ras, ras3, stopOnError = FALSE))

  ## and with RasterLayer
  ras <- prepInputs(
    url = url, destinationPath = td,
    fun = "raster::raster", targetFile = targetFile
  )
  ras2 <- raster::projectRaster(ras, crs = terra::crs("+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +datum=NAD83 +units=m +no_defs", proj = TRUE))
  testthat::expect_true(.compareRas(ras, ras))
  testthat::expect_true(.compareRas(ras, ras, ras))
  testthat::expect_error(.compareRas(ras, ras, ras2))
  testthat::expect_false(.compareRas(ras, ras, ras2, stopOnError = FALSE))

  ras3 <- raster::extend(ras, 10)
  testthat::expect_true(.compareRas(ras, ras3, ext = FALSE, rowcol = FALSE))
  testthat::expect_error(.compareRas(ras, ras3))
  testthat::expect_false(.compareRas(ras, ras3, stopOnError = FALSE))
})

testthat::test_that("test .compareRas, .compareCRS -- vectors only", {
  f <- system.file("ex/lux.shp", package = "terra")
  v <- terra::vect(f)
  v2 <- terra::project(v, terra::crs(v))

  testthat::expect_true(.compareRas(v, v))
  testthat::expect_true(.compareRas(v, v, v))
  testthat::expect_true(.compareRas(v, v, v2))

  v2 <- terra::project(v2, "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +datum=NAD83 +units=m +no_defs")
  testthat::expect_error(.compareRas(v, v2))
  testthat::expect_false(.compareRas(v, v2, stopOnError = FALSE))

  v3 <- terra::buffer(v, 10)
  testthat::expect_error(.compareRas(v, v3))
  testthat::expect_false(.compareRas(v, v3, stopOnError = FALSE))
  testthat::expect_true(.compareRas(v, v3, ext = FALSE))

  ## with SPDF
  v4 <- raster::shapefile(f)
  v5 <- raster::buffer(v4, 10)
  testthat::expect_true(.compareRas(v, v4))
  testthat::expect_true(.compareRas(v, v3, v4, ext = FALSE))
  testthat::expect_error(.compareRas(v, v3, v4))
  testthat::expect_error(.compareRas(v, v5, v4))
  testthat::expect_false(.compareRas(v, v5, v4, stopOnError = FALSE))

  testthat::expect_true(.compareRas(v3, v5))

  ## with sf
  v6 <- sf::st_read(f, quiet = TRUE)
  v7 <- sf::st_as_sf(v5)
  testthat::expect_true(.compareRas(v, v4, v6))
  testthat::expect_error(.compareRas(v6, v7))
  testthat::expect_false(.compareRas(v6, v7, stopOnError = FALSE))

  testthat::expect_error(.compareRas(v, v4, v6, v7))
  testthat::expect_false(.compareRas(v, v4, v6, v7, stopOnError = FALSE))
  testthat::expect_true(.compareRas(v, v4, v6, v7, ext = FALSE))

  testthat::expect_true(.compareRas(v3, v5, v7))
})

testthat::test_that("test .compareRas, .compareCRS -- vectors and rasters", {
  f <- system.file("ex/lux.shp", package = "terra")
  v <- terra::vect(f)

  ras <- terra::rast(v, resolution = 0.1)

  testthat::expect_true(.compareRas(v, ras, ext = FALSE))
  testthat::expect_error(.compareRas(v, ras))
  testthat::expect_false(.compareRas(v, ras, stopOnError = FALSE))

  v2 <- terra::project(v, "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +datum=NAD83 +units=m +no_defs")
  testthat::expect_error(.compareRas(ras, v2, ext = FALSE))
  testthat::expect_false(.compareRas(ras, v2, ext = FALSE, stopOnError = FALSE))
  testthat::expect_true(.compareRas(ras, v2, crs = FALSE, ext = FALSE))

  ## with raster and SPDF
  ras2 <- raster::raster(ras)
  v3 <- raster::shapefile(f)
  testthat::expect_true(.compareRas(v, v3, ras, ras2, ext = FALSE))
  testthat::expect_error(.compareRas(v, v3, ras, ras2))
  testthat::expect_false(.compareRas(v, v3, ras, ras2, stopOnError = FALSE))
})
