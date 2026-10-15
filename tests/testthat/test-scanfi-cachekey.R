## The Cache() of the per-species prepInputs() calls in loadSCANFISpeciesLayers() must be keyed
## on the spatial targets (to, cropTo, projectTo, maskTo), not only on the file names.
## Offline: the folder listing and prepInputs() are mocked.
test_that("loadSCANFISpeciesLayers does not reuse a cache entry across different cropTo", {
  skip_if_not_installed("terra")
  skip_if_not_installed("reproducible")

  cacheDir <- withr::local_tempdir()
  dPath <- withr::local_tempdir()
  withr::local_options(reproducible.cachePath = cacheDir, reproducible.useMemoise = FALSE)

  fakeListing <- data.table::data.table(
    name = c("SCANFI_spsCC_abiebals_2020_v2_20260119.tif", "SCANFI_spsCC_betupapy_2020_v2_20260119.tif")
  )
  fakeListing[, url := paste0("https://example.invalid/", name)]
  testthat::local_mocked_bindings(
    listGoogleDriveFolder = function(...) fakeListing,
    .package = "reproducible"
  )
  testthat::local_mocked_bindings(
    prepInputs = function(..., cropTo = NULL) {
      cropTo
    }
  )

  sppEquiv <- data.table::data.table(SCANFI = c("abiebals", "betupapy"), Species = c("Abie_bal", "Betu_pap"))
  rast <- function(nrows) {
    r <- terra::rast(nrows = nrows, ncols = 6, xmin = 0, xmax = 600, ymin = 0, ymax = 600,
                     crs = "EPSG:3978")
    terra::values(r) <- 50
    r
  }
  load <- function(cropTo) {
    LandR::loadSCANFISpeciesLayers(
      dPath = dPath, cropTo = cropTo, sppEquiv = sppEquiv, studyAreaName = "x",
      cachePath = cacheDir
    )
  }

  ## the mocked prepInputs() returns its cropTo, so the result shows which grid it came from
  expect_equal(terra::nrow(suppressMessages(load(rast(6)))), 6)
  expect_equal(terra::nrow(suppressMessages(load(rast(7)))), 7)

  ## identical targets (a fresh but equal object) still hit the cache
  expect_message(load(rast(7)), "Loaded! Cached result")
})
