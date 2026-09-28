## SCANFI v3 land cover: crosswalk, synthetic-raster recode, forest-land handling of burn
## scars, URL construction and access-error messaging. The windowed /vsicurl read itself is
## network-dependent and skipped offline/non-interactively, like the other SCANFI network
## tests in this package (see test-prepInputs_SCANFI_structure.R).

test_that("scanfiV3ToCanadaLCC covers codes 1-20 and 255, exactly once each", {
  expect_identical(sort(scanfiV3ToCanadaLCC$scanfiV3), sort(c(1:20, 255L)))
  expect_false(any(duplicated(scanfiV3ToCanadaLCC$scanfiV3)))

  ## burn scars (4) keep their own code, not a Canada LCC code, and stay flammable
  expect_identical(scanfiV3ToCanadaLCC$lcc[scanfiV3ToCanadaLCC$scanfiV3 == 4], 60)
  expect_true(is.na(scanfiV3ToCanadaLCC$lcc[scanfiV3ToCanadaLCC$scanfiV3 == 255]))

  ## low and tall shrubs (7, 8) both collapse to the Canada LCC shrub code
  expect_identical(scanfiV3ToCanadaLCC$lcc[scanfiV3ToCanadaLCC$scanfiV3 %in% c(7, 8)], c(50, 50))

  ## cropland, urban, road (17-19) have no Canada LCC analogue here
  expect_identical(scanfiV3ToCanadaLCC$lcc[scanfiV3ToCanadaLCC$scanfiV3 %in% 17:19], c(0, 0, 0))

  ## 60 is not already a treed class, so a burn scar off forest land stays non-forest
  expect_false(60 %in% LandR:::.treedLCCClasses)
})

test_that(".applySCANFIv3Crosswalk() recodes a synthetic v3 raster to Canada LCC codes", {
  r <- terra::rast(nrows = 1, ncols = 21, xmin = 0, xmax = 21000, ymin = 0, ymax = 1000,
                   crs = "EPSG:3979")
  terra::values(r) <- c(1:20, 255L)

  out <- LandR:::.applySCANFIv3Crosswalk(r)

  expect_identical(
    as.vector(terra::values(out)),
    c(20, 30, 33, 60, 40, 100, 50, 50, 220, 230,
      210, 210, 210, 210, 210, 210, 0, 0, 0, 31, NA)
  )
})

test_that("burn scars become the disturbed code on forest land but stay 60 off it", {
  r <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2000, ymin = 0, ymax = 1000,
                   crs = "EPSG:3979")
  ## a burn scar that is forest land (grew trees in another scanned year), and one that
  ## is not (e.g. a bog margin the record never shows as treed)
  lcc <- r
  terra::values(lcc) <- c(60, 60)
  mask <- r
  terra::values(mask) <- c(1, 0)

  out <- LandR:::.applyForestLand(lcc, mask, disturbedCode = 240)

  expect_identical(as.vector(terra::values(out)), c(240, 60))
})

test_that(".scanfiV3Url() builds the year's COG url and rejects years v3 doesn't publish", {
  expect_identical(
    LandR:::.scanfiV3Url(2020),
    paste0(
      "https://download-telecharger.services.geo.ca/pub/nrcan_rncan/Forests_Foret/",
      "SCANFI/v3/cog_SCANFI_landcover_2020_v3_20260528.tif"
    )
  )
  expect_error(LandR:::.scanfiV3Url(1900), "does not exist")
})

test_that("prepInputs_SCANFI_LCC_FAO() rejects a year V3 does not publish", {
  expect_error(
    prepInputs_SCANFI_LCC_FAO(year = 1900, dataVersion = "V3"),
    "does not exist"
  )
})

test_that(".isSCANFIv3AccessError() recognizes https access failures but not other errors", {
  expect_true(LandR:::.isSCANFIv3AccessError(simpleError("HTTP response code said 403")))
  expect_true(LandR:::.isSCANFIv3AccessError(simpleError("curl error: Could not open connection")))
  expect_true(LandR:::.isSCANFIv3AccessError(simpleError("Operation timed out")))

  expect_false(LandR:::.isSCANFIv3AccessError(simpleError("subscript out of bounds")))
  expect_false(LandR:::.isSCANFIv3AccessError(simpleError("non-numeric argument")))
})

test_that(".withSCANFIv3Access explains an access failure and keeps the original error", {
  err <- tryCatch(
    LandR:::.withSCANFIv3Access(
      stop(simpleError("HTTP response code said 403")),
      what = "the SCANFI v3 land cover map",
      dataYear = 2020
    ),
    error = function(e) conditionMessage(e)
  )
  flat <- gsub("[[:space:]]+", " ", err)

  expect_match(flat, "Could not download the SCANFI v3 land cover map (V3 2020)", fixed = TRUE)
  expect_match(flat, "User-Agent", fixed = TRUE)
  expect_match(flat, "ftp://ftp.maps.canada.ca", fixed = TRUE)
  expect_match(flat, "HTTP response code said 403", fixed = TRUE)
})

test_that(".withSCANFIv3Access leaves unrelated errors alone", {
  err <- tryCatch(
    LandR:::.withSCANFIv3Access(stop("something else broke"), what = "the SCANFI v3 land cover map"),
    error = function(e) conditionMessage(e)
  )
  expect_identical(err, "something else broke")
})

test_that("prepInputs_SCANFI_LCC_FAO(dataVersion = 'V3') reads a windowed study area", {
  testthat::skip_if_not(interactive(), "network-dependent, run interactively")
  testthat::skip_if_offline()

  ## a small window near Kenora, ON, in the SCANFI v3 CRS (EPSG:3979)
  to <- terra::rast(
    xmin = -1200000, xmax = -1195000, ymin = 700000, ymax = 705000,
    resolution = 30, crs = "EPSG:3979"
  )

  t0 <- Sys.time()
  out <- prepInputs_SCANFI_LCC_FAO(
    year = 2020, dataVersion = "V3", to = to, forestLandFrom = "lccYears",
    forestLandYears = 2020
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  message("prepInputs_SCANFI_LCC_FAO(dataVersion = 'V3') windowed read took ", elapsed, "s")

  expect_s4_class(out, "SpatRaster")
  expect_true(all(terra::values(out) %in% c(scanfiV3ToCanadaLCC$lcc, 240, NA)))
})
