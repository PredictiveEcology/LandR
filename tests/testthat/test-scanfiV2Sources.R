## SCANFI v2 files are fetched from the arbutus mirror by url, not by Google Drive id.
## No network except the last test.

arbutusPrefix <- "^https://object-arbutus[.]alliancecan[.]ca/529feda02d404836b0e91b1d360b4df0:predictiveecology/SCANFI_v2/"

test_that("every (layer, year) of SCANFI v2 resolves to an arbutus url named after its file", {
  for (layer in names(LandR:::.scanfiV2Layers)) {
    for (y in LandR:::.scanfi_v2_years) {
      src <- LandR:::.scanfiV2Source(layer, y)
      expect_match(src$url, paste0(arbutusPrefix, y, "/", src$targetFile, "$"))
      expect_match(src$targetFile, paste0("_", y, "_v2_[0-9]+[.]tif$"))
    }
  }
  expect_error(LandR:::.scanfiV2Source("age", 2021), "available for 1985")
})

test_that("the species files of each SCANFI v2 year come from the manifest", {
  for (y in LandR:::.scanfi_v2_years) {
    sp <- LandR:::.scanfiV2SpeciesFiles(y)
    expect_gt(nrow(sp), 10)
    expect_true(all(grepl(paste0("_", y, "_v2_"), sp$name)))
    expect_match(sp$url, arbutusPrefix)
  }
})

test_that("the default SCANFI v2 requests pass arbutus urls to prepInputs, without Google Drive", {
  cacheDir <- withr::local_tempdir()
  withr::local_options(reproducible.cachePath = cacheDir, reproducible.useMemoise = FALSE)
  seen <- character()
  testthat::local_mocked_bindings(
    prepInputs = function(..., url = NULL) {
      seen <<- c(seen, url)
      terra::rast(nrows = 2, ncols = 2, vals = 50)
    }
  )
  testthat::local_mocked_bindings(
    listGoogleDriveFolder = function(...) stop("must not list a Drive folder"),
    .package = "reproducible"
  )
  dPath <- withr::local_tempdir()
  suppressMessages(prepRawBiomassMap(dataSource = "SCANFI", dataYear = 2010, destinationPath = dPath))
  suppressMessages(prepInputsStandAgeMap(dataSource = "SCANFI", dataYear = 2010, destinationPath = dPath,
                                         fireURL = NA))
  suppressMessages(prepInputs_SCANFI_structure("height", year = 2010, destinationPath = dPath))
  expect_length(seen, 3)
  expect_match(seen, arbutusPrefix)
  expect_match(seen, "/2010/")

  sppEquiv <- data.table::data.table(SCANFI = "PICE_MAR", Species = "Pice_mar")
  seen <- character()
  suppressMessages(loadSCANFISpeciesLayers(dPath = dPath, sppEquiv = sppEquiv, year = 2010,
                                           cachePath = cacheDir))
  expect_match(seen, paste0(arbutusPrefix, "2010/SCANFI_spsCC_PICE_MAR_2010_v2_"))
})

test_that("the land cover source of V2 is the arbutus file and V1 is unchanged", {
  expect_match(LandR:::.scanfiLCCFAOSource(2020, "V2")$url, arbutusPrefix)
  expect_equal(LandR:::.scanfiLCCFAOSource(2020, "V2")$targetFile,
               "SCANFI_att_nfiLandcover_CanadaLCCclassCodes_2020_v2_20260119.tif")
  expect_match(LandR:::.scanfiLCCFAOSource(2020, "V1")$url, "^https://drive[.]google[.]com/file/d/")
})

test_that("an arbutus url answers a HEAD request", {
  skip_if_offline("object-arbutus.alliancecan.ca")
  skip_on_cran()
  url <- LandR:::.scanfiV2Source("age", 2020)$url
  resp <- httr2::req_perform(httr2::req_method(httr2::request(url), "HEAD"))
  expect_equal(httr2::resp_status(resp), 200L)
})
