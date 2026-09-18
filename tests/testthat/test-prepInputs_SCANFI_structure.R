test_that("prepInputs_SCANFI_structure rejects years and versions it cannot serve", {
  ## No network: these are the argument checks, which run before any download.
  expect_error(prepInputs_SCANFI_structure("height", year = 2021),
               "available for 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025 only")
  expect_error(prepInputs_SCANFI_structure("closure", year = 1984),
               "available for 1985")
  expect_error(prepInputs_SCANFI_structure("height", year = 2020, dataVersion = "V1"),
               "published for V2 only")
  expect_error(prepInputs_SCANFI_structure("biomass", year = 2020), "'arg' should be one of")
})

test_that("prepInputs_SCANFI_structure has an id for every V2 year of both attributes", {
  ids <- LandR:::.scanfiStructureIds
  expect_setequal(names(ids), c("height", "closure"))
  for (att in names(ids)) {
    expect_setequal(names(ids[[att]]), as.character(LandR:::.scanfi_v2_years))
    expect_true(all(nzchar(ids[[att]])))
  }
  ## ids must be distinct: a copy-paste between years would silently serve the wrong layer
  expect_equal(anyDuplicated(unlist(ids, use.names = FALSE)), 0L)
})

test_that("prepInputs_SCANFI_structure downloads height and closure", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_not_installed("withr")
  skip_if_not_installed("googledrive")
  testthat::skip_if_not(googledrive::drive_has_token(), "No Drive token")

  withr::local_package("terra")
  dPath <- withr::local_tempdir("inputs_")
  withr::local_options(list(
    reproducible.destinationPath = dPath,
    reproducible.rasterRead = "terra::rast",
    reproducible.useTerra = TRUE
  ))
  sa <- LandR::randomStudyArea(size = 1e8, seed = 5)
  rtm <- terra::mask(terra::rast(terra::vect(sa), res = 250, vals = 1), terra::vect(sa))

  h <- prepInputs_SCANFI_structure("height", year = 2020, to = rtm)
  expect_s4_class(h, "SpatRaster")
  expect_true(LandR::.compareRas(h, rtm))
  expect_true(all(terra::values(h, mat = FALSE) >= 0, na.rm = TRUE))

  cc <- prepInputs_SCANFI_structure("closure", year = 2020, to = rtm)
  expect_s4_class(cc, "SpatRaster")
  ## closure is a percent
  expect_true(all(terra::values(cc, mat = FALSE) <= 100, na.rm = TRUE))
})
