## Forest land is what separates a temporarily open forest (a burn, a cutblock, a fire-scoured
## barren) from land that is permanently open (a bog, a rock barren that never carried trees).
## These tests use small in-memory rasters, so nothing is downloaded.

lccRas <- function(values) {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3000, ymin = 0, ymax = 2000,
                   crs = "EPSG:3978")
  terra::values(r) <- values
  r
}

##              burn   bog   lake  mature  cutblock  barren
lccY     <- c(  50,    50,    20,    210,     100,      33)
faoLayer <- c(   2,     0,     0,      1,       1,       0)
lcc2015  <- c( 210,    50,    20,    210,     220,      33)

test_that("forestLandMask(): FAO codes 1 and 2 are both forest land", {
  mask <- forestLandMask(faoRas = lccRas(faoLayer))
  expect_identical(as.vector(terra::values(mask)), c(1, 0, 0, 1, 1, 0))
})

test_that("forestLandMask(): a pixel treed in any scanned year is forest land", {
  mask <- forestLandMask(lccList = list(lccRas(lccY), lccRas(lcc2015)))
  ## the burn is treed in 2015 and the cutblock too; the bog, lake and barren never are
  expect_identical(as.vector(terra::values(mask)), c(1, 0, 0, 1, 1, 0))
})

test_that("forestLandMask(): the union is taken over both sources", {
  ## a stand that had recovered by the FAO layer's year is code 1, which the old rule
  ## (code 2 only) discarded; either source alone can also come back empty
  faoOnly <- forestLandMask(faoRas = lccRas(rep(0, 6)))
  expect_identical(as.vector(terra::values(faoOnly)), rep(0, 6))

  both <- forestLandMask(lccList = list(lccRas(lcc2015)), faoRas = lccRas(faoLayer))
  expect_identical(as.vector(terra::values(both)), c(1, 0, 0, 1, 1, 0))
})

test_that("forestLandMask(): NAs are not forest land, and do not propagate", {
  mask <- forestLandMask(lccList = list(lccRas(c(210, NA, NA, 210, NA, NA))),
                         faoRas = lccRas(c(NA, NA, 0, 1, 2, NA)))
  expect_false(any(is.na(terra::values(mask))))
  expect_identical(as.vector(terra::values(mask)), c(1, 0, 0, 1, 1, 0))
})

test_that("forestLandMask(): needs at least one source", {
  expect_error(forestLandMask(), "at least one of")
})

test_that(".applyForestLand(): non-treed forest land becomes the disturbed code", {
  mask <- forestLandMask(lccList = list(lccRas(lcc2015)), faoRas = lccRas(faoLayer))
  out <- .applyForestLand(lccRas(lccY), mask, disturbedCode = 240)
  ## burn and cutblock are forest land without trees; bog, lake and barren are not forest land
  expect_identical(as.vector(terra::values(out)), c(240, 50, 20, 210, 240, 33))
})

test_that(".applyForestLand(): barren and rock are eligible when the record says forest land", {
  ## a barren pixel that the time series shows as treed in other years is exposed by
  ## something temporary, such as a severe fire
  mask <- forestLandMask(lccList = list(lccRas(c(210, 50, 20, 210, 220, 210))))
  out <- .applyForestLand(lccRas(c(33, 50, 20, 210, 100, 32)), mask, disturbedCode = 240)
  expect_identical(as.vector(terra::values(out)), c(240, 50, 20, 210, 240, 240))
})

test_that(".applyForestLand(): treed pixels keep their class, and NAs stay unset", {
  mask <- forestLandMask(faoRas = lccRas(rep(1, 6)))
  out <- .applyForestLand(lccRas(c(210, 220, 230, 81, 50, NA)), mask, disturbedCode = 240)
  vals <- as.vector(terra::values(out))
  expect_identical(vals[1:5], c(210, 220, 230, 81, 240))
  expect_true(is.na(vals[6])) ## terra writes an unset cell as NaN
})

test_that(".applyForestLand(): convertibleClasses narrows what may change", {
  mask <- forestLandMask(faoRas = lccRas(rep(1, 6)))
  out <- .applyForestLand(lccRas(c(33, 50, 20, 210, 100, 32)), mask,
                          convertibleClasses = c(40, 50, 80, 100), disturbedCode = 240)
  expect_identical(as.vector(terra::values(out)), c(33, 240, 20, 210, 240, 32))
})

test_that(".defaultForestLandYears(): one per decade plus the most recent", {
  expect_identical(.defaultForestLandYears(1984:2022), c(1985L, 1995L, 2005L, 2015L, 2022L))
  ## SCANFI V2 publishes every five years
  expect_identical(.defaultForestLandYears(seq(1985L, 2025L, 5L)), c(1985L, 1995L, 2005L, 2015L, 2025L))
  ## SCANFI V1 has only three years, and they are what gets used
  expect_identical(.defaultForestLandYears(c(2000L, 2010L, 2020L)), c(2000L, 2010L, 2020L))
})

test_that("prepInputs_FAO_forest(): only published years are accepted", {
  expect_error(prepInputs_FAO_forest(year = 2010), "published for 2019 and 2022")
})
