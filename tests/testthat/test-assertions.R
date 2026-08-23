## a categorical ecoregionMap plus the tables whose ecoregionGroups must agree with it
ergFixture <- function(ergs = c("1_11", "1_12", "2_11"), withCats = TRUE,
                       ergCol = "ecoregionGroup") {
  ras <- terra::rast(
    nrows = 6, ncols = 6, xmin = 0, xmax = 600, ymin = 0, ymax = 600, crs = "EPSG:3978"
  )
  terra::values(ras) <- rep(seq_along(ergs), length.out = terra::ncell(ras))
  if (withCats) {
    levs <- data.frame(ID = seq_along(ergs), ergs)
    colnames(levs)[2] <- ergCol
    levels(ras) <- levs
  }
  list(
    ecoregionMap = ras,
    cohortData = data.table(ecoregionGroup = factor(rep(ergs, 2L)), B = 1),
    speciesEcoregion = data.table(ecoregionGroup = factor(ergs)),
    minRelativeB = data.table(ecoregionGroup = factor(ergs))
  )
}

test_that("assertERGs passes when the ecoregionGroups all match", {
  f <- ergFixture()
  expect_no_error(assertERGs(
    f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
    doAssertion = TRUE
  ))
})

test_that("assertERGs detects mismatched ecoregionGroups", {
  f <- ergFixture()
  f$speciesEcoregion <- data.table(ecoregionGroup = factor(c("1_11", "1_12", "3_11")))
  expect_error(
    capture.output(suppressMessages(assertERGs(
      f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
      doAssertion = TRUE
    ))),
    "exactly the same"
  )
})

## https://github.com/PredictiveEcology/LandR/issues/190
test_that("assertERGs errors informatively when ecoregionMap has no categories", {
  f <- ergFixture(withCats = FALSE)
  expect_snapshot(error = TRUE, assertERGs(
    f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
    doAssertion = TRUE
  ))
})

test_that("assertERGs errors informatively when ecoregionMap has no ecoregionGroup column", {
  f <- ergFixture(ergCol = "ecoregion")
  expect_snapshot(error = TRUE, assertERGs(
    f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
    doAssertion = TRUE
  ))
})

test_that("assertERGs does nothing when doAssertion is FALSE", {
  f <- ergFixture(withCats = FALSE)
  expect_no_error(assertERGs(
    f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
    doAssertion = FALSE
  ))
})
