## standAgeMapGenerator() and vegTypeMapGenerator() must leave the caller's cohortData alone.
## NRV_summary passes sim$cohortData straight in, so a `:=` inside standAgeMapGenerator()
## added a `weightedAge` column to the live simulation's cohortData.

mapsCohortData <- function() {
  data.table::data.table(
    pixelGroup = c(1L, 1L, 2L, 2L, 2L, 3L),
    speciesCode = factor(c("Pice_gla", "Popu_tre", "Pice_gla", "Pinu_ban", "Popu_tre", "Pinu_ban")),
    age = c(95L, 42L, 130L, 61L, 17L, 8L),
    B = c(3000L, 1000L, 500L, 2500L, 2000L, 150L)
  )
}

mapsPixelGroupMap <- function() {
  ## pixelGroup 4 has no cohorts, and one pixel is NA
  r <- terra::rast(terra::ext(0, 3, 0, 3), resolution = 1)
  terra::values(r) <- c(1, 1, 2, 2, 3, 3, NA, 4, 2)
  r
}

mapsSppEquiv <- function() {
  data.table::data.table(
    LandR = c("Pice_gla", "Pinu_ban", "Popu_tre"),
    Type = c("Conifer", "Conifer", "Deciduous")
  )
}

test_that("standAgeMapGenerator() does not modify cohortData or pixelGroupMap", {
  for (w in list("biomass", NA)) {
    cd <- mapsCohortData()
    cdBefore <- data.table::copy(cd)
    pgm <- mapsPixelGroupMap()
    pgmNames <- names(pgm)
    pgmValues <- terra::values(pgm, mat = FALSE)

    standAgeMapGenerator(cd, pgm, weight = w)

    expect_identical(cd, cdBefore, label = paste0("cohortData after weight = ", w))
    expect_identical(names(pgm), pgmNames)
    expect_identical(terra::values(pgm, mat = FALSE), pgmValues)
  }
})

test_that("standAgeMapGenerator() returns ages rounded down to the decade", {
  ## biomass-weighted: pg1 (95*3000 + 42*1000) / 4000 = 81.75 -> 80
  ##                   pg2 (130*500 + 61*2500 + 17*2000) / 5000 = 50.3 -> 50
  ##                   pg3 8 -> 0
  ## unweighted max:   pg1 95 -> 90; pg2 130 -> 130; pg3 8 -> 0
  ## pixelGroup 4 has no cohorts, so it is NA like the NA pixel.
  byB <- standAgeMapGenerator(mapsCohortData(), mapsPixelGroupMap(), weight = "biomass")
  byMax <- standAgeMapGenerator(mapsCohortData(), mapsPixelGroupMap(), weight = NA)

  expect_s4_class(byB, "SpatRaster")
  expect_identical(terra::values(byB, mat = FALSE), c(80, 80, 50, 50, 0, 0, NA, NA, 50))
  expect_identical(terra::values(byMax, mat = FALSE), c(90, 90, 130, 130, 0, 0, NA, NA, 130))

  ## a stale `weightedAge` column in the input (as the old by-reference version left behind)
  ## does not change the result
  cd <- mapsCohortData()
  cd[["weightedAge"]] <- 999
  expect_identical(
    terra::values(standAgeMapGenerator(cd, mapsPixelGroupMap(), weight = "biomass"), mat = FALSE),
    terra::values(byB, mat = FALSE)
  )
})

test_that("vegTypeMapGenerator() does not modify cohortData or pixelGroupMap", {
  for (mt in 0:2) {
    cd <- mapsCohortData()
    cdBefore <- data.table::copy(cd)
    pgm <- mapsPixelGroupMap()
    terra::values(pgm) <- c(1, 1, 2, 2, 3, 3, NA, 1, 2) ## every pixelGroup has cohorts
    pgmNames <- names(pgm)
    pgmValues <- terra::values(pgm, mat = FALSE)

    vtm <- suppressMessages(vegTypeMapGenerator(
      cd,
      pixelGroupMap = pgm,
      vegLeadingProportion = 0.75,
      mixedType = mt,
      sppEquiv = mapsSppEquiv(),
      sppEquivCol = "LandR",
      doAssertion = TRUE
    ))

    expect_s4_class(vtm, "SpatRaster")
    expect_identical(cd, cdBefore, label = paste0("cohortData after mixedType = ", mt))
    expect_identical(names(pgm), pgmNames)
    expect_identical(terra::values(pgm, mat = FALSE), pgmValues)
  }
})
