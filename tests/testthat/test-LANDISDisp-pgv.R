## Supplying `pgv` gives the same result as reading it from pixelGroupMap.

test_that("LANDISDisp(pgv =) matches the raster path, including NA cells", {
  skip_if_not_installed("terra")
  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L, successionTimestep = 10L)
  pgm <- fix$pixelGroupMap
  ## mask some cells to NA so the NA path is exercised
  v <- as.integer(terra::values(pgm, mat = FALSE))
  v[seq(7L, length(v), by = 13L)] <- NA_integer_
  terra::values(pgm) <- v
  fix$pixelGroupMap <- pgm

  for (useCpp in c(TRUE, FALSE)) {
    viaRaster <- runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = useCpp)
    viaPgv <- runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = useCpp, pgv = v)
    expect_gt(nrow(viaRaster), 0L)
    expect_equal(viaPgv, viaRaster, info = paste("useCpp =", useCpp))
  }
})

test_that("LANDISDisp errors when pgv has the wrong length", {
  skip_if_not_installed("terra")
  fix <- makeLANDISDispFixture(size = "tiny", fixtureSeed = 11L)
  expect_error(
    runLANDISDispOnFixture(fix, runSeed = 1L, pgv = c(1L, 2L)),
    "ncell"
  )
})
