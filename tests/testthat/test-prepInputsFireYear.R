testthat::test_that("prepInputs fire year works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  td <- withr::local_tempdir("dest_")

  badExtent <- terra::ext(-80.7421995, -75.816238, 44.3156278, 46.971853)
  badPoly <- terra::vect(badExtent)
  terra::crs(badPoly) <- "epsg:4269"
  badRTM <- terra::rast(
    terra::ext(badPoly),
    crs = terra::crs(badPoly),
    resolution = c(0.01, 0.01),
    vals = 1
  )
  furl <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
  testthat::expect_no_error({
    suppressWarnings({
      prepInputsFireYear(
        rasterToMatch = badRTM, maskTo = badPoly,
        destinationPath = td,
        url = furl,
        fireField = "YEAR"
      )
    })
  })

  goodPoly <- randomStudyArea(size = 6e8, seed = 5)
  goodRas <- terra::rast(goodPoly, vals = 1, res = c(250, 250))
  goodRas <- terra::mask(goodRas, goodPoly)

  goodFire <- suppressWarnings({
    prepInputsFireYear(
      rasterToMatch = goodRas,
      earliestYear = 2000, #limit postProcessing
      maskTo = goodPoly,
      destinationPath = td,
      url = furl,
      fireField = "YEAR"
    )
  })
  rasNAs <- as.vector(goodRas)
  fireNAs <- as.vector(goodFire)
  expect_true(all(is.na(fireNAs[is.na(rasNAs)]))) #ensures postProcess was correct
  #previously was not masking correct

})
