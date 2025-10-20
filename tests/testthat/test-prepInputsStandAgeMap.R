testthat::test_that("prepInputsStandAgeMap work", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  td <- withr::local_tempdir("dest_")

  goodPoly <- randomStudyArea(size = 6e8, seed = 5)
  goodRas <- terra::rast(goodPoly, vals = 1, res = c(250, 250))
  goodRas <- terra::mask(goodRas, goodPoly)


  out <- prepInputsStandAgeMap(dataSource = "SCANFI",
                               dataYear = 2000, studyArea = goodPoly,
                               rasterToMatch = goodRas, destinationPath = td,
                               overwrite = TRUE)
  expect_true(length(attr(out, "imputedPix")) > 0)

  out2 <- prepInputsStandAgeMap(dataSource = "SCANFI",
                                dataYear = 2020, studyArea = goodPoly,
                                rasterToMatch = goodRas, destinationPath = td,
                                overwrite = TRUE)
  expect_true(out2 != out)

  out3 <- prepInputsStandAgeMap(dataSource = "SCANFI",
                                dataYear = 2020, studyArea = goodPoly,
                                rasterToMatch = goodRas, destinationPath = td,
                                overwrite = TRUE, firePerimeters = NULL,
                                fireURL = NULL)


})
