testthat::test_that("prepInputsStandAgeMap works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("googledrive")
  td <- withr::local_tempdir("dest_")
  testthat::skip_if_not(interactive())

  #this default randomStudyArea is not centred on a place with much fire...
  goodPoly <- randomStudyArea(size = 5e10, seed = 7)
  goodRas <- terra::rast(goodPoly, vals = 1, res = c(250, 250))
  goodRas <- terra::mask(goodRas, goodPoly)


  out <- prepInputsStandAgeMap(dataSource = "SCANFI",
                               dataYear = 2000, studyArea = goodPoly,
                               rasterToMatch = goodRas, destinationPath = td,
                               overwrite = TRUE)
  expect_true(length(attr(out, "imputedPix")) > 0)

  out2 <- prepInputsStandAgeMap(dataSource = "KNN",
                                dataYear = 2011, studyArea = goodPoly,
                                rasterToMatch = goodRas, destinationPath = td,
                                overwrite = TRUE, firePerimeters = NULL,
                                fireURL = NULL)
  expect_true(length(attr(out2, "imputedPix")) != length(attr(out, "imputedPix")))

  out3 <- prepInputsStandAgeMap(dataSource = "KNN",
                                dataYear = 2011, studyArea = goodPoly,
                                rasterToMatch = goodRas, destinationPath = td)
  expect_true(length(attr(out3, "imputedPix")) != length(attr(out2, "imputedPix")))

  out4 <- prepInputsStandAgeMap(dataSource = "SCANFI",
                                dataYear = 2020, studyArea = goodPoly,
                                rasterToMatch = goodRas, destinationPath = td,
                                overwrite = TRUE)
  expect_true(length(attr(out4, "imputedPix")) != length(attr(out, "imputedPix")))


})
