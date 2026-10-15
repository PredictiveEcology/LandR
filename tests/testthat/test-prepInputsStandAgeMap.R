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

## Offline: the data year limits the fires used to correct stand age. A pixel burned in 1980 and
## 2015 holds 2015 in the unlimited fire raster, so for dataYear 2011 it kept its satellite age.
test_that("prepInputsStandAgeMap uses the latest fire at or before dataYear", {
  skip_if_not_installed("terra")
  crs <- "EPSG:3978"
  rtm <- terra::rast(terra::ext(0, 2400, 0, 2400), resolution = 240, crs = crs, vals = 1)
  src <- withr::local_tempdir("ageSrc_")
  terra::writeRaster(terra::rast(rtm, vals = 80L), file.path(src, "age.tif"), datatype = "INT2U")
  square <- "POLYGON ((600 600, 1200 600, 1200 1200, 600 1200, 600 600))"
  fires <- terra::vect(c(square, square), crs = crs)
  fires$YEAR <- c(1980, 2015)
  terra::writeVector(fires, file.path(src, "fires.shp"))
  zipped <- withr::with_dir(src, utils::zip("fires.zip", list.files(src, "^fires\\."), flags = "-q"))
  skip_if_not(identical(zipped, 0L), "no zip utility")

  out <- suppressWarnings(suppressMessages(
    prepInputsStandAgeMap(dataSource = "KNN", dataYear = 2011, rasterToMatch = rtm,
                          ageURL = paste0("file://", file.path(src, "age.tif")),
                          fireURL = paste0("file://", file.path(src, "fires.zip")),
                          destinationPath = withr::local_tempdir("dest_"),
                          maskWithRTM = FALSE, useCache = FALSE)))
  expect_equal(out[terra::cellFromXY(rtm, cbind(900, 900))][[1]], 31L) # 2011 - 1980
})
