## Regression test (offline): the fire-year raster must have exactly the geometry of
## rasterToMatch. replaceAgeInFires() indexes the stand age map (on rasterToMatch) with it; with
## `studyArea` (a polygon) in the dots, the final postProcess cropped to the polygon's bounding
## box and lost the grid's edge rows: "[`[<-`] lengths of cells and values do not match"
## (BC TSA13, FireSense ELF 14.4).
test_that("prepInputsFireYear returns rasterToMatch's grid when given a studyArea polygon", {
  skip_if_not_installed("terra")
  td <- withr::local_tempdir("fireGrid_")
  crs <- "EPSG:3978"
  ## a 240 m grid with whole rows and columns beyond the polygon's bounding box
  rtm <- terra::rast(terra::ext(0, 2400, 0, 2400), resolution = 240, crs = crs, vals = 1)
  sa <- terra::vect("POLYGON ((500 500, 1900 550, 1850 1900, 550 1850, 500 500))", crs = crs)
  fires <- terra::vect(c("POLYGON ((600 600, 1200 600, 1200 1200, 600 1200, 600 600))",
                         "POLYGON ((1300 1300, 1800 1300, 1800 1800, 1300 1800, 1300 1300))"), crs = crs)
  fires$YEAR <- c(1995, 2010)
  src <- withr::local_tempdir("fireSrc_")
  terra::writeVector(fires, file.path(src, "fires.shp"))
  zipped <- withr::with_dir(src, utils::zip("fires.zip", list.files(src, "^fires\\."), flags = "-q"))
  skip_if_not(identical(zipped, 0L), "no zip utility")

  out <- suppressWarnings(suppressMessages(
    prepInputsFireYear(rasterToMatch = rtm, studyArea = sa, destinationPath = td,
                       url = paste0("file://", file.path(src, "fires.zip")), fireField = "YEAR",
                       useCache = FALSE)))
  expect_true(terra::compareGeom(out, rtm, stopOnError = FALSE))
  ## still masked to the study area, and the fires are there
  expect_true(all(is.na(terra::values(terra::mask(out, sa, inverse = TRUE)))))
  expect_setequal(stats::na.omit(unique(terra::values(out, mat = FALSE))), c(1995, 2010))
})

## A pixel that burned twice keeps only its latest fire year (rasterize fun = max). For an earlier
## data year the later fire hid the earlier one, so replaceAgeInFires() skipped the pixel and it kept
## its satellite age. `latestYear` drops the fires after the data year before rasterizing.
test_that("prepInputsFireYear(latestYear) keeps the latest fire at or before that year", {
  skip_if_not_installed("terra")
  td <- withr::local_tempdir("fireLatest_")
  crs <- "EPSG:3978"
  rtm <- terra::rast(terra::ext(0, 2400, 0, 2400), resolution = 240, crs = crs, vals = 1)
  fires <- terra::vect(c("POLYGON ((600 600, 1200 600, 1200 1200, 600 1200, 600 600))",
                         "POLYGON ((600 600, 1200 600, 1200 1200, 600 1200, 600 600))",
                         "POLYGON ((1300 1300, 1800 1300, 1800 1800, 1300 1800, 1300 1300))"), crs = crs)
  fires$YEAR <- c(1980, 2015, 2010)
  src <- withr::local_tempdir("fireSrc_")
  terra::writeVector(fires, file.path(src, "fires.shp"))
  zipped <- withr::with_dir(src, utils::zip("fires.zip", list.files(src, "^fires\\."), flags = "-q"))
  skip_if_not(identical(zipped, 0L), "no zip utility")
  makeFires <- function(...) {
    suppressWarnings(suppressMessages(
      prepInputsFireYear(rasterToMatch = rtm, destinationPath = td, fireField = "YEAR",
                         url = paste0("file://", file.path(src, "fires.zip")), useCache = FALSE, ...)))
  }

  twice <- terra::cellFromXY(rtm, cbind(900, 900))
  noLimit <- makeFires()
  expect_equal(noLimit[twice][[1]], 2015)
  limited <- makeFires(latestYear = 1990)
  expect_equal(limited[twice][[1]], 1980)
  ## the 2010 fire is later than the limit, so its pixels have no fire
  expect_true(is.na(limited[terra::cellFromXY(rtm, cbind(1500, 1500))][[1]]))
  expect_equal(noLimit[terra::cellFromXY(rtm, cbind(1500, 1500))][[1]], 2010)

  ## stand age at 1990 inside the 1980 burn is 10
  age <- replaceAgeInFires(terra::rast(rtm, vals = 80L), limited, startTime = 1990)
  expect_equal(age[twice][[1]], 10L)
  ## no fires at all before the limit gives an empty raster, not an error
  expect_true(all(is.na(terra::values(makeFires(latestYear = 1960)))))
})
