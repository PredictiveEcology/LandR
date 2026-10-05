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
