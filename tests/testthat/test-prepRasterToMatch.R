test_that("prepRasterToMatch works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  testthat::skip_if_not_installed("withr")

  td <- withr::local_tempdir("dest_")

  withr::local_package("terra")

  # get SAs from PredictiveEcology.org example Chapter 14
  originalcrs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  smallExtent <- c(xmin = -104.757, xmax = -104.48835, ymin = 55.68663, ymax = 55.94491)
  studyAreaS <- terra::vect(terra::ext(smallExtent))
  studyAreaS <- terra::vect(terra::geom(studyAreaS), "polygons",
    crs = originalcrs,
    atts = data.frame(id = 1:length(studyAreaS))
  )

  largeExtent <- c(xmin = -104.757, xmax = -104.2197, ymin = 55.68663, ymax = 56.20319)

  studyAreaL <- terra::vect(terra::ext(largeExtent))
  studyAreaL <- terra::vect(terra::geom(studyAreaL), "polygons",
    crs = originalcrs,
    atts = data.frame(id = 1:length(studyAreaL))
  )

  Biomass_corecrs <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
                           "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  studyAreaL <- terra::project(studyAreaL, Biomass_corecrs)
  studyAreaS <- terra::project(studyAreaS, Biomass_corecrs)

  rtmE <- terra::rast(terra::ext(c(-105, -104, 55, 56)),
                      resolution = c(0.01, 0.01), vals = 1)

  ## warn that rasterToMatchLarge and rasterToMatch are both missing
  testthat::expect_warning(
    RTMs <- prepRasterToMatch(
      studyArea = studyAreaS,
      studyAreaLarge = studyAreaL,
      rasterToMatch = NULL,
      rasterToMatchLarge = NULL,
      templateRas = rtmE,
      destinationPath = td,
      studyAreaName = "foo"
    )
  )
})
