testthat::test_that("test prepRawBiomassMap", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()

  skip_if_not_installed("withr")
  skip_if_not_installed("googledrive")
  testthat::skip_if_not(googledrive::drive_has_token(), "No Drive token")

  withr::local_package("reproducible")
  withr::local_package("SpaDES.tools")
  withr::local_package("terra")
  withr::local_package("sf")

  dPath <- withr::local_tempdir("inputs_")
  cPath <- withr::local_tempdir("cache_")

  withr::local_options(list(
    reproducible.destinationPath = dPath,
    reproducible.inputPaths = NULL,
    reproducible.overwrite = TRUE,
    reproducible.rasterRead = "terra::rast",
    reproducible.useTerra = TRUE
  ))

  biomassURL <- paste0(
    "http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
    "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
    "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
  )
  studyArea <- randomStudyArea()
  RTM <- rast(resolution = 1, crs = crs(studyArea), extent = ext(studyArea))
  RTM[] <- 1L
  RTM <- terra::mask(RTM, studyArea)
  RTM[sample(1:ncell(RTM), 50)] <- NA

  ## use SA for cropping/masking, not proj
  ## new args
  reproducible::clearCache(userTags = "test", ask = FALSE)
  testthat::expect_warning({
    rawBiomassMap <- prepRawBiomassMap(
      url = biomassURL,
      cropTo = studyArea,
      maskTo = studyArea,
      projectTo = NA
    )
  }, regexp = "CRS do not match") ## prepInputs crs warning

  ## old args
  testthat::expect_warning({
    rawBiomassMap2 <- prepRawBiomassMap(
      url = biomassURL,


      studyArea = studyArea
    )
  }, regexp = "CRS do not match") ## prepInputs crs warning

  testthat::expect_false(crs(studyArea) == crs(rawBiomassMap))
  testthat::expect_true(compareGeom(rawBiomassMap, rawBiomassMap2, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  testthat::expect_false(any(rawBiomassMap[] != rawBiomassMap2[], na.rm = TRUE))

  ## use SA for cropping/masking, proj with RTM
  ## new args
  reproducible::clearCache(userTags = "test", ask = FALSE)
  rawBiomassMap <- prepRawBiomassMap(
    url = biomassURL,
    cropTo = studyArea,
    maskTo = studyArea,
    projectTo = RTM
  )

  ## old args
  rawBiomassMap2 <- prepRawBiomassMap(
    url = biomassURL,
    studyArea = studyArea,
    rasterToMatch = RTM,
    maskWithRTM = FALSE
  )

  testthat::expect_true(st_crs(studyArea) == st_crs(rawBiomassMap))
  testthat::expect_true(compareGeom(rawBiomassMap, rawBiomassMap2, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  testthat::expect_false(any(rawBiomassMap[] != rawBiomassMap2[], na.rm = TRUE))
  testthat::expect_false(all(is.na(rawBiomassMap[]) == is.na(RTM[])))
  testthat::expect_false(all(is.na(rawBiomassMap2[]) == is.na(RTM[])))

  ## use RTM for everything
  ## new args
  # rawBiomassMap <- prepRawBiomassMap(url = biomassURL,
  #                                    to = RTM,
  #                                    projectTo = crs(studyArea))   ## see reproducible #331

  reproducible::clearCache(userTags = "test", ask = FALSE)
  rawBiomassMap <- prepRawBiomassMap(
    url = biomassURL,
    to = RTM
  ) ## for some reason when not interactive the masking doesn't happen if only supplying `to`

  rawBiomassMap <- prepRawBiomassMap(
    url = biomassURL,
    to = RTM
  ) ## for some reason when not interactive the masking doesn't happen if only supplying `to`
  testthat::expect_true(all(is.na(rawBiomassMap[]) == is.na(RTM[])))

  ## old args
  reproducible::clearCache(userTags = "test", ask = FALSE)
  rawBiomassMap2 <- prepRawBiomassMap(
    url = biomassURL,
    studyArea = studyArea,
    rasterToMatch = RTM,
    maskWithRTM = TRUE # ,
    # useSAcrs = TRUE    ## due to reproducible #331 we can't reproduce this.
  )

  testthat::expect_true(compareGeom(rawBiomassMap, rawBiomassMap2, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  testthat::expect_false(any(rawBiomassMap[] != rawBiomassMap2[], na.rm = TRUE))
  testthat::expect_true(all(is.na(rawBiomassMap2[]) == is.na(RTM[]))) ## see reproducible #330

  ## testing w/o URL
  studyTest = {
    targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
                       "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
    ecod <- reproducible::prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip")
    ecod <- ecod[ecod$ECODISTRIC == "332",]
    ecod <- sf::st_transform(ecod, targetCRS)
  }

  ## KNN
  knn2001 <- prepRawBiomassMap(to = studyTest,
                               dataSource = "KNN", dataYear = 2001)
  knn2001_mean <- terra::global(knn2001, mean, na.rm = TRUE)
  knn2011 <- prepRawBiomassMap(to = studyTest,
                               dataSource = "KNN", dataYear = 2011)
  knn2011_mean <- terra::global(knn2011, mean, na.rm = TRUE)
  testthat::expect_true(compareGeom(knn2001, knn2011, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  testthat::expect_true(knn2011_mean < knn2001_mean)

  ## SCANFI
  SCANFI2000 <- prepRawBiomassMap(to = studyTest,
                                  dataSource = "SCANFI", dataYear = 2000)
  SCANFI2000_mean <- terra::global(SCANFI2000, mean, na.rm = TRUE)
  SCANFI2020 <- prepRawBiomassMap(to = studyTest,
                                  dataSource = "SCANFI", dataYear = 2020)
  SCANFI2020_mean <- terra::global(SCANFI2020, mean, na.rm = TRUE)
  testthat::expect_true(compareGeom(SCANFI2000, SCANFI2000, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  testthat::expect_true(SCANFI2020_mean < SCANFI2000_mean)
  testthat::expect_true(SCANFI2000_mean > knn2001_mean)

  ## cache
  knn2001_c1 <- prepRawBiomassMap(to = studyTest, userTags = "cTest",
                                  dataSource = "KNN", dataYear = 2001)
  mess1 <- capture_messages({
    knn2001_c2 <- prepRawBiomassMap(to = studyTest, userTags = "cTest",
                                    dataSource = "KNN", dataYear = 2001)
  })
  expect_true(any(grepl("Loaded! Cached result", mess1)))

  ## cache using file-backed raster
  knn2001_c1 <- prepRawBiomassMap(to = studyTest, userTags = "cTest",
                                  dataSource = "KNN", dataYear = 2001,
                                  writeTo = "testcacheRast.tif")
  mess2 <- capture_messages({
    knn2001_c2 <- prepRawBiomassMap(to = studyTest, userTags = "cTest",
                                    dataSource = "KNN", dataYear = 2001)
  })
  expect_true(any(grepl("Loaded! Cached result", mess2)))
})
