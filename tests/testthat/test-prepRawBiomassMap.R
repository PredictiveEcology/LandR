testthat::test_that("test prepRawBiomassMap", {
  testthat::skip_on_cran() ## not necessary, but here in case LandR goes to CRAN
  testthat::skip_on_ci()

  withr::local_package("reproducible")
  withr::local_package("SpaDES.tools")
  withr::local_package("terra")
  withr::local_package("sf")

  dPath <- file.path(tempdir(), "inputs")
  opts <- options("reproducible.inputPaths" = NULL,
                  "reproducible.overwrite" = TRUE,
                  "reproducible.useTerra" = TRUE,
                  "reproducible.rasterRead" = "terra::rast",
                  "reproducible.destinationPath" = dPath)

  on.exit(options(opts), add = TRUE)

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
  rawBiomassMap <- suppressWarnings({
    prepRawBiomassMap(url = biomassURL,
                      studyAreaName = "test",
                      cacheTags = "test",
                      cropTo = studyArea,
                      maskTo = studyArea,
                      projectTo = NA)
  })

  ## old args
  rawBiomassMap2 <- suppressWarnings({
    prepRawBiomassMap(url = biomassURL,
                      studyAreaName = "test",
                      cacheTags = "test",
                      studyArea = studyArea)
  })

  expect_false(crs(studyArea) == crs(rawBiomassMap))
  expect_true(compareGeom(rawBiomassMap, rawBiomassMap2, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  expect_false(any(rawBiomassMap[] != rawBiomassMap2[], na.rm = TRUE))

  ## use SA for cropping/masking, proj with RTM
  ## new args
  reproducible::clearCache(userTags = "test", ask = FALSE)
  rawBiomassMap <- prepRawBiomassMap(url = biomassURL,
                                     studyAreaName = "test",
                                     cacheTags = "test",
                                     cropTo = studyArea,
                                     maskTo = studyArea,
                                     projectTo = RTM)

  ## old args
  rawBiomassMap2 <- prepRawBiomassMap(url = biomassURL,
                                      studyAreaName = "test",
                                      cacheTags = "test",
                                      studyArea = studyArea,
                                      rasterToMatch = RTM,
                                      maskWithRTM = FALSE)

  expect_true(st_crs(studyArea) == st_crs(rawBiomassMap))
  expect_true(compareGeom(rawBiomassMap, rawBiomassMap2, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  expect_false(any(rawBiomassMap[] != rawBiomassMap2[], na.rm = TRUE))
  expect_false(all(is.na(rawBiomassMap[]) == is.na(RTM[])))
  expect_false(all(is.na(rawBiomassMap2[]) == is.na(RTM[])))

  ## use RTM for everything
  ## new args
  # rawBiomassMap <- prepRawBiomassMap(url = biomassURL,
  #                                    studyAreaName = "test",
  #                                    cacheTags = "test",
  #                                    to = RTM,
  #                                    projectTo = crs(studyArea))   ## this is failing; reported issue #331-reproducible

  reproducible::clearCache(userTags = "test", ask = FALSE)
  rawBiomassMap <- prepRawBiomassMap(url = biomassURL,
                                     studyAreaName = "test",
                                     cacheTags = "test",
                                     to = RTM)   ## for some reason when not interactive the masking doesn't happen if only supplying `to`

  rawBiomassMap <- prepRawBiomassMap(url = biomassURL,
                                     studyAreaName = "test",
                                     cacheTags = "test",
                                     to = RTM)   ## for some reason when not interactive the masking doesn't happen if only supplying `to`
  expect_true(all(is.na(rawBiomassMap[]) == is.na(RTM[])))

  ## old args
  reproducible::clearCache(userTags = "test", ask = FALSE)
  rawBiomassMap2 <- prepRawBiomassMap(url = biomassURL,
                                      studyAreaName = "test",
                                      cacheTags = "test",
                                      studyArea = studyArea,
                                      rasterToMatch = RTM,
                                      maskWithRTM = TRUE #,
                                      # useSAcrs = TRUE    ## due to issue #331-reproducible we can't reproduce this.
                                      )

  expect_true(compareGeom(rawBiomassMap, rawBiomassMap2, rowcol = TRUE, res = TRUE, stopOnError = FALSE))
  expect_false(any(rawBiomassMap[] != rawBiomassMap2[], na.rm = TRUE))
  # expect_true(all(is.na(rawBiomassMap2[]) == is.na(RTM[])))    ### May 25 2023, reported as issue #330 on reproducible
})
