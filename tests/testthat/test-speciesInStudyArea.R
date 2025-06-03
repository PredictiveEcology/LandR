testthat::test_that("speciesInStudyArea works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("googledrive")

  googledrive::drive_deauth()

  td <- withr::local_tempdir("dest_")

  targetCRS <- paste(
    "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0",
    "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  )
  ecod <- reproducible::prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    destinationPath = td, fun = "sf::st_read"
  )
  ecod <- ecod[ecod$ECODISTRIC == "1008", ]
  ecod <- sf::st_transform(ecod, targetCRS)

  speciesInStudy <- speciesInStudyArea(ecod, dPath = td)

  testthat::expect_true("Pinu_Con" %in% speciesInStudy$speciesList)
  testthat::expect_false("Abie_Ama" %in% speciesInStudy$speciesList)
  testthat::expect_false("Pseu_Men_Men" %in% speciesInStudy$speciesList)

  speciesInStudyLandR <- speciesInStudyArea(ecod, dPath = td, sppEquivCol = "LandR")

  testthat::expect_true("Pinu_con" %in% speciesInStudyLandR$speciesList)

  testthat::expect_warning(speciesInStudyArea(ecod, dPath = td, dataSource = "NTEMS"))

  speciesInStudyNTEMS <- suppressWarnings(speciesInStudyArea(ecod, dPath = td, dataSource = "NTEMS"))

  testthat::expect_true("Pinu_con" %in% speciesInStudyNTEMS$speciesList)
  testthat::expect_true("Abie_las" %in% speciesInStudyNTEMS$speciesList)
  testthat::expect_false("Abie_ama" %in% speciesInStudyNTEMS$speciesList)
})
