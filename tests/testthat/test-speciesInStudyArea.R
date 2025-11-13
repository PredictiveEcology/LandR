testthat::test_that("speciesInStudyArea works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("googledrive")

  # googledrive::drive_deauth()

  td <- withr::local_tempdir("dest_")

  targetCRS <- paste(
    "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77",
    "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  )
  ecod <- reproducible::prepInputs(
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    destinationPath = td,
    fun = "sf::st_read"
  )
  ecod <- ecod[ecod$ECODISTRIC == "1008", ]
  ecod <- sf::st_transform(ecod, targetCRS)

  speciesInStudy <- speciesInStudyArea(ecod, dataSource = "KNN", dPath = td)

  testthat::expect_true("Pinu_Con" %in% speciesInStudy$speciesList)
  testthat::expect_false("Abie_Ama" %in% speciesInStudy$speciesList)
  testthat::expect_false("Pseu_Men_Men" %in% speciesInStudy$speciesList)

  speciesInStudyLandR <- speciesInStudyArea(ecod, dPath = td, sppEquivCol = "LandR")

  testthat::expect_true("Pinu_con" %in% speciesInStudyLandR$speciesList)

  testthat::expect_warning(speciesInStudyArea(ecod, dPath = td, dataSource = "NTEMS"))

  speciesInStudyNTEMS <- suppressWarnings(speciesInStudyArea(
    ecod,
    dPath = td,
    dataSource = "NTEMS"
  ))

  testthat::expect_true("Pinu_con" %in% speciesInStudyNTEMS$speciesList)
  testthat::expect_true("Abie_las" %in% speciesInStudyNTEMS$speciesList)
  testthat::expect_false("Abie_ama" %in% speciesInStudyNTEMS$speciesList)

  speciesInStudySCANFI <- speciesInStudyArea(ecod, dPath = td, dataSource = "SCANFI")

  testthat::expect_true("PINU_CON_LAT" %in% speciesInStudySCANFI$speciesList)
  testthat::expect_false("ABIE_AMA" %in% speciesInStudySCANFI$speciesList)
  testthat::expect_true("PSEU_MEN" %in% speciesInStudySCANFI$speciesList)
  testthat::expect_false("PINU_STR" %in% speciesInStudySCANFI$speciesList)
})
