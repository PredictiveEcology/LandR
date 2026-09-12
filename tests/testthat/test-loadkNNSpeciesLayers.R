## TEMPORARY to deal with intermittent NFI server SSL issue
.sslVerify <- local({
  tryURL <- paste0(
    "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
    "canada-forests-attributes_attributs-forests-canada/2001",
    "-attributes_attributs-2001/"
  )
  sslIssue <- isTRUE(tryCatch(
    {
      suppressWarnings(RCurl::getURL(tryURL)) ## partial match warning in Rcurl
    },
    error = function(e) {
      grepl("SSL certificate problem", paste(e))
    }
  ))

  ifelse(sslIssue, 0L, unname(curl::curl_options("^ssl_verifypeer$")))
})

testthat::test_that("test download kNN SpeciesLayers with kNN website - all species", {
  testthat::skip_on_cran()
  testthat::skip_if_not(interactive())
  testthat::skip_if_not_installed(c("curl", "googledrive", "httr", "RCurl", "withr", "XML"))

  cPath <- withr::local_tempdir("cache_")
  dPath <- withr::local_tempdir("inputs_")

  withr::local_options(list(reproducible.cachePath = cPath))
  withr::local_seed(123)

  SA <- randomStudyArea(size = 10000000)
  RTM <- terra::rast(SA, resolution = 250, vals = 1L)
  sppEquiv <- sppEquivalencies_CA
  sppEquivCol <- "KNN"

  sppEquiv <- sppEquiv[KNN != ""]

  ## get all available species for 2001
  testthat::expect_warning({
    httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
      speciesLayers <- loadkNNSpeciesLayers(
        dPath = dPath,
        rasterToMatch = RTM,
        studyArea = SA,
        year = 2001,
        url = paste0(
          "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
          "canada-forests-attributes_attributs-forests-canada/",
          "2001-attributes_attributs-2001/"
        ),
        sppEquiv = sppEquiv,
        sppEquivCol = sppEquivCol,
        thresh = 0
      )
      })
  }, regexp = "Will use remaining matching species, but check if this is correct")
  ## TODO: WARNING: [readValues] raster has no values

  expectedSpp <- c(
    "Abie_Las", "Betu_Pap", "Pice_Eng",
    "Pice_Gla", "Pice_Mar", "Pinu_Alb",
    "Pinu_Con", "Pinu_Spp", "Popu_Bal",
    "Popu_Tre", "Pseu_Men", "Broadleaf_Spp",
    "Needleleaf_Spp"
  )

  testthat::expect_true(all(expectedSpp %in% names(speciesLayers)))
  testthat::expect_true(.compareRas(RTM, speciesLayers, res = TRUE, stopiffalse = FALSE))

  ## get all available species for 2011
  testthat::expect_warning({
    httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
      speciesLayers2011 <- loadkNNSpeciesLayers(
        dPath = dPath,
        rasterToMatch = RTM,
        studyArea = SA,
        year = 2011,
        url = paste0(
          "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
          "canada-forests-attributes_attributs-forests-canada/",
          "2011-attributes_attributs-2011/"
        ),
        sppEquiv = sppEquiv,
        sppEquivCol = sppEquivCol,
        thresh = 0
      )
    })
  }, regexp = "Will use remaining matching species, but check if this is correct")
  ## TODO: WARNING: [readValues] raster has no values

  expectedSpp <- c(
    "Abie_Las", "Abie_Spp", "Betu_Pap",
    "Lari_Occ", "Pice_Eng", "Pice_Gla",
    "Pice_Mar", "Pice_Spp", "Pinu_Con",
    "Pinu_Spp", "Popu_Bal", "Popu_Tre",
    "Pseu_Men", "Broadleaf_Spp", "Needleleaf_Spp"
  )
  testthat::expect_true(all(expectedSpp %in% names(speciesLayers2011)))
  testthat::expect_true(.compareRas(RTM, speciesLayers2011, res = TRUE, stopOnError = FALSE)) ## TODO: fails
})

testthat::test_that("test download kNN SpeciesLayers with kNN website - three species", {
  testthat::skip_on_cran()
  testthat::skip_if_not(interactive())
  testthat::skip_if_not_installed(c("googledrive", "RCurl", "withr", "XML"))
  testthat::skip_if_offline()

  cPath <- withr::local_tempdir("cache_")
  dPath <- withr::local_tempdir("inputs_")

  withr::local_options(reproducible.cachePath = cPath)

  withr::local_seed(123)

  SA <- randomStudyArea(size = 10000000)
  RTM <- terra::rast(SA, resolution = 250, vals = 1L)
  sppEquiv <- sppEquivalencies_CA
  sppEquivCol <- "KNN"

  sppEquiv <- sppEquiv[KNN %in% c("Pinu_Spp", "Popu_Tre", "Pice_Mar")]

  ## get all available species for 2001
  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers <- loadkNNSpeciesLayers(
      dPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2001,
      url = paste0(
        "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
        "canada-forests-attributes_attributs-forests-canada/",
        "2001-attributes_attributs-2001/"
      ),
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })
  ## TODO: WARNING: [readValues] raster has no values

  testthat::expect_true(all(sppEquiv$KNN %in% names(speciesLayers)))
  testthat::expect_true(.compareRas(RTM, speciesLayers, res = TRUE, stopOnError = FALSE))

  ## get all available species for 2011
  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers2011 <- loadkNNSpeciesLayers(
      dPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2011,
      url = paste0(
        "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
        "canada-forests-attributes_attributs-forests-canada/",
        "2011-attributes_attributs-2011/"
      ),
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })
  ## TODO: WARNING: [readValues] raster has no values

  testthat::expect_true(all(sppEquiv$KNN %in% names(speciesLayers2011)))
  testthat::expect_true(.compareRas(RTM, speciesLayers2011, res = TRUE, stopOnError = FALSE))
})

testthat::test_that("test download kNN SpeciesLayers bad website - three species", {
  testthat::skip_if_not(interactive())
  testthat::skip_if_not_installed(c("googledrive", "RCurl", "withr", "XML"))
  testthat::skip_if_offline()

  cPath <- withr::local_tempdir("cache_")
  dPath <- withr::local_tempdir("inputs_")

  withr::local_seed(123)

  SA <- randomStudyArea(size = 10000000)
  RTM <- terra::rast(SA, resolution = 250, vals = 1L)
  sppEquiv <- sppEquivalencies_CA
  sppEquivCol <- "KNN"

  sppEquiv <- sppEquiv[KNN %in% c("Pinu_Spp", "Popu_Tre", "Pice_Mar")]

  ## get all available species for 2001
  withr::local_options(list(reproducible.cachePath = cPath))

  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers <- prepSpeciesLayers_KNN(
      destinationPath = dPath,
      outputPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2001,
      url = "dvnsebvebvwebv.cajey/aebcbeh/", ## triggers drive_auth(); needs interactive
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })
  testthat::expect_true(all(sppEquiv$KNN %in% names(speciesLayers)))
  testthat::expect_true(LandR::.compareRas(RTM, speciesLayers, res = TRUE, stopOnError = FALSE))

  ## get all available species for 2011
  speciesLayers2011 <- prepSpeciesLayers_KNN(
    destinationPath = dPath,
    outputPath = dPath,
    rasterToMatch = RTM,
    studyArea = SA,
    year = 2011,
    url = "dvnsebvebvwebv.cajey/aebcbeh/", ## triggers drive_auth(); needs interactive
    sppEquiv = sppEquiv,
    sppEquivCol = sppEquivCol,
    thresh = 0
  )
  testthat::expect_true(all(sppEquiv$KNN %in% names(speciesLayers2011)))
  testthat::expect_true(.compareRas(RTM, speciesLayers2011, res = TRUE, stopOnError = FALSE))
})
