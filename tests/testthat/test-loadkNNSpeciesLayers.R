## TEMPORARY to deal with intermittent NFI server SSL issue
.sslVerify <- local({
  tryURL <- paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                   "canada-forests-attributes_attributs-forests-canada/2001",
                   "-attributes_attributs-2001/")
  sslIssue <- isTRUE(tryCatch({
    suppressWarnings(RCurl::getURL(tryURL)) ## partial match warning in Rcurl
  }, error = function(e) {
    grepl("SSL certificate problem", paste(e))
  }))

  ifelse(sslIssue, 0L, unname(curl::curl_options("^ssl_verifypeer$")))
})

test_that("test download kNN SpeciesLayers with kNN website - all species", {
  skip_if_not(interactive())
  skip_if_not_installed(c("curl", "googledrive", "httr", "RCurl", "XML"))

  mainDir <- tempdir()
  dPath <- file.path(mainDir, "inputs")
  cPath <- file.path(mainDir, "cache")
  set.seed(123)
  SA <- randomStudyArea(size = 10000000)
  RTM <- rast(SA, res = 250)
  sppEquiv <- sppEquivalencies_CA
  sppEquivCol <- "KNN"

  sppEquiv <- sppEquiv[KNN != ""]

  ## get all available species for 2001
  opt_cachePath <- getOption("reproducible.cachePath")
  options(reproducible.cachePath = cPath)

  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers <- loadkNNSpeciesLayers(
      dPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2001,
      url = paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                   "canada-forests-attributes_attributs-forests-canada/",
                   "2001-attributes_attributs-2001/"),
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })

  expectedSpp <- c("Abie_Las", "Betu_Pap", "Pice_Eng",
                   "Pice_Gla", "Pice_Mar", "Pinu_Alb",
                   "Pinu_Con", "Pinu_Spp", "Popu_Bal",
                   "Popu_Tre", "Pseu_Men", "Broadleaf_Spp",
                   "Needleleaf_Spp")

  expect_true(all(expectedSpp %in% names(speciesLayers)))
  expect_true(.compareRas(RTM, speciesLayers, res = TRUE, stopiffalse = FALSE))

  ## get all available species for 2011
  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers2011 <- loadkNNSpeciesLayers(
      dPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2011,
      url = paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                   "canada-forests-attributes_attributs-forests-canada/",
                   "2011-attributes_attributs-2011/"),
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })
  expectedSpp <- c("Abie_Las", "Abie_Spp", "Betu_Pap",
                   "Lari_Occ", "Pice_Eng", "Pice_Gla",
                   "Pice_Mar", "Pice_Spp", "Pinu_Con",
                   "Pinu_Spp", "Popu_Bal", "Popu_Tre",
                   "Pseu_Men", "Broadleaf_Spp", "Needleleaf_Spp")
  expect_true(all(expectedSpp %in% names(speciesLayers2011)))
  expect_true(.compareRas(RTM, speciesLayers2011, res = TRUE, stopOnError = FALSE))

  options(reproducible.cachePath = opt_cachePath)
  unlink(mainDir, recursive = TRUE)
})

test_that("test download kNN SpeciesLayers with kNN website - three species", {
  skip_if_not(interactive())
  skip_if_not_installed(c("googledrive", "RCurl", "XML"))

  mainDir <- tempdir()
  dPath <- file.path(mainDir, "inputs")
  cPath <- file.path(mainDir, "cache")
  set.seed(123)
  SA <- randomStudyArea(size = 10000000)
  RTM <- rast(SA, res = 250)
  sppEquiv <- sppEquivalencies_CA
  sppEquivCol <- "KNN"

  sppEquiv <- sppEquiv[KNN %in% c("Pinu_Spp", "Popu_Tre", "Pice_Mar")]

  ## get all available species for 2001
  opt_cachePath <- getOption("reproducible.cachePath")
  options(reproducible.cachePath = cPath)

  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers <- loadkNNSpeciesLayers(
      dPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2001,
      url = paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                   "canada-forests-attributes_attributs-forests-canada/",
                   "2001-attributes_attributs-2001/"),
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })
  expect_true(all(sppEquiv$KNN %in% names(speciesLayers)))
  expect_true(.compareRas(RTM, speciesLayers, res = TRUE, stopOnError = FALSE))

  ## get all available species for 2011
  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers2011 <- loadkNNSpeciesLayers(
      dPath = dPath,
      rasterToMatch = RTM,
      studyArea = SA,
      year = 2011,
      url = paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                   "canada-forests-attributes_attributs-forests-canada/",
                   "2011-attributes_attributs-2011/"),
      sppEquiv = sppEquiv,
      sppEquivCol = sppEquivCol,
      thresh = 0
    )
  })
  expect_true(all(sppEquiv$KNN %in% names(speciesLayers2011)))
  expect_true(.compareRas(RTM, speciesLayers2011, res = TRUE, stopOnError = FALSE))

  options(reproducible.cachePath = opt_cachePath)
  unlink(mainDir, recursive = TRUE)
})

test_that("test download kNN SpeciesLayers bad website - three species", {
  skip_if_not(interactive())
  skip_if_not_installed(c("googledrive", "RCurl", "XML"))

  mainDir <- tempdir()
  dPath <- file.path(mainDir, "inputs")
  cPath <- file.path(mainDir, "cache")
  set.seed(123)
  SA <- randomStudyArea(size = 10000000)
  RTM <- rast(SA, res = 250)
  sppEquiv <- sppEquivalencies_CA
  sppEquivCol <- "KNN"

  sppEquiv <- sppEquiv[KNN %in% c("Pinu_Spp", "Popu_Tre", "Pice_Mar")]

  ## get all available species for 2001
  opt_cachePath <- getOption("reproducible.cachePath")
  options(reproducible.cachePath = cPath)

  httr::with_config(config = httr::config(ssl_verifypeer = .sslVerify), {
    speciesLayers <- prepSpeciesLayers_KNN(destinationPath = dPath,
                                           outputPath = dPath,
                                           rasterToMatch = RTM,
                                           studyArea = SA,
                                           year = 2001,
                                           url = "dvnsebvebvwebv.cajey/aebcbeh/",
                                           sppEquiv = sppEquiv,
                                           sppEquivCol = sppEquivCol,
                                           thresh = 0)
  })
  expect_true(all(sppEquiv$KNN %in% names(speciesLayers)))
  expect_true(LandR::.compareRas(RTM, speciesLayers, res = TRUE, stopOnError = FALSE))

  ## get all available species for 2011
  speciesLayers2011 <- prepSpeciesLayers_KNN(destinationPath = dPath,
                                             outputPath = dPath,
                                             rasterToMatch = RTM,
                                             studyArea = SA,
                                             year = 2011,
                                             url = "dvnsebvebvwebv.cajey/aebcbeh/",
                                             sppEquiv = sppEquiv,
                                             sppEquivCol = sppEquivCol,
                                             thresh = 0)
  expect_true(all(sppEquiv$KNN %in% names(speciesLayers2011)))
  expect_true(.compareRas(RTM, speciesLayers2011, res = TRUE, stopOnError = FALSE))

  options(reproducible.cachePath = opt_cachePath)
  unlink(mainDir, recursive = TRUE)
})
