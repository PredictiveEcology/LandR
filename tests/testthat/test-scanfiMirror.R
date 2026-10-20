## SCANFI Drive ids that 404 for anonymous users (2020 land cover, 2020 stand age) are served
## from the arbutus mirror through reproducible's urlRemap hook. No network here: these check
## the mapping, not the download.

test_that("the SCANFI remap sends the 404ing Drive ids to arbutus", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  remap <- scanfiUrlRemap()
  expect_true(is.function(remap))

  lcc2020 <- remap("https://drive.google.com/file/d/1EGp7LUA7cXMR6KpXDmu617xsjwGM6aIx",
                   filename = NA_character_)
  expect_match(lcc2020, "^https://object-arbutus\\.alliancecan\\.ca/")
  expect_match(lcc2020, "CanadaLCCclassCodes_2020_v2_20260119\\.tif$")

  age2020 <- remap("https://drive.google.com/file/d/1nXPS3bpFUESYieNfXO25OKlZJEgqtRnD",
                   filename = NA_character_)
  expect_match(age2020, "SCANFI_age_median_2020_v2_20260119\\.tif$")
})

test_that("the SCANFI species folders are remapped too, so listing needs no login", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  dirs <- attr(scanfiUrlRemap(), "byDir")
  ## the 2020 V2 species folder used by prepSpeciesLayers_SCANFI()
  expect_true("15T4HIFeqzwp0TuOuxmYoexuXdLFnCZBi" %in% names(dirs))
  expect_length(dirs, 9L)   # one per SCANFI v2 year, 1985-2025
})

test_that("an unrelated URL is left alone", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  expect_null(scanfiUrlRemap()("https://example.org/some/file.tif", filename = "file.tif"))
})

test_that("every manifest url is on the current arbutus host", {
  m <- utils::read.csv(system.file("extdata", "arbutus_manifest_SCANFI_v2.csv", package = "LandR"),
                       stringsAsFactors = FALSE)
  expect_true(all(startsWith(m$url, "https://object-arbutus.alliancecan.ca/")))
  expect_false(any(grepl("cloud.computecanada.ca", m$url, fixed = TRUE)))
})

test_that("LandR does not set reproducible.urlRemap when it loads", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  expect_false(exists(".setScanfiMirror", envir = asNamespace("LandR"), inherits = FALSE))
})

test_that("the default remap is applied for the call only, and the global option is untouched", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  withr::local_options(list(reproducible.urlRemap = NULL, LandR.scanfiMirror = TRUE))
  inside <- .withScanfiRemap(getOption("reproducible.urlRemap"))
  expect_true(is.function(inside))
  expect_null(getOption("reproducible.urlRemap"))

  ## also when the call fails
  try(.withScanfiRemap(stop("boom")), silent = TRUE)
  expect_null(getOption("reproducible.urlRemap"))
})

test_that("a remap the user set always wins", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  mine <- function(url, filename) NULL
  withr::local_options(list(reproducible.urlRemap = mine, LandR.scanfiMirror = TRUE))
  expect_identical(.withScanfiRemap(getOption("reproducible.urlRemap")), mine)
  expect_identical(getOption("reproducible.urlRemap"), mine)
})

test_that("the default remap can be switched off", {
  withr::local_options(list(reproducible.urlRemap = NULL, LandR.scanfiMirror = FALSE))
  expect_null(.withScanfiRemap(getOption("reproducible.urlRemap")))
})

test_that("SCANFI downloads run under the default remap, and leave no global option", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  withr::local_options(list(reproducible.urlRemap = NULL, LandR.scanfiMirror = TRUE,
                            reproducible.cachePath = withr::local_tempdir(),
                            reproducible.useMemoise = FALSE))
  seen <- NULL
  testthat::local_mocked_bindings(
    prepInputs = function(...) {
      seen <<- getOption("reproducible.urlRemap")
      terra::rast(nrows = 2, ncols = 2, vals = 50)
    }
  )
  suppressMessages(prepInputsStandAgeMap(dataSource = "SCANFI", dataYear = 2020, dataVersion = "V2",
                                         destinationPath = withr::local_tempdir(), fireURL = NA))
  expect_true(is.function(seen))
  expect_null(getOption("reproducible.urlRemap"))
})

test_that("the forest-land land cover reads run under the default remap too", {
  ## prepInputs_SCANFI_LCC_FAO() reads the land cover of each forest-land year besides `year`'s
  ## own. Those reads sat outside the per-download remap, so they went to Google Drive.
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  withr::local_options(list(reproducible.urlRemap = NULL, LandR.scanfiMirror = TRUE,
                            LandR.forestLandUseCache = NULL, reproducible.useCache = FALSE,
                            reproducible.cachePath = withr::local_tempdir(),
                            reproducible.useMemoise = FALSE))
  lcc <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000,
                     crs = "EPSG:3978", vals = c(210, 50, 20, 210))
  remapped <- character()
  testthat::local_mocked_bindings(
    prepInputs = function(url = NULL, ...) {
      if (is.function(getOption("reproducible.urlRemap"))) remapped <<- c(remapped, url)
      lcc
    }
  )
  suppressMessages(prepInputs_SCANFI_LCC_FAO(year = 2020, dataVersion = "V2",
                                             forestLandFrom = "lccYears",
                                             destinationPath = withr::local_tempdir()))
  forestLandYears <- setdiff(.defaultForestLandYears(.scanfi_v2_years), 2020)
  expect_length(forestLandYears, 5)
  for (y in forestLandYears) {
    expect_true(.scanfiLCCFAOSource(y, "V2")$url %in% remapped, info = paste("year", y))
  }
  expect_null(getOption("reproducible.urlRemap"))
})

test_that("a supplied ageURL reaches prepInputs for SCANFI V1 and V2", {
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useMemoise = FALSE)
  seen <- character()
  testthat::local_mocked_bindings(
    prepInputs = function(..., url = NULL) {
      seen <<- c(seen, url)
      terra::rast(nrows = 2, ncols = 2, vals = 50)
    }
  )
  dPath <- withr::local_tempdir()
  mine <- c(V1 = "https://example.invalid/my_age_v1.tif", V2 = "https://example.invalid/my_age_v2.tif")
  for (v in c("V1", "V2")) {
    suppressMessages(prepInputsStandAgeMap(dataSource = "SCANFI", dataYear = 2020, dataVersion = v,
                                           ageURL = mine[[v]], destinationPath = dPath, fireURL = NA))
  }
  expect_equal(seen, unname(mine))
})

test_that("access errors name the mirror when it failed, Google Drive when Drive did", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  withr::local_options(list(reproducible.urlRemap = NULL, LandR.scanfiMirror = TRUE))
  msg <- function(url, version = "V2", err = "Client error: (404) Not Found") {
    out <- tryCatch(
      LandR:::.withSCANFIAccess(stop(simpleError(err)), what = "the SCANFI stand age map", dataYear = 2020,
                                dataVersion = version, urlArg = "ageURL", url = url),
      error = function(e) conditionMessage(e)
    )
    gsub("[[:space:]]+", " ", out)
  }
  v2id <- "https://drive.google.com/file/d/1nXPS3bpFUESYieNfXO25OKlZJEgqtRnD" ## 2020 age
  mirror <- msg(v2id)
  expect_match(mirror, "from the SCANFI mirror", fixed = TRUE)
  expect_false(grepl("Google Drive", mirror, fixed = TRUE))
  expect_false(grepl("drive_auth", mirror, fixed = TRUE))

  for (other in list(msg(v2id, "V1"), msg("https://drive.google.com/file/d/abc"),
                     ## the mirror was meant to apply, but the error says Drive refused
                     msg(v2id, err = paste0("Could not access the Google Drive resource:\n  ", v2id)))) {
    expect_match(other, "from Google Drive", fixed = TRUE)
    expect_false(grepl("SCANFI mirror.", other, fixed = TRUE))
  }
  withr::local_options(list(reproducible.urlRemap = function(url, filename) NULL))
  expect_match(msg(v2id), "from Google Drive", fixed = TRUE)
})
