## SCANFI Drive ids that 404 for anonymous users (2020 land cover, 2020 stand age) are served
## from the arbutus mirror through reproducible's urlRemap hook. No network here: these check
## the mapping, not the download.

test_that("the SCANFI remap sends the 404ing Drive ids to arbutus", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  remap <- scanfiUrlRemap()
  expect_true(is.function(remap))

  lcc2020 <- remap("https://drive.google.com/file/d/1EGp7LUA7cXMR6KpXDmu617xsjwGM6aIx",
                   filename = NA_character_)
  expect_match(lcc2020, "^https://object-arbutus\\.cloud\\.computecanada\\.ca/")
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

test_that("loading LandR never replaces a remap the user set", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  mine <- function(url, filename) NULL
  withr::local_options(list(reproducible.urlRemap = mine, LandR.scanfiMirror = TRUE))
  expect_false(.setScanfiMirror())
  expect_identical(getOption("reproducible.urlRemap"), mine)
})

test_that("the mirror is installed when nothing is set, and can be switched off", {
  skip_if_not("makeUrlRemap" %in% getNamespaceExports("reproducible"))
  withr::local_options(list(reproducible.urlRemap = NULL, LandR.scanfiMirror = FALSE))
  expect_false(.setScanfiMirror())
  expect_null(getOption("reproducible.urlRemap"))

  options(LandR.scanfiMirror = TRUE)
  expect_true(.setScanfiMirror())
  expect_true(is.function(getOption("reproducible.urlRemap")))
})
