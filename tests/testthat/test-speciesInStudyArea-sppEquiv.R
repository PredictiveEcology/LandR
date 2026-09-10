## `speciesInStudyArea()` returns the species table (`sppEquiv`) a study area needs, built the
## way fireSense_ELFs used to build it by hand: no `_Spp` genus entries, only species with
## LANDIS traits, and Engelmann spruce's two SCANFI entries merged into one `Pice_eng` row.
## These use a small species-presence raster, so nothing is downloaded.

speciesPresence <- function(categories) {
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4000, ymin = 0, ymax = 4000,
                   crs = "EPSG:3978")
  terra::values(r) <- rep(seq_along(categories), length.out = terra::ncell(r))
  levels(r) <- data.frame(ID = seq_along(categories), category = categories)
  r
}

studyAreaFor <- function(r) terra::as.polygons(terra::ext(r), crs = terra::crs(r))

test_that("a supplied speciesPresentRas is used (no 'bb' not found)", {
  r <- speciesPresence(c("ABIE_AMA__PSEU_MEN", "THUJ_PLI__TSUG_HET"))
  out <- speciesInStudyArea(studyAreaFor(r), speciesPresentRas = r)
  expect_setequal(out$speciesList, c("ABIE_AMA", "PSEU_MEN", "THUJ_PLI", "TSUG_HET"))
})

test_that("sppEquiv: no Engelmann spruce gives a filtered table, not NULL", {
  withr::local_package("data.table")
  r <- speciesPresence(c("ABIE_AMA__PSEU_MEN__THUJ_PLI", "POPU_GRA__TSUG_HET"))
  out <- speciesInStudyArea(studyAreaFor(r), speciesPresentRas = r)

  expect_s3_class(out$sppEquiv, "data.table")
  expect_true(all(c("Abie_ama", "Pseu_men", "Thuj_pli", "Tsug_het") %in% out$sppEquiv$LandR))
  ## species without LANDIS traits are left out
  expect_true(all(out$sppEquiv$LANDIS_traits != ""))
  expect_false("Popu_gra" %in% out$sppEquiv$LandR)
})

test_that("sppEquiv: Engelmann spruce's two entries become one Pice_eng", {
  withr::local_package("data.table")
  r <- speciesPresence(c("PICE_ENG__PSEU_MEN", "PICE_ENG_GLA"))
  out <- speciesInStudyArea(studyAreaFor(r), speciesPresentRas = r)

  expect_true("Pice_eng" %in% out$sppEquiv$LandR)
  expect_false("Pice_eng_gla" %in% out$sppEquiv$LandR)
  expect_identical(anyDuplicated(out$sppEquiv), 0L)
})

test_that("sppEquiv is keyed on sppEquivCol when one is given", {
  withr::local_package("data.table")
  r <- speciesPresence(c("ABIE_AMA__PSEU_MEN"))
  out <- speciesInStudyArea(studyAreaFor(r), speciesPresentRas = r, sppEquivCol = "LandR")
  expect_true(all(c("Abie_ama", "Pseu_men") %in% out$sppEquiv$LandR))
})
