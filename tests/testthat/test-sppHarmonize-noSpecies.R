## A zero-row sppEquiv is a statement ("no tree species here"), not a missing table. It used
## to be replaced by all of LandR::sppEquivalencies_CA, which turned a non-forested study area
## into a 195-species one downstream (fireSense_dataPrepFit, 2026-09-12).
test_that("sppHarmonize keeps a supplied zero-row sppEquiv", {
  noSpp <- LandR::sppEquivalencies_CA[0]
  out <- suppressMessages(sppHarmonize(sppEquiv = noSpp, sppNameVector = NULL,
                                       sppEquivCol = "LandR", sppColorVect = NULL))
  expect_identical(NROW(out$sppEquiv), 0L)
  expect_identical(out$sppNameVector, character(0))
  expect_identical(out$sppEquivCol, "LandR")
})

test_that("sppHarmonize still falls back to the full table when sppEquiv is NULL", {
  out <- suppressMessages(sppHarmonize(sppEquiv = NULL, sppNameVector = c("Pice_mar", "Pinu_ban"),
                                       sppEquivCol = "LandR", sppColorVect = NULL))
  expect_setequal(out$sppEquiv$LandR, c("Pice_mar", "Pinu_ban"))
})
