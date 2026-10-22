## Several species names can map to one species in sppEquivCol (PSEU_MEN and PSEU_MEN_GLA are both
## Pseu_men in "LandR"). A module that calls sppHarmonize() with the output of an earlier module's
## call (Biomass_speciesData, then Biomass_borealDataPrep) must get the same species back.

test_that("sppHarmonize() accepts its own output when several names merge into one species", {
  nms <- c("PSEU_MEN", "PSEU_MEN_GLA", "POPU_TRE")
  o1 <- suppressWarnings(sppHarmonize(NULL, nms, "LandR", NULL, vegLeadingProportion = 0.8))
  expect_setequal(names(o1$sppColorVect), c("Pseu_men", "Popu_tre", "Mixed"))

  o2 <- suppressWarnings(sppHarmonize(o1$sppEquiv, o1$sppNameVector, "LandR", o1$sppColorVect,
                                      vegLeadingProportion = 0.8))
  expect_identical(o2$sppColorVect, o1$sppColorVect)
  expect_identical(o2$sppEquiv, o1$sppEquiv)
})
