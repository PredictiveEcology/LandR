test_that("each LandR species has a single FuelClass in sppEquivalencies_CA", {
  ## Several SCANFI varieties share one LandR name (e.g. PSEU_MEN, PSEU_MEN_GLA and
  ## PSEU_MEN_MEN are all `Pseu_men`). Joins keyed on the LandR column see every one of
  ## those rows, so if they disagree on FuelClass the species is assigned two fuel classes
  ## and fireSenseUtils::cohortsToFuelClasses() stops. Coastal Douglas-fir once carried
  ## "CedrMplOther" while the other two carried "DgFrPoPine".
  sppEquiv <- data.table::as.data.table(LandR::sppEquivalencies_CA)
  withFuel <- unique(sppEquiv[nzchar(LandR) & nzchar(FuelClass), c("LandR", "FuelClass")])
  multi <- withFuel[, list(classes = paste(sort(FuelClass), collapse = ", ")), by = "LandR"][
    grepl(",", classes, fixed = TRUE)]

  expect_identical(nrow(multi), 0L,
                   info = paste(multi$LandR, "->", multi$classes, collapse = "; "))
})
