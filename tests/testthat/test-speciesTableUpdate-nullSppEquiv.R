test_that("speciesTableUpdate() uses sppEquivalencies_CA when sppEquiv is NULL", {
  ## `data.table(utils::data("sppEquivalencies_CA", ...))` built a one-column table holding
  ## the *name* of the dataset, not the dataset, so the NULL default could never match a
  ## species and the trait updates below were silently or noisily lost.
  traitCols <- LandR:::.speciesTableColNames

  species <- data.table::data.table(
    species = "Abie_bal", Area = "BSW", longevity = 150L, sexualmature = 25L,
    shadetolerance = 5, firetolerance = 1L, seeddistance_eff = 30L,
    seeddistance_max = 160L, resproutprob = 0, resproutage_min = 0L,
    resproutage_max = 0L, postfireregen = "none", leaflongevity = 3L,
    wooddecayrate = 0.1, mortalityshape = 15L, growthcurve = 0, leafLignin = 0.2,
    hardsoft = "soft"
  )
  speciesTable <- data.table::copy(species)[, species := "ABIE.BAL"]
  data.table::setcolorder(species, traitCols)
  data.table::setcolorder(speciesTable, traitCols)

  updated <- speciesTableUpdate(
    species = data.table::copy(species), speciesTable = data.table::copy(speciesTable),
    sppEquiv = NULL, sppEquivCol = "LandR"
  )

  ## Burton & Cumming (1995) values for balsam fir, applied by speciesTableUpdate()
  expect_identical(updated$longevity, 200L)
  expect_identical(updated$shadetolerance, 3)

  ## the same as passing the table explicitly
  explicit <- speciesTableUpdate(
    species = data.table::copy(species), speciesTable = data.table::copy(speciesTable),
    sppEquiv = data.table::as.data.table(LandR::sppEquivalencies_CA), sppEquivCol = "LandR"
  )
  expect_identical(updated, explicit)
})
