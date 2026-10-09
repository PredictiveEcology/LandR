test_that("modifySpeciesAndSpeciesEcoregionTable adds no row for a species without traits", {
  ## Biomass_core#124: the join used to be a right join, so a species in speciesTable with no
  ## speciesEcoregion row gained a row that was NA in every column, year included.
  se <- data.table::data.table(
    ecoregionGroup = factor(c("1_01", "1_02")), speciesCode = factor("Pice_mar"),
    establishprob = 0.5, maxB = 5000L, maxANPP = 160L, year = 2020
  )
  spp <- data.table::data.table(
    species = c("Pice_mar", "Abie_las"), hardsoft = "soft", longevity = 200L,
    growthcurve = 0.5, mortalityshape = 15L, inflationFactor = c(1.1, 1.2), mANPPproportion = 3.3
  )
  out <- suppressMessages(modifySpeciesAndSpeciesEcoregionTable(se, spp))$newSpeciesEcoregion
  expect_equal(nrow(out), 2L)
  expect_false(anyNA(out))
  expect_equal(unique(as.character(out$speciesCode)), "Pice_mar")
  expect_equal(out$maxB, rep(5500L, 2))
})
