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

test_that("each elm row carries its own KNN name", {
  ## The KNN column was shifted up one row across the Ulmus block: Ulmus pumila carried
  ## `Ulmu_Rub`, Ulmus rubra carried `Ulmu_Spp` and Ulmus spp. carried `Ulmu_Tho`, none of
  ## which are kNN layers. Rock elm, the row those names belong to, had no LandR name at all,
  ## so equivalentName("Ulmu_Tho", column = "LandR") returned the elm genus, `Ulmu_spp`.
  sppEquiv <- data.table::as.data.table(LandR::sppEquivalencies_CA)
  elms <- sppEquiv[startsWith(LandR, "Ulmu"), c("Latin_full", "LandR", "KNN")]

  expect_identical(elms[Latin_full == "Ulmus pumila", KNN], "")
  expect_identical(elms[Latin_full == "Ulmus rubra", KNN], "Ulmu_Rub")
  expect_identical(elms[Latin_full == "Ulmus spp.", KNN], "Ulmu_Spp")
  expect_identical(elms[Latin_full == "Ulmus thomasii", KNN], "Ulmu_Tho")

  expect_identical(equivalentName("Ulmu_Tho", sppEquiv, column = "LandR"), "Ulmu_tho")
  expect_identical(equivalentName("Ulmu_Rub", sppEquiv, column = "LandR"), "Ulmu_rub")
})

test_that("rock elm and pagoda dogwood have LandR names", {
  ## Both are single species with a full set of other names; a blank LandR name made them
  ## unmatchable, as LandR is the column rows are keyed on (e.g. in speciesInStudyArea()).
  sppEquiv <- data.table::as.data.table(LandR::sppEquivalencies_CA)

  expect_identical(sppEquiv[NFI == "ULMU_THO", LandR], "Ulmu_tho")
  expect_identical(sppEquiv[NFI == "CORN_ALT", LandR], "Corn_alt")
})
