## leadingByPixelGroup() (behind vegTypeGenerator() and vegTypeMapGenerator()), .initiateNewCohorts(),
## calcSeverityB() and rmMissingCohorts() were made faster without changing a result. The expected
## values are what the code before the change returned on the inputs in helper-fastCohortOps.R.

expected <- readRDS(test_path("fixtures", "fastCohortOps-expected.rds"))
got <- fastCohortOpsRun()

test_that("vegTypeGenerator() and vegTypeMapGenerator() give what they gave before (with the RNG state)", {
  skipIfNoLongDouble()
  for (nm in names(expected$leading)) {
    expect_identical(got$leading[[nm]], expected$leading[[nm]], label = nm)
  }
})

test_that("updateCohortData() and .initiateNewCohorts() give what they gave before, for initialB given, NA and NULL", {
  for (nm in names(expected$regen)) {
    expect_identical(got$regen[[nm]], expected$regen[[nm]], label = nm)
  }
})

test_that("calcSeverityB() and rmMissingCohorts() give what they gave before", {
  expect_identical(got$severity, expected$severity)
  expect_identical(got$rmMissing, expected$rmMissing)
})

test_that("leadingByPixelGroup() is the same as the original algorithm, per pixelGroup", {
  cd <- fastOpsCohorts()
  for (cdi in list(cd, fastOpsCohorts(bigPixelGroup = TRUE))) {
    expect_true(.canSumLeadingByIntegerKey(cdi, 2, fastOpsSppEquiv(), "LandR", "pixelGroup", "B"))
    fast <- leadingByPixelGroup(cdi, 0.75, 2, fastOpsSppEquiv(), "LandR")
    slow <- .leadingByPixelGroupSlow(cdi, 0.75, 2, fastOpsSppEquiv(), "LandR", "pixelGroup", "B")
    expect_identical(fast$pixelGroup, slow$pixelGroup)
    expect_identical(fast$leading, slow$leading)
    expect_identical(fast$mixed, slow$mixed)
    expect_identical(fast$speciesGroupB, slow$speciesGroupB)
    expect_identical(names(fast), names(slow))
  }
  ## the factor is in character order, "Mixed" first
  expect_identical(levels(leadingByPixelGroup(cd, 0.75, 2, fastOpsSppEquiv(), "LandR")$leading),
                   c("Acer_rub", "Betu_pap", "Lari_lar", "Mixed", "Pice_gla", "Popu_tre"))
})

test_that("a mixedwood decided at the threshold by two or more deciduous species is not decided by summing order", {
  skipIfNoLongDouble()
  ## pixelGroups 8, 9 and 13: the deciduous proportions add up to the threshold in base::sum() (what
  ## .isMixedwood() uses), but a sum in double precision lands on the other side of it
  b <- function(pg, spp) fastOpsCohorts()[pixelGroup == pg & speciesCode %in% spp]$B
  sumBase <- function(pg) {
    cd <- fastOpsCohorts()[pixelGroup == pg]
    cd <- cd[order(as.character(speciesCode))]
    p <- cd$B / sum(cd$B)
    dec <- p[as.character(cd$speciesCode) %in% c("Acer_rub", "Betu_pap", "Popu_tre")]
    list(base = sum(dec), double = Reduce(`+`, dec))
  }
  s8 <- sumBase(8L)
  expect_true(s8$base == 0.75 && s8$double < 0.75)
  s13 <- sumBase(13L)
  expect_true(s13$base == 0.25 && s13$double > 0.25)

  lead <- leadingByPixelGroup(fastOpsCohorts(), 0.75, 2, fastOpsSppEquiv(), "LandR")
  setkey(lead, pixelGroup)
  expect_identical(lead[.(c(8L, 9L, 13L)), as.character(leading)], c("Acer_rub", "Pice_gla", "Pice_gla"))
  expect_false(any(lead[.(c(8L, 9L, 13L)), mixed]))
  expect_true(lead[.(10L), mixed]) ## three deciduous species at 15% each
})

test_that("leadingByPixelGroup() takes the original algorithm when the fast one does not apply", {
  sppEq <- fastOpsSppEquiv()
  cd <- fastOpsCohorts()
  expect_false(.canSumLeadingByIntegerKey(cd, 1, sppEq, "LandR", "pixelGroup", "B"))
  expect_false(.canSumLeadingByIntegerKey(cd, 2, sppEq, "LandR", "pixelGroup", "cover"))
  expect_false(.canSumLeadingByIntegerKey(copy(cd)[, speciesCode := as.character(speciesCode)],
                                          2, sppEq, "LandR", "pixelGroup", "B"))
  expect_false(.canSumLeadingByIntegerKey(copy(cd)[1, B := NA_integer_], 2, sppEq, "LandR", "pixelGroup", "B"))
  expect_false(.canSumLeadingByIntegerKey(copy(cd)[1, pixelGroup := 0L], 2, sppEq, "LandR", "pixelGroup", "B"))
  twoTypes <- rbind(sppEq, data.table(LandR = "Pice_gla", Type = "Deciduous"))
  expect_false(.canSumLeadingByIntegerKey(cd, 2, twoTypes, "LandR", "pixelGroup", "B"))
  expect_error(leadingByPixelGroup(cd, 0.75, 2, sppEq, "noSuchColumn"), "not in sppEquiv")
})

test_that("leadingByPixelGroup() works on a pixelGroup column with another name", {
  cd <- fastOpsCohorts()
  setnames(cd, "pixelGroup", "pg")
  a <- leadingByPixelGroup(cd, 0.75, 2, fastOpsSppEquiv(), "LandR", pixelGroupColName = "pg")
  b <- leadingByPixelGroup(fastOpsCohorts(), 0.75, 2, fastOpsSppEquiv(), "LandR")
  expect_identical(a$pg, b$pixelGroup)
  expect_identical(a$leading, b$leading)
})

test_that("the speed-ups leave the RNG state alone", {
  i <- fastOpsRegenInputs()
  burned <- data.table(pixelIndex = 5L, pixelGroup = 2L)
  post <- i$cohortData[burned, on = "pixelGroup", allow.cartesian = TRUE][, B := 0L]
  set.seed(123)
  seed <- get(".Random.seed", envir = globalenv())
  suppressMessages({
    leadingByPixelGroup(fastOpsCohorts(), 0.75, 2, fastOpsSppEquiv(), "LandR")
    vegTypeGenerator(fastOpsCohorts(), 0.75, sppEquiv = fastOpsSppEquiv(), sppEquivCol = "LandR",
                     doAssertion = FALSE)
    vegTypeMapGenerator(fastOpsCohorts(), fastOpsPixelGroupMap(), 0.75, sppEquiv = fastOpsSppEquiv(),
                        sppEquivCol = "LandR", doAssertion = FALSE)
    calcSeverityB(i$cohortData, post)
    rmMissingCohorts(i$cohortData, i$pixelGroupMap, doAssertion = FALSE)
  })
  expect_identical(get(".Random.seed", envir = globalenv()), seed)
})

test_that(".initiateNewCohorts() computes the B of the older cohorts only when initialB is NA or NULL", {
  i <- fastOpsRegenInputs()
  new <- function() {
    data.table(pixelGroup = c(1L, 3L), speciesCode = i$spp[c(2, 2)], ecoregionGroup = i$ecos[c(1, 1)], age = 1L)
  }
  init <- function(initialB) {
    cd <- copy(i$cohortData)
    out <- .initiateNewCohorts(new(), cd, i$pixelGroupMap, currentTime = 7,
                               speciesEcoregion = i$speciesEcoregion, successionTimestep = 10,
                               initialB = initialB)
    expect_false("oldSumB" %in% names(cd))
    out[age == 1L & speciesCode == "Pinu_ban"]
  }
  ## given: B is initialB
  expect_identical(init(10)$B, c(10L, 10L))
  ## NA: B falls with the B already there (pixelGroup 3: 800 of eco1's maxB of 9000) ...
  expect_identical(init(NA)$B[2], 173L)
  ## ... except in pixelGroup 1, whose first row is young: its older B (4500) counts as 0, so B is
  ## maxANPP (200). Known separate issue (unique(by = "pixelGroup") keeps the first row); unchanged here.
  expect_identical(init(NA)$B[1], 200L)
  expect_identical(init(NULL), init(NA))
})
