## Two thresholds, two questions (see ?leadingProportions):
##
##   LandR.mixedwoodProp       all conifers vs all broadleaves -- the NTEMS/EOSD/NFI rule
##   LandR.leadingSpeciesProp  one species' share; unset, it takes the mixedwood value
##
## Only `LandROptions()` writes a number down; every function reads it through the two accessors.

vegFns <- function() list(vegTypeGenerator, vegTypeMapGenerator.data.table)

test_that(".onLoad sets the mixedwood option, and leadingSpeciesProp inherits from it", {
  expect_identical(getOption("LandR.mixedwoodProp"), 0.75)
  ## Unset on purpose: set at load, it could not track a later change to LandR.mixedwoodProp.
  expect_null(getOption("LandR.leadingSpeciesProp"))
  expect_identical(mixedwoodProp(), 0.75)
  expect_identical(leadingSpeciesProp(), 0.75)
})

test_that("moving LandR.mixedwoodProp moves the leading threshold with it", {
  withr::local_options(list(LandR.mixedwoodProp = 0.6, LandR.leadingSpeciesProp = NULL))
  expect_identical(mixedwoodProp(), 0.6)
  expect_identical(leadingSpeciesProp(), 0.6)
})

test_that("LandR.leadingSpeciesProp can be set on its own, e.g. 'just a majority'", {
  withr::local_options(list(LandR.leadingSpeciesProp = 0.51))
  expect_identical(leadingSpeciesProp(), 0.51)
  expect_identical(mixedwoodProp(), 0.75)   ## unchanged
})

test_that("the accessors survive both options being cleared", {
  withr::local_options(list(LandR.mixedwoodProp = NULL, LandR.leadingSpeciesProp = NULL))
  expect_identical(mixedwoodProp(), 0.75)
  expect_identical(leadingSpeciesProp(), 0.75)
})

test_that("each caller takes the threshold its question needs", {
  ## mixedType = 2 asks the mixedwood question; every other mixedType asks about one species.
  withr::local_options(list(LandR.mixedwoodProp = 0.75, LandR.leadingSpeciesProp = NULL))
  expect_identical(.leadingProp(NULL, mixedType = 2), 0.75)
  expect_identical(.leadingProp(NULL, mixedType = 1), 0.75)

  withr::local_options(list(LandR.leadingSpeciesProp = 0.51))
  expect_identical(.leadingProp(NULL, mixedType = 1), 0.51)
  expect_identical(.leadingProp(NULL, mixedType = 0), 0.51)
  expect_warning(expect_identical(.leadingProp(NULL, mixedType = 2), 0.75), "mixedType = 2")

  ## An explicit argument is never ambiguous, so it neither warns nor consults an option.
  expect_silent(expect_identical(.leadingProp(0.9, mixedType = 2), 0.9))
})

test_that("lccMapGenerator and plotVTM default to the mixedwood threshold", {
  withr::local_options(list(LandR.mixedwoodProp = 0.65))
  expect_identical(eval(formals(lccMapGenerator)$vegLeadingProportion), 0.65)
  expect_identical(eval(formals(plotVTM)$vegLeadingProportion), 0.65)
})

test_that("the vegetation-typing functions take no threshold in their formals", {
  ## They resolve it in the body, because which option applies depends on mixedType.
  for (f in vegFns()) {
    expect_null(eval(formals(f)$vegLeadingProportion))
  }
})

## ---- the classification itself, not just the plumbing ----------------------------------------

sppEquivTest <- function() {
  data.table::data.table(
    LandR = c("Pice_gla", "Lari_lar", "Popu_tre", "Betu_pap", "Acer_rub"),
    Type = c("Conifer", "Conifer", "Deciduous", "Deciduous", "Deciduous"),
    EN_generic_short = c("Pice_gla", "Lari_lar", "Popu_tre", "Betu_pap", "Acer_rub")
  )
}

leadingOf <- function(B, spp, ...) {
  x <- data.table::data.table(pixelGroup = 1L, speciesCode = spp, B = B)
  out <- vegTypeGenerator(x, sppEquiv = sppEquivTest(), sppEquivCol = "LandR",
                          mixedType = 2, doAssertion = FALSE, ...)
  as.character(unique(out$leading))
}

test_that("mixedType = 2 sums the broadleaf group rather than testing one species at a time", {
  withr::local_package("data.table")
  ## 55% conifer, 45% broadleaf split over THREE species: no single one reaches 0.25, but the
  ## group does, so this is mixedwood. Before the fix it was called pure Pice_gla.
  expect_identical(
    leadingOf(c(55, 15, 15, 15), c("Pice_gla", "Popu_tre", "Betu_pap", "Acer_rub")),
    "Mixed"
  )
  ## and the group rule still calls a genuinely conifer-dominated stand conifer
  expect_identical(
    leadingOf(c(80, 8, 6, 6), c("Pice_gla", "Popu_tre", "Betu_pap", "Acer_rub")),
    "Pice_gla"
  )
})

test_that("Larix is a conifer: it never pushes a stand toward mixedwood", {
  withr::local_package("data.table")
  ## 50/50 white spruce and tamarack is 0% broadleaf -- conifer, not mixed.
  expect_false(identical(leadingOf(c(50, 50), c("Pice_gla", "Lari_lar")), "Mixed"))
  ## and sppEquivalencies_CA agrees, which is what the merge on `Type` relies on
  sppEq <- data.table::as.data.table(sppEquivalencies_CA)
  expect_true(all(sppEq[grepl("^Lari_", LandR), Type] == "Conifer"))
})

test_that("a species missing from sppEquiv is not counted as broadleaf", {
  withr::local_package("data.table")
  ## Pinu_ban is absent from the test sppEquiv, so its Type is NA after the merge. It must not
  ## poison the sum: 60/40 conifer/absent is not mixedwood on broadleaf share alone.
  expect_false(identical(leadingOf(c(60, 40), c("Pice_gla", "Pinu_ban")), "Mixed"))
})

test_that("changing LandR.mixedwoodProp changes what is called mixedwood", {
  withr::local_package("data.table")
  ## 70% conifer / 30% broadleaf: mixedwood at 0.75, conifer-leading at 0.65
  atDefault <- withr::with_options(
    list(LandR.mixedwoodProp = 0.75), leadingOf(c(70, 30), c("Pice_gla", "Popu_tre"))
  )
  lower <- withr::with_options(
    list(LandR.mixedwoodProp = 0.65), leadingOf(c(70, 30), c("Pice_gla", "Popu_tre"))
  )
  expect_identical(atDefault, "Mixed")
  expect_identical(lower, "Pice_gla")
})

test_that("mixedType = 1 uses the single-species threshold", {
  withr::local_package("data.table")
  x <- data.table(pixelGroup = 1L, speciesCode = c("Pice_gla", "Popu_tre"), B = c(60, 40))
  vtg1 <- function() {
    as.character(unique(vegTypeGenerator(copy(x), sppEquiv = sppEquivTest(),
                                         sppEquivCol = "LandR", mixedType = 1,
                                         doAssertion = FALSE)$leading))
  }
  ## 0.6 < 0.75, so the top species is not leading: Mixed
  expect_identical(withr::with_options(list(LandR.leadingSpeciesProp = NULL), vtg1()), "Mixed")
  ## "just a majority" makes it leading
  expect_identical(withr::with_options(list(LandR.leadingSpeciesProp = 0.51), vtg1()), "Pice_gla")
})

## ---- subsample size -------------------------------------------------------------------------

test_that("subsetDataSize() is 500 and comes from LandROptions()", {
  expect_identical(subsetDataSize(), 500L)
  expect_identical(LandROptions()[["LandR.subsetDataSize"]], 500L)
})

test_that("options(LandR.subsetDataSize) moves what subsetDT() keeps", {
  withr::local_package("data.table")
  dt <- data.table(grp = rep(letters[1:2], each = 40), x = 1:80)
  withr::local_options(list(LandR.subsetDataSize = 10))
  expect_identical(subsetDataSize(), 10)
  out <- suppressMessages(subsetDT(dt, by = "grp"))   ## doSubset = TRUE -> the option
  expect_identical(nrow(out), 20L)                    ## 10 per group
  ## an explicit number still wins
  expect_identical(nrow(suppressMessages(subsetDT(dt, by = "grp", doSubset = 5))), 10L)
})
