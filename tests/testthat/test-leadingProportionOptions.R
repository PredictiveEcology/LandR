## The "leading" threshold -- the share of a stand held by one type above which it stops being
## called mixed -- is read as a nested pair of options:
##
##   getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", <default>))
##
## The outer option is one knob for all of them; unset, each family keeps the default it has
## always had (0.8 for vegetation typing, 0.75 for the land-cover legend, the NTEMS/EOSD value).
## The two inner defaults differ by history, not by concept.

vegFns <- function() list(vegTypeGenerator, vegTypeMapGenerator.data.table, plotVTM)

test_that(".onLoad sets the inner options but NOT the outer one", {
  expect_identical(getOption("LandR.vegLeadingProportion"), 0.8)
  expect_identical(getOption("LandR.lccLeadingProportion"), 0.75)
  ## Setting the outer option at load would consume the outer slot and the fallthrough could
  ## never happen, so it must stay unset.
  expect_null(getOption("NTEMS.mixedwoodProp"))
})

test_that("with the outer option unset, each family falls through to its own default", {
  withr::local_options(list(NTEMS.mixedwoodProp = NULL))
  for (f in vegFns()) {
    expect_identical(eval(formals(f)$vegLeadingProportion), 0.8)
  }
  expect_identical(eval(formals(lccMapGenerator)$vegLeadingProportion), 0.75)
})

test_that("NTEMS.mixedwoodProp overrides every family at once", {
  withr::local_options(list(NTEMS.mixedwoodProp = 0.65))
  for (f in vegFns()) {
    expect_identical(eval(formals(f)$vegLeadingProportion), 0.65)
  }
  expect_identical(eval(formals(lccMapGenerator)$vegLeadingProportion), 0.65)
})

test_that("the outer option wins over an inner one that is also set", {
  withr::local_options(list(NTEMS.mixedwoodProp = 0.65,
                            LandR.vegLeadingProportion = 0.55,
                            LandR.lccLeadingProportion = 0.95))
  expect_identical(eval(formals(vegTypeGenerator)$vegLeadingProportion), 0.65)
  expect_identical(eval(formals(lccMapGenerator)$vegLeadingProportion), 0.65)
})

test_that("an inner option moves only its own family", {
  withr::local_options(list(NTEMS.mixedwoodProp = NULL, LandR.vegLeadingProportion = 0.6))
  expect_identical(eval(formals(vegTypeGenerator)$vegLeadingProportion), 0.6)
  expect_identical(eval(formals(lccMapGenerator)$vegLeadingProportion), 0.75)
})

## Not just the formals: the option has to reach the classification itself.
test_that("changing the option changes what vegTypeGenerator calls leading", {
  withr::local_package("data.table")
  ## 70% Pice_gla / 30% Popu_tre: mixed at 0.8, conifer-leading at 0.6
  x <- data.table(pixelGroup = 1L, speciesCode = c("Pice_gla", "Popu_tre"), B = c(70, 30))
  sppEquiv <- data.table(LandR = c("Pice_gla", "Popu_tre"), Type = c("Conifer", "Deciduous"),
                         EN_generic_short = c("Pice_gla", "Popu_tre"))
  vtg <- function(...) {
    vegTypeGenerator(copy(x), sppEquiv = sppEquiv, sppEquivCol = "LandR",
                     mixedType = 2, doAssertion = FALSE, ...)
  }

  atDefault <- withr::with_options(list(NTEMS.mixedwoodProp = NULL), vtg())
  viaOuter <- withr::with_options(list(NTEMS.mixedwoodProp = 0.6), vtg())
  expect_false(identical(as.character(atDefault$leading), as.character(viaOuter$leading)))
})
