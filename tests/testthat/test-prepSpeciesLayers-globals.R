## `year` is not a formal of prepSpeciesLayers_SCANFI(); its Google Drive fallback
## (taken when RCurl::url.exists() fails, e.g. during a network blip) referenced it
## anyway, so it resolved to data.table::year and the call died with
## "cannot coerce type 'closure' to vector of type 'character'".

test_that("prepSpeciesLayers_SCANFI does not reference an unbound `year`", {
  skip_if_not_installed("codetools")
  vars <- codetools::findGlobals(prepSpeciesLayers_SCANFI, merge = FALSE)$variables
  expect_false("year" %in% vars)
  expect_true("dataYear" %in% names(formals(prepSpeciesLayers_SCANFI)))
})
