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

## `projectTo = rasterToMatch` was passed twice to loadSCANFISpeciesLayers(). R allows
## duplicate names in `...`, so nothing errored and the second copy simply sat there --
## exactly the kind of thing that survives review. Checked statically because calling the
## function needs Drive access.
test_that("the loadSCANFISpeciesLayers() call passes no argument twice", {
  calls <- as.list(body(prepSpeciesLayers_SCANFI))
  theCall <- NULL
  findCall <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("loadSCANFISpeciesLayers"))) theCall <<- x
      lapply(as.list(x), findCall)
    }
    invisible(NULL)
  }
  lapply(calls, findCall)
  expect_false(is.null(theCall))

  argNames <- names(as.list(theCall))[-1]
  argNames <- argNames[nzchar(argNames)]
  expect_identical(anyDuplicated(argNames), 0L)
})

## The cache tag is what Cache() searches on, so tagging SCANFI layers "KNN" made them
## indistinguishable from kNN ones.
test_that("prepSpeciesLayers_SCANFI tags its cache entry SCANFI, not KNN", {
  src <- paste(deparse(body(prepSpeciesLayers_SCANFI)), collapse = " ")
  expect_true(grepl('"speciesLayers", "SCANFI"', src, fixed = TRUE))
  expect_false(grepl('"speciesLayers", "KNN"', src, fixed = TRUE))
})
