test_that("LandROptions() defaults are what .onLoad() sets", {
  opts <- LandROptions()
  isNull <- vapply(opts, is.null, logical(1))

  for (nm in names(opts)[!isNull]) {
    expect_identical(getOption(nm), opts[[nm]], info = nm)
  }
})

test_that("LandR.leadingSpeciesProp is a member of LandROptions() but stays unset", {
  ## Unset, leadingSpeciesProp() falls through to LandR.mixedwoodProp, so the two track each
  ## other and moving the mixedwood value moves both. Setting it at load would fix it there. A
  ## NULL default documents the option while leaving it unset, because options() ignores a NULL.
  opts <- LandROptions()

  expect_true("LandR.leadingSpeciesProp" %in% names(opts))
  expect_null(opts[["LandR.leadingSpeciesProp"]])
  expect_false("LandR.leadingSpeciesProp" %in% names(options()))
})

test_that("LandR.forestLandUseCache is a member of LandROptions() but stays unset", {
  ## Unset, the forest-land inputs follow reproducible.useCache, like any other Cache() call.
  opts <- LandROptions()

  expect_true("LandR.forestLandUseCache" %in% names(opts))
  expect_null(opts[["LandR.forestLandUseCache"]])
  expect_false("LandR.forestLandUseCache" %in% names(options()))

  withr::local_options(LandR.forestLandUseCache = NULL, reproducible.useCache = FALSE)
  expect_identical(LandR:::.forestLandUseCache(), FALSE)
  withr::local_options(LandR.forestLandUseCache = "always")
  expect_identical(LandR:::.forestLandUseCache(), "always")
})
