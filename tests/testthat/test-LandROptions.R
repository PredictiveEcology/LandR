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
