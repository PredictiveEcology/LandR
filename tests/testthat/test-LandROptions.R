test_that("LandROptions() defaults are what .onLoad() sets", {
  opts <- LandROptions()
  isNull <- vapply(opts, is.null, logical(1))

  for (nm in names(opts)[!isNull]) {
    expect_identical(getOption(nm), opts[[nm]], info = nm)
  }
})

test_that("NTEMS.mixedwoodProp is a member of LandROptions() but stays unset", {
  ## It is the outer option of
  ##   getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", <default>))
  ## so setting it would consume that slot and no inner default could ever be reached. A NULL
  ## default documents the option while leaving it unset, because options() ignores a NULL.
  opts <- LandROptions()

  expect_true("NTEMS.mixedwoodProp" %in% names(opts))
  expect_null(opts[["NTEMS.mixedwoodProp"]])
  expect_false("NTEMS.mixedwoodProp" %in% names(options()))
})
