## The SCANFI species-layer functions have moved to the `*to` family. `rasterToMatch` and
## `studyArea` still work -- they arrive through `...` and are translated by .legacyToTo() --
## but they are no longer formals, so the ambiguity of "which one did you pass?" is gone from
## the signature itself. Static checks: exercising these functions needs Drive access.

test_that("both SCANFI functions expose the *to family and no longer the legacy pair", {
  for (f in list(loadSCANFISpeciesLayers, prepSpeciesLayers_SCANFI)) {
    nms <- names(formals(f))
    expect_true(all(c("to", "cropTo", "projectTo", "maskTo") %in% nms))
    expect_false("rasterToMatch" %in% nms)
    expect_false("studyArea" %in% nms)
    expect_true("..." %in% nms) ## the legacy pair has to have somewhere to land
  }
})

test_that("the legacy pair still resolves, by the documented table", {
  ## a caller on the old API passing both: the raster gives the geometry, the polygon masks
  res <- .legacyToTo(rasterToMatch = "rtm", studyArea = "sa")
  expect_identical(res$cropTo, "rtm")
  expect_identical(res$projectTo, "rtm")
  expect_identical(res$maskTo, "sa")

  ## and the old shim's mapping, which had it the other way round, is not what we do
  expect_false(identical(res$projectTo, "sa"))
  expect_false(identical(res$maskTo, "rtm"))
})

test_that("loadSCANFISpeciesLayers passes the resolved *to set through to prepInputs", {
  src <- paste(deparse(body(loadSCANFISpeciesLayers)), collapse = " ")
  ## the call used to be `to = rasterToMatch`, which could only ever express one of the cases
  expect_false(grepl("to = rasterToMatch", src, fixed = TRUE))
  expect_true(grepl(".legacyToTo", src, fixed = TRUE))
  for (a in c("cropTo = cropTo", "projectTo = projectTo", "maskTo = maskTo")) {
    expect_true(grepl(a, src, fixed = TRUE))
  }
})

test_that("prepSpeciesLayers_SCANFI forwards the resolved set, not the old inverted pair", {
  src <- paste(deparse(body(prepSpeciesLayers_SCANFI)), collapse = " ")
  expect_true(grepl(".legacyToTo", src, fixed = TRUE))
  ## the previous shim mapped projectTo -> rasterToMatch and to -> studyArea
  expect_false(grepl("projectTo = rasterToMatch", src, fixed = TRUE))
  expect_false(grepl("to = studyArea", src, fixed = TRUE))
})
