## The legacy pair is ambiguous -- what it asks for depends on which of the two was supplied --
## so the translation follows the table in ?reproducible::postProcess rather than any of the
## hand-rolled variants that had accumulated around the package. The helper never inspects the
## objects, so plain placeholders are enough here.

rtm <- "the-raster"
sa <- "the-polygon"

test_that("rasterToMatch alone is simply `to`", {
  res <- .legacyToTo(rasterToMatch = rtm)
  expect_identical(res$to, rtm)
  expect_null(res$cropTo)
  expect_null(res$projectTo)
  expect_null(res$maskTo)
})

test_that("maskWithRTM = FALSE omits the mask that `to` would otherwise apply", {
  res <- .legacyToTo(rasterToMatch = rtm, maskWithRTM = FALSE)
  expect_identical(res$to, rtm)
  expect_true(is.na(res$maskTo))
})

test_that("studyArea alone crops and masks but does not reproject", {
  res <- .legacyToTo(studyArea = sa)
  expect_null(res$to)
  expect_identical(res$cropTo, sa)
  expect_identical(res$maskTo, sa)
  ## the documented "projection: no" -- NA is reproducible's "omit this step"
  expect_true(is.na(res$projectTo))
})

test_that("useSAcrs lets the study area supply the CRS", {
  res <- .legacyToTo(studyArea = sa, useSAcrs = TRUE)
  expect_identical(res$projectTo, sa)
  expect_identical(res$cropTo, sa)
  expect_identical(res$maskTo, sa)
})

test_that("both: the raster gives the geometry, the polygon gives the mask", {
  res <- .legacyToTo(rasterToMatch = rtm, studyArea = sa)
  expect_null(res$to)
  expect_identical(res$cropTo, rtm)
  expect_identical(res$projectTo, rtm)
  expect_identical(res$maskTo, sa)
  ## the case most easily got wrong: masking to the raster instead of the study area
  expect_false(identical(res$maskTo, rtm))
})

test_that("neither leaves everything untouched", {
  res <- .legacyToTo()
  expect_true(all(vapply(res, is.null, logical(1))))
  expect_named(res, c("to", "cropTo", "projectTo", "maskTo"))
})

## A caller that has already migrated must not be second-guessed by a legacy argument that
## arrives through `...` from somewhere further up.
test_that("an explicit *to argument always beats a legacy one", {
  res <- .legacyToTo(to = "explicit", rasterToMatch = rtm, studyArea = sa)
  expect_identical(res$to, "explicit")

  res2 <- .legacyToTo(maskTo = "explicit", rasterToMatch = rtm, studyArea = sa)
  expect_identical(res2$maskTo, "explicit")
  expect_identical(res2$cropTo, rtm) ## the rest still translated
})

test_that("a fully explicit caller is returned unchanged", {
  res <- .legacyToTo(to = "a", cropTo = "b", projectTo = "c", maskTo = "d",
                     rasterToMatch = rtm, studyArea = sa)
  expect_identical(unlist(res), c(to = "a", cropTo = "b", projectTo = "c", maskTo = "d"))
})
