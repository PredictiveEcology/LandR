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

## `$` on a LIST partial-matches, so `dots$studyArea` silently returns a `studyAreaName` that was
## passed through `...`. Biomass_speciesData passes `studyAreaName` and no `studyArea`, so every
## SCANFI species-layer call took the studyArea-only branch and set cropTo to the character
## "6.2.1" -- "cropTo must be a Raster*, Spat*, sf or Spatial object", after ~3 h of simInit, for
## all ten ELFs of a fitting wave (2026-09-17). Reading legacy names exactly is the fix.
test_that("a name that merely PREFIX-matches a legacy argument is not treated as one", {
  expect_null(.legacyDot(list(studyAreaName = "6.2.1"), "studyArea"))
  expect_null(.legacyDot(list(rasterToMatchLarge = "big"), "rasterToMatch"))
  ## the trap itself, so the reason for the helper stays visible
  expect_identical(list(studyAreaName = "6.2.1")$studyArea, "6.2.1")
  ## an exact name still comes through, and an absent one is NULL
  expect_identical(.legacyDot(list(studyArea = sa), "studyArea"), sa)
  expect_null(.legacyDot(list(), "studyArea"))
  ## an explicit NULL is still absent, not an error
  expect_null(.legacyDot(list(studyArea = NULL), "studyArea"))
})

test_that("passing only studyAreaName leaves the *to arguments untouched", {
  ## what the SCANFI species-layer functions compute from their dots
  dots <- list(studyAreaName = "6.2.1", outputPath = "somewhere")
  res <- .legacyToTo(to = sa, cropTo = NULL, projectTo = rtm, maskTo = NULL,
                     rasterToMatch = .legacyDot(dots, "rasterToMatch"),
                     studyArea = .legacyDot(dots, "studyArea"))
  expect_identical(res$to, sa)
  expect_identical(res$projectTo, rtm)
  expect_null(res$cropTo)   ## was "6.2.1"
  expect_null(res$maskTo)   ## was "6.2.1"
})

test_that("both SCANFI species-layer functions read their legacy dots exactly", {
  files <- file.path("..", "..", "R", c("maps.R", "prepSpeciesLayers.R"))
  files <- files[file.exists(files)]
  skip_if_not(length(files) == 2L, "package source not available (installed-package test run)")
  code <- unlist(lapply(files, readLines, warn = FALSE))
  callSites <- grep("\\.legacyToTo\\(", code)
  expect_gte(length(callSites), 2L)
  near <- unlist(lapply(callSites, function(i) code[i:min(i + 6L, length(code))]))
  expect_false(any(grepl("dots\\$(studyArea|rasterToMatch)\\b", near)))
})
