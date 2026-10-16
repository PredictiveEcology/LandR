## updateCohortData()'s wholesale re-grouping, calcSiteShade() and rmMissingCohorts() were made faster
## without changing a result: the same cohortData (rows in the same order, same columns), the same
## pixelGroup numbers on the map, the same RNG state. The expected values are what the code before the
## change returned on the inputs in helper-regroup.R.

expected <- readRDS(test_path("fixtures", "regroup-expected.rds"))
got <- regroupRun()

test_that("updateCohortData() re-groups as before: pixelGroup numbers, row order, columns, RNG state", {
  for (nm in names(expected$update)) {
    expect_identical(got$update[[nm]], expected$update[[nm]], label = nm)
  }
})

test_that("the re-grouping merges identical communities and numbers the pixelGroups without gaps", {
  i <- regroupInputs(1)
  before <- terra::values(i$pixelGroupMap, mat = FALSE)
  after <- got$update$plain$pixelGroupMap
  expect_identical(is.na(after), is.na(before))
  ## pixelGroups 5 and 6 are the same community, so their pixels without new cohorts share one number afterwards
  untouched <- setdiff(which(before %in% 5:6), regroupNewCohorts(i)$pixelIndex)
  expect_gt(length(untouched), 1L)
  expect_equal(length(unique(after[untouched])), 1L)
  treed <- sort(unique(after[!is.na(after) & after > 0]))
  expect_identical(treed, seq_along(treed))
  expect_identical(sort(unique(got$update$plain$cohortData$pixelGroup)), treed)
})

test_that("calcSiteShade() gives what it gave before, and leaves the same prevMortality and sumB in cohortData", {
  expect_identical(got$siteShade, expected$siteShade)
})

test_that("rmMissingCohorts() gives what it gave before", {
  expect_identical(got$rmMissing, expected$rmMissing)
})
