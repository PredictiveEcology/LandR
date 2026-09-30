## Thread count must not change LANDISDisp results.

test_that("LANDISDisp output is identical for threads = 1, 2, 4 and the R reference", {
  skip_if_not_installed("terra")
  ## the OpenMP path only engages above ~20000 active receivers per step
  fix <- makeLANDISDispFixture(size = "large", fixtureSeed = 11L, successionTimestep = 10L)
  ompAvail <- LandR:::landisDispHasOpenMP()

  run <- function(nThr, useCpp = TRUE) {
    withr::with_options(list(LandR.LANDISDisp.threads = nThr),
                        runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = useCpp))
  }
  ref <- run(1L)
  expect_gt(nrow(ref), 0L)
  expect_identical(run(1L, useCpp = FALSE), ref)
  for (n in c(2L, 4L)) {
    skip_if(!ompAvail, "OpenMP not available")
    expect_identical(run(n), ref, info = paste("threads =", n))
  }
})
