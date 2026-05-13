## Seed-locked golden tests for LANDISDisp.
##
## These tests are bit-for-bit comparisons against snapshots captured by
## .dev_landisdisp/build-baselines.R. They guard against any silent change to
## the LANDISDisp algorithm — including the upcoming Rcpp port, where the C++
## implementation must reproduce the *exact* same RNG-driven output as the
## previous R reference implementation.
##
## To rebuild snapshots intentionally:
##   Rscript .dev_landisdisp/build-baselines.R
## (only do this when the algorithm is genuinely changing).

test_that("LANDISDisp matches seed-locked snapshots", {
  skip_if_not_installed("withr")
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")

  fixturesDir <- testthat::test_path("fixtures")
  if (!dir.exists(fixturesDir)) skip("fixtures directory not found")

  specs <- list(
    list(size = "tiny",   fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "tiny",   fixtureSeed = 11L, runSeed = 1729L, ts =  1L),
    list(size = "tiny",   fixtureSeed = 23L, runSeed =   42L, ts =  1L),
    list(size = "small",  fixtureSeed = 11L, runSeed =   42L, ts = 10L),
    list(size = "small",  fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
    list(size = "medium", fixtureSeed = 11L, runSeed =   42L, ts = 10L)
  )

  for (s in specs) {
    info <- sprintf("size=%s fixSeed=%d runSeed=%d ts=%d",
                    s$size, s$fixtureSeed, s$runSeed, s$ts)
    fname <- sprintf("LANDISDisp_%s_fix%d_run%d_ts%d.rds",
                     s$size, s$fixtureSeed, s$runSeed, s$ts)
    fpath <- file.path(fixturesDir, fname)
    if (!file.exists(fpath)) {
      skip(sprintf("missing baseline %s — regenerate with build-baselines.R",
                   fname))
    }
    fix <- makeLANDISDispFixture(size = s$size, fixtureSeed = s$fixtureSeed,
                                 successionTimestep = s$ts)
    actual <- runLANDISDispOnFixture(fix, runSeed = s$runSeed)
    expected <- readRDS(fpath)
    expect_equal(actual, expected, info = info)
  }
})

test_that("LANDISDisp is deterministic for a fixed seed (run twice = same)", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")

  fix <- makeLANDISDispFixture(size = "tiny", fixtureSeed = 11L)
  out1 <- runLANDISDispOnFixture(fix, runSeed = 42L)
  out2 <- runLANDISDispOnFixture(fix, runSeed = 42L)
  expect_equal(out1, out2)
})

test_that("LANDISDisp returns different output for different seeds", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
  out1 <- runLANDISDispOnFixture(fix, runSeed = 42L)
  out2 <- runLANDISDispOnFixture(fix, runSeed = 1729L)
  expect_false(isTRUE(all.equal(out1, out2)))
})

test_that("LANDISDisp matches snapshots on rcv-heavy xlarge fixtures (slow)", {
  ## Gated: takes ~30s. Run with: LANDR_SLOW_TESTS=1 R CMD check ...
  if (!identical(Sys.getenv("LANDR_SLOW_TESTS"), "1")) {
    skip("set LANDR_SLOW_TESTS=1 to enable rcv-heavy xlarge tests")
  }
  skip_if_not_installed("withr")
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")

  fixturesDir <- testthat::test_path("fixtures")
  if (!dir.exists(fixturesDir)) skip("fixtures directory not found")

  specs <- list(
    list(size = "xlarge_dense", fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "xlarge_dense", fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
    list(size = "xxlarge",      fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "xxxlarge",     fixtureSeed = 11L, runSeed =   42L, ts =  1L)
  )

  for (s in specs) {
    info <- sprintf("size=%s fixSeed=%d runSeed=%d ts=%d",
                    s$size, s$fixtureSeed, s$runSeed, s$ts)
    fname <- sprintf("LANDISDisp_%s_fix%d_run%d_ts%d.rds",
                     s$size, s$fixtureSeed, s$runSeed, s$ts)
    fpath <- file.path(fixturesDir, fname)
    if (!file.exists(fpath)) {
      skip(sprintf("missing baseline %s — regenerate with build-baselines.R",
                   fname))
    }
    fix <- makeLANDISDispFixture(size = s$size, fixtureSeed = s$fixtureSeed,
                                 successionTimestep = s$ts)
    actual <- runLANDISDispOnFixture(fix, runSeed = s$runSeed)
    expected <- readRDS(fpath)
    expect_equal(actual, expected, info = info)
  }
})

test_that("Rcpp and R implementations produce bit-identical output", {
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")
  skip_if_not_installed("Rcpp")
  ## Skip if the package's compiled library isn't loaded (e.g., in a
  ## bare-source workflow without sourceCpp). Detection: useCpp path will
  ## error if spiralLoopCpp doesn't exist.
  haveCpp <- tryCatch(
    is.function(get("spiralLoopCpp", envir = globalenv(), inherits = TRUE)),
    error = function(e) FALSE
  )
  if (!haveCpp) skip("spiralLoopCpp not available (compiled library not loaded)")

  specs <- list(
    list(size = "tiny",   fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "tiny",   fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
    list(size = "tiny",   fixtureSeed = 23L, runSeed =   42L, ts =  1L),
    list(size = "small",  fixtureSeed = 11L, runSeed =   42L, ts = 10L),
    list(size = "small",  fixtureSeed = 11L, runSeed = 1729L, ts =  1L),
    list(size = "medium", fixtureSeed = 11L, runSeed =   42L, ts = 10L)
  )
  if (identical(Sys.getenv("LANDR_SLOW_TESTS"), "1")) {
    specs <- c(specs, list(
      list(size = "xlarge_dense", fixtureSeed = 11L, runSeed =   42L, ts =  1L),
      list(size = "xxlarge",      fixtureSeed = 11L, runSeed =   42L, ts =  1L),
      list(size = "xxxlarge",     fixtureSeed = 11L, runSeed =   42L, ts =  1L)
    ))
  }

  for (s in specs) {
    info <- sprintf("size=%s fixSeed=%d runSeed=%d ts=%d",
                    s$size, s$fixtureSeed, s$runSeed, s$ts)
    fix <- makeLANDISDispFixture(size = s$size, fixtureSeed = s$fixtureSeed,
                                 successionTimestep = s$ts)
    outR   <- runLANDISDispOnFixture(fix, runSeed = s$runSeed, useCpp = FALSE)
    outCpp <- runLANDISDispOnFixture(fix, runSeed = s$runSeed, useCpp = TRUE)
    expect_equal(outCpp, outR, info = info)
  }
})
