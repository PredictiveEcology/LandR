## Bit-identical parity tests: R reference vs Rcpp port.
##
## These tests run BOTH implementations in the same R session against the same
## seeded inputs and assert the outputs are identical. They are independent of
## the goldens in test-LANDISDisp-seedLocked.R: a regression that affects both
## paths equally would still pass these, but a regression that affects only the
## Cpp path would be caught here. A regression that affects only the R path
## would be caught by the goldens.
##
## Each test_that() block targets a different axis (RNG seed coverage, raster
## edge cases, parameter sweep, multi-call invariance, factor speciesCode,
## etc.) so a failure narrows the suspect surface.

skip_if_no_cpp <- function() {
  haveCpp <- tryCatch(
    is.function(get("spiralLoopCpp", envir = globalenv(), inherits = TRUE)),
    error = function(e) FALSE
  )
  if (!haveCpp) skip("spiralLoopCpp not available (compiled library not loaded)")
}

parityCheck <- function(fix, runSeed, info = "", ...) {
  outR   <- runLANDISDispOnFixture(fix, runSeed = runSeed, useCpp = FALSE, ...)
  outCpp <- runLANDISDispOnFixture(fix, runSeed = runSeed, useCpp = TRUE, ...)
  expect_equal(outCpp, outR, info = info)
}

# ---------------------------------------------------------------------------
# 1. Wide RNG-seed coverage. The current seed-locked tests use ~10 seeds; this
# fuzzes 50 random seeds against tiny + small. Catches divergence that only
# fires in narrow RNG-stream slices.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree across 50 random seeds (tiny + small fixtures)", {
  skip_if_no_cpp()
  skip_if_not_installed("SpaDES.tools")

  ## reproducible seed list
  set.seed(20260501L)
  seeds <- sample.int(.Machine$integer.max, 50L)

  for (sz in c("tiny", "small")) {
    fix <- makeLANDISDispFixture(size = sz, fixtureSeed = 11L)
    for (sd in seeds) {
      parityCheck(fix, runSeed = sd,
                  info = sprintf("size=%s seed=%d", sz, sd))
    }
  }
})

# ---------------------------------------------------------------------------
# 2. pixelGroupMap with NA (masked) cells. Cpp path uses
# `if (pg == NA_INTEGER) continue`, which only fires when the receiver's
# spiral target lands on a masked cell. None of the regular fixtures have NA
# pixels — this introduces them.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree when pixelGroupMap has NA cells", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
  pgm <- terra::rast(fix$pixelGroupMap) # deep copy
  terra::values(pgm) <- terra::values(fix$pixelGroupMap)
  ## punch holes: a stripe and a few scattered cells
  v <- terra::values(pgm, mat = FALSE)
  v[1:50] <- NA_integer_
  v[seq(500, 2000, by = 17)] <- NA_integer_
  terra::values(pgm) <- as.integer(v)
  fix$pixelGroupMap <- pgm

  for (sd in c(1L, 42L, 1729L)) {
    parityCheck(fix, runSeed = sd, info = sprintf("seed=%d", sd))
  }
})

# ---------------------------------------------------------------------------
# 3. Non-contiguous pixelGroup IDs. Bitmask is a flat vector of length
# maxPg + 1; large pg IDs allocate a larger vector but the result must
# match the matrix-based R reference exactly.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree with sparse (non-contiguous) pixelGroup IDs", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)

  ## remap pixelGroup values: 1->10, 2->20, 3->300, ... (gaps up to 30k)
  oldPgs <- sort(unique(c(fix$dtSrc$pixelGroup, fix$dtRcv$pixelGroup,
                          as.vector(fix$pixelGroupMap[]))))
  oldPgs <- oldPgs[!is.na(oldPgs)]
  set.seed(7L)
  newIds <- sort(sample.int(30000L, length(oldPgs)))
  remap <- setNames(newIds, as.character(oldPgs))

  pgmNew <- terra::rast(fix$pixelGroupMap)
  vv <- as.vector(fix$pixelGroupMap[])
  vNew <- remap[as.character(vv)]
  vNew[is.na(vNew)] <- NA_integer_
  terra::values(pgmNew) <- as.integer(vNew)
  fix$pixelGroupMap <- pgmNew
  fix$dtSrc[, pixelGroup := remap[as.character(pixelGroup)]]
  fix$dtRcv[, pixelGroup := remap[as.character(pixelGroup)]]
  fix$dtRcvFull[, pixelGroup := remap[as.character(pixelGroup)]]

  parityCheck(fix, runSeed = 42L, info = "sparse-pg")
})

# ---------------------------------------------------------------------------
# 4. Single-species fixture (numSp = 1). Exercises the smallest bitmask and
# the tightest active-species pruning.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree with a single species", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
  fix$dtSrc <- fix$dtSrc[speciesCode == 6L]                     # Popu_tre only
  fix$dtRcv <- fix$dtRcv[speciesCode == 6L]
  fix$dtRcvFull <- fix$speciesTable[fix$dtRcv, on = "speciesCode"]

  parityCheck(fix, runSeed = 42L, info = "single species")
})

# ---------------------------------------------------------------------------
# 5. dtRcv contains species absent from dtSrc. The species-filter inside
# LANDISDisp should drop those rows; both paths must agree on the result.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree when dtRcv has species not in dtSrc", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
  ## pull all species 1..7 in dtSrc, but only species 1, 2 for dtRcv — opposite
  ## of fixture-default. Use a recipe that keeps dtRcv rows for species the
  ## fixture might not have produced.
  fix$dtRcv <- fix$dtRcv[speciesCode %in% c(1L, 2L)]
  fix$dtRcvFull <- fix$speciesTable[fix$dtRcv, on = "speciesCode"]

  parityCheck(fix, runSeed = 42L, info = "rcv-species-subset")
})

# ---------------------------------------------------------------------------
# 6. Multi-call invariance: five LANDISDisp calls in sequence, each consuming
# the previous's RNG state. Both paths must produce the same five outputs.
# Catches RNG-state leaks (e.g., Get/PutRNGstate misuse in the C++ port).
# ---------------------------------------------------------------------------
test_that("R and Cpp agree across 5 sequential LANDISDisp calls", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)

  callN <- function(useCpp, n = 5L) {
    set.seed(42L)
    lapply(seq_len(n), function(i) {
      out <- LANDISDisp(
        dtSrc = fix$dtSrc, dtRcv = fix$dtRcvFull,
        pixelGroupMap = fix$pixelGroupMap, speciesTable = fix$speciesTable,
        successionTimestep = fix$successionTimestep, verbose = 1, useCpp = useCpp
      )
      data.table::setattr(out, "ReasonForStop", NULL)
      out[, .(pixelIndex, speciesCode, DistOfSuccess, species)][
        order(pixelIndex, speciesCode)]
    })
  }

  rList <- callN(useCpp = FALSE)
  cList <- callN(useCpp = TRUE)
  expect_equal(cList, rList)
  ## sanity: the 5 outputs should differ from each other (RNG advances)
  expect_false(isTRUE(all.equal(rList[[1]], rList[[2]])))
})

# ---------------------------------------------------------------------------
# 7. Factor speciesCode path. By default the fixture uses integer codes;
# LANDISDisp re-factorises them. This explicitly tests the path where the
# caller passes factor speciesCode columns (the more common upstream case).
# ---------------------------------------------------------------------------
test_that("R and Cpp agree with factor speciesCode columns", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
  ## Convert to factor with the speciesTable's species-name labels
  spLevels <- fix$speciesTable$species
  toFactor <- function(dt) {
    set(dt, NULL, "speciesCode",
        factor(spLevels[dt$speciesCode], levels = spLevels))
    dt
  }
  fix$dtSrc <- toFactor(data.table::copy(fix$dtSrc))
  fix$dtRcv <- toFactor(data.table::copy(fix$dtRcv))
  st2 <- data.table::copy(fix$speciesTable)
  set(st2, NULL, "speciesCode", factor(st2$species, levels = spLevels))
  fix$speciesTable <- st2
  fix$dtRcvFull <- st2[fix$dtRcv, on = "speciesCode"]

  parityCheck(fix, runSeed = 42L, info = "factor speciesCode")
})

# ---------------------------------------------------------------------------
# 8. Parameter sweep: vary Ward kernel parameters k, b and successionTimestep.
# Each combination is its own seed-stream; both paths must match exactly.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree across k/b/successionTimestep sweep", {
  skip_if_no_cpp()

  grid <- expand.grid(
    k  = c(0.5, 0.95, 0.99),
    b  = c(0.001, 0.01, 0.1),
    ts = c(1L, 5L, 10L, 25L)
  )
  fix <- makeLANDISDispFixture(size = "tiny", fixtureSeed = 11L)

  for (i in seq_len(nrow(grid))) {
    fixI <- fix
    fixI$successionTimestep <- grid$ts[i]
    info <- sprintf("k=%g b=%g ts=%d", grid$k[i], grid$b[i], grid$ts[i])
    parityCheck(fixI, runSeed = 42L, info = info,
                k = grid$k[i], b = grid$b[i])
  }
})

# ---------------------------------------------------------------------------
# 9. Non-default cellSize (raster resolution). Ward kernel scales with
# distance in raw units; both paths must agree under different cell sizes.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree at non-default cellSize (250m, 1000m)", {
  skip_if_no_cpp()

  for (cs in c(250, 1000)) {
    fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
    pgmNew <- terra::rast(
      xmin = 0, xmax = ncol(fix$pixelGroupMap) * cs,
      ymin = 0, ymax = nrow(fix$pixelGroupMap) * cs,
      resolution = c(cs, cs),
      vals = as.vector(fix$pixelGroupMap[])
    )
    fix$pixelGroupMap <- pgmNew
    parityCheck(fix, runSeed = 42L, info = sprintf("cellSize=%g", cs))
  }
})

# ---------------------------------------------------------------------------
# 10. Multi-fixture stress: sweep 3 fixture seeds × 3 RNG seeds × 2 timesteps
# under both paths. 18 combinations — quick smoke test that catches drift.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree across fixture×RNG×timestep grid", {
  skip_if_no_cpp()

  for (fsd in c(11L, 23L, 99L)) {
    for (rsd in c(1L, 42L, 9999L)) {
      for (ts in c(1L, 10L)) {
        fix <- makeLANDISDispFixture(size = "small", fixtureSeed = fsd,
                                     successionTimestep = ts)
        info <- sprintf("fsd=%d rsd=%d ts=%d", fsd, rsd, ts)
        parityCheck(fix, runSeed = rsd, info = info)
      }
    }
  }
})

# ---------------------------------------------------------------------------
# 11. Verbose-level invariance. Both impls record the per-iteration
# DistOfSuccess only when verbose >= 1, but the Success status (which
# pixelIndex × speciesCode pairs got seed) must be identical across verbose
# levels in each impl. Catches any verbose-only branch that affects the RNG
# sequence (none should, but lock it in).
# ---------------------------------------------------------------------------
test_that("Success status is invariant across verbose levels (R and Cpp)", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)

  successKeys <- function(useCpp, verbose) {
    set.seed(42L)
    out <- LANDISDisp(
      dtSrc = fix$dtSrc, dtRcv = fix$dtRcvFull,
      pixelGroupMap = fix$pixelGroupMap, speciesTable = fix$speciesTable,
      successionTimestep = fix$successionTimestep,
      verbose = verbose, useCpp = useCpp
    )
    data.table::setattr(out, "ReasonForStop", NULL)
    sort(paste(out$pixelIndex, out$speciesCode, sep = "/"))
  }

  for (useCpp in c(FALSE, TRUE)) {
    k0 <- successKeys(useCpp, verbose = 0)
    k1 <- successKeys(useCpp, verbose = 1)
    expect_equal(k1, k0,
                 info = sprintf("useCpp=%s verbose 0 vs 1", useCpp))
  }
  ## And the Success keys must agree across impls at each verbose level
  expect_equal(successKeys(useCpp = TRUE,  verbose = 0),
               successKeys(useCpp = FALSE, verbose = 0),
               info = "verbose=0  Cpp == R")
  expect_equal(successKeys(useCpp = TRUE,  verbose = 1),
               successKeys(useCpp = FALSE, verbose = 1),
               info = "verbose=1  Cpp == R")
})

# ---------------------------------------------------------------------------
# 12. Long-soak: 200 sequential LANDISDisp calls on the same fixture, both
# paths consuming the SAME RNG stream. Catches any cumulative RNG-state leak
# in the C++ port (Get/PutRNGstate misuse, accidental rng draws outside the
# loop, etc.). Gated by LANDR_SLOW_TESTS=1 because it is ~30 s.
# ---------------------------------------------------------------------------
test_that("Long soak: 200 sequential calls produce identical outputs (slow)", {
  if (!identical(Sys.getenv("LANDR_SLOW_TESTS"), "1")) {
    skip("set LANDR_SLOW_TESTS=1 to enable long-soak parity test")
  }
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "tiny", fixtureSeed = 11L)
  N <- 200L

  hashSeq <- function(useCpp) {
    set.seed(20260501L)
    vapply(seq_len(N), function(i) {
      out <- LANDISDisp(
        dtSrc = fix$dtSrc, dtRcv = fix$dtRcvFull,
        pixelGroupMap = fix$pixelGroupMap, speciesTable = fix$speciesTable,
        successionTimestep = fix$successionTimestep,
        verbose = 1, useCpp = useCpp
      )
      data.table::setattr(out, "ReasonForStop", NULL)
      cols <- intersect(c("pixelIndex", "speciesCode", "DistOfSuccess", "species"),
                        colnames(out))
      out <- out[, ..cols]
      data.table::setorderv(out, c("pixelIndex", "speciesCode"))
      digest::digest(out, algo = "sha256")
    }, character(1))
  }
  hR <- hashSeq(useCpp = FALSE)
  hC <- hashSeq(useCpp = TRUE)
  expect_equal(hC, hR)
  ## Sanity: the 200 outputs should not be all the same hash (RNG advancing)
  expect_gt(length(unique(hR)), N / 2)
})

# ---------------------------------------------------------------------------
# 13. Larger species table (16 species). Stays well within the 64-bit bitmask
# cap. Catches any accidental hard-coded numSp assumption in the C++ port.
# ---------------------------------------------------------------------------
test_that("R and Cpp agree with a 16-species fixture", {
  skip_if_no_cpp()

  for (rsd in c(1L, 42L, 9999L)) {
    for (ts in c(1L, 10L)) {
      fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L,
                                   successionTimestep = ts, nSpecies = 16L)
      info <- sprintf("16 species  ts=%d  rsd=%d", ts, rsd)
      parityCheck(fix, runSeed = rsd, info = info)
    }
  }
})

# ---------------------------------------------------------------------------
# 14. Sparse pixelGroup IDs in the millions. The C++ bitmask is a flat vector
# of length maxPg + 1; with maxPg ~ 5e6 that is ~40 MB. This test confirms it
# allocates and gives correct results, AND that the result matches the R
# reference (which uses srcPixelMatrix, a separate code path that also scales
# with maxPg).
# ---------------------------------------------------------------------------
test_that("R and Cpp agree with pixelGroup IDs in the millions", {
  skip_if_no_cpp()

  fix <- makeLANDISDispFixture(size = "small", fixtureSeed = 11L)
  oldPgs <- sort(unique(c(fix$dtSrc$pixelGroup, fix$dtRcv$pixelGroup,
                          as.vector(fix$pixelGroupMap[]))))
  oldPgs <- oldPgs[!is.na(oldPgs)]
  set.seed(13L)
  newIds <- sort(sample.int(5000000L, length(oldPgs)))
  remap <- setNames(newIds, as.character(oldPgs))

  pgmNew <- terra::rast(fix$pixelGroupMap)
  vv <- as.vector(fix$pixelGroupMap[])
  vNew <- remap[as.character(vv)]
  vNew[is.na(vNew)] <- NA_integer_
  terra::values(pgmNew) <- as.integer(vNew)
  fix$pixelGroupMap <- pgmNew
  fix$dtSrc[, pixelGroup := remap[as.character(pixelGroup)]]
  fix$dtRcv[, pixelGroup := remap[as.character(pixelGroup)]]
  fix$dtRcvFull[, pixelGroup := remap[as.character(pixelGroup)]]

  parityCheck(fix, runSeed = 42L,
              info = sprintf("maxPg=%d", max(newIds)))
})

# ---------------------------------------------------------------------------
# 15. Cross-implementation parity on the hash-manifest scenario grid.
# Each spec runs LANDISDisp twice in-session — once with useCpp=TRUE and once
# with useCpp=FALSE under the same seed — and asserts the outputs match.
# Replaces the previous saved-SHA256 manifest test, which was R-version
# brittle: data.table's internal index bytes differ between R minor versions
# even when content compares equal under expect_equal/all.equal.
# ---------------------------------------------------------------------------
test_that("Cpp and R implementations agree on hash-manifest scenario grid", {
  skip_if_no_cpp()

  specs <- list(
    list(size = "tiny",         fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "tiny",         fixtureSeed = 11L, runSeed = 1729L, ts =  1L),
    list(size = "tiny",         fixtureSeed = 23L, runSeed =   42L, ts =  1L),
    list(size = "small",        fixtureSeed = 11L, runSeed =   42L, ts = 10L),
    list(size = "small",        fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
    list(size = "medium",       fixtureSeed = 11L, runSeed =   42L, ts = 10L),
    list(size = "xlarge_dense", fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "xlarge_dense", fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
    list(size = "xxlarge",      fixtureSeed = 11L, runSeed =   42L, ts =  1L),
    list(size = "xxxlarge",     fixtureSeed = 11L, runSeed =   42L, ts =  1L)
  )

  slowSizes <- c("xlarge_dense", "xxlarge", "xxxlarge")
  isSlowEnabled <- identical(Sys.getenv("LANDR_SLOW_TESTS"), "1")

  for (s in specs) {
    if (s$size %in% slowSizes && !isSlowEnabled) next
    info <- sprintf("size=%s fixSeed=%d runSeed=%d ts=%d",
                    s$size, s$fixtureSeed, s$runSeed, s$ts)
    fix <- makeLANDISDispFixture(size = s$size, fixtureSeed = s$fixtureSeed,
                                 successionTimestep = s$ts)
    parityCheck(fix, runSeed = s$runSeed, info = info)
  }
})

