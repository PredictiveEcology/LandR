#!/usr/bin/env Rscript
## Build seed-locked golden baselines for LANDISDisp from the *current* (R-only)
## implementation. Re-run only when intentionally rebaselining; tests compare
## against these RDS snapshots.
##
## Outputs (committed to fixture dir for tests):
##   tests/testthat/fixtures/LANDISDisp_<size>_fix<fixSeed>_run<runSeed>.rds
##
## Run from package root:
##   Rscript .dev_landisdisp/build-baselines.R

suppressPackageStartupMessages({
  library(data.table); library(terra); library(SpaDES.tools); library(reproducible); library(digest)
})

.pkgEnv <- new.env(parent = emptyenv())
source("R/seedDispersalLANDIS.R")
source("tests/testthat/helper-LANDISDisp-fixtures.R")

fixturesDir <- "tests/testthat/fixtures"
dir.create(fixturesDir, showWarnings = FALSE, recursive = TRUE)

## (size, fixtureSeed, runSeed, successionTimestep) tuples covering both the
## ts==1 and ts>1 code branches, while keeping per-species saturation < 60%
## of eligible receivers (so the goldens have real discriminative power).
specs <- list(
  list(size = "tiny",         fixtureSeed = 11L, runSeed =   42L, ts =  1L),
  list(size = "tiny",         fixtureSeed = 11L, runSeed = 1729L, ts =  1L),
  list(size = "tiny",         fixtureSeed = 23L, runSeed =   42L, ts =  1L),
  list(size = "small",        fixtureSeed = 11L, runSeed =   42L, ts = 10L),
  list(size = "small",        fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
  list(size = "medium",       fixtureSeed = 11L, runSeed =   42L, ts = 10L),
  ## rcv-heavy x-large fixtures — nominally slow (regenerate ~30s total)
  list(size = "xlarge_dense", fixtureSeed = 11L, runSeed =   42L, ts =  1L),
  list(size = "xlarge_dense", fixtureSeed = 11L, runSeed = 1729L, ts = 10L),
  list(size = "xxlarge",      fixtureSeed = 11L, runSeed =   42L, ts =  1L),
  ## landscape-scale fixture (~9M cells, ~4% rcv) — adds ~25s to regeneration
  list(size = "xxxlarge",     fixtureSeed = 11L, runSeed =   42L, ts =  1L)
)

manifest <- data.frame(
  file = character(0), sha256 = character(0),
  rows = integer(0), generated_at = character(0),
  stringsAsFactors = FALSE
)
for (s in specs) {
  cat(sprintf("[baseline] size=%-13s fixSeed=%d runSeed=%d ts=%2d ... ",
              s$size, s$fixtureSeed, s$runSeed, s$ts))
  fix <- makeLANDISDispFixture(size = s$size, fixtureSeed = s$fixtureSeed,
                               successionTimestep = s$ts)
  ## Goldens are anchored to the R reference implementation; Cpp must match.
  t0 <- Sys.time()
  out <- runLANDISDispOnFixture(fix, runSeed = s$runSeed, useCpp = FALSE)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  fname <- sprintf("LANDISDisp_%s_fix%d_run%d_ts%d.rds",
                   s$size, s$fixtureSeed, s$runSeed, s$ts)
  saveRDS(out, file.path(fixturesDir, fname), version = 2L)
  ## Hash sentinel: thumbprint that travels with the goldens. Anyone reviewing
  ## a regenerate-baselines diff can see at a glance which goldens changed.
  h <- digest::digest(out, algo = "sha256")
  manifest <- rbind(manifest, data.frame(
    file = fname, sha256 = h, rows = nrow(out),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  ))
  cat(sprintf("rows=%d  %.2fs  -> %s  sha256=%s...\n",
              nrow(out), dt, fname, substr(h, 1, 12)))
}

manifestPath <- file.path(fixturesDir, "MANIFEST.csv")
write.csv(manifest, manifestPath, row.names = FALSE, quote = TRUE)
cat(sprintf("\nWrote manifest: %s (%d entries)\n", manifestPath, nrow(manifest)))
cat("[baseline] done. Snapshots in", fixturesDir, "\n")
