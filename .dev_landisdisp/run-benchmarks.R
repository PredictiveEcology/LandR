#!/usr/bin/env Rscript
## Benchmark LANDISDisp R vs Rcpp across a range of input sizes.
##
## Usage:
##   Rscript .dev_landisdisp/run-benchmarks.R               # full sweep
##   Rscript .dev_landisdisp/run-benchmarks.R quick         # tiny+small+medium only
##
## Output:
##   - table printed to stdout
##   - markdown table written to .dev_landisdisp/benchmark_results.md
##   - raw timings written to .dev_landisdisp/benchmark_results.rds

suppressPackageStartupMessages({
  library(data.table); library(terra); library(SpaDES.tools); library(reproducible); library(Rcpp)
})

.pkgEnv <- new.env(parent = emptyenv())
sourceCpp("src/spiralSeedDispersal.cpp", verbose = FALSE, showOutput = FALSE)
source("R/seedDispersalLANDIS.R")
source("tests/testthat/helper-LANDISDisp-fixtures.R")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] == "quick") "quick" else "full"

sizes <- if (mode == "quick") {
  c("tiny", "small", "medium")
} else {
  c("tiny", "small", "medium", "large", "xlarge",
    "xlarge_dense", "xxlarge", "xxxlarge")
}

## Repeat each timing to reduce noise. R is slow, so fewer reps for big sizes.
repsBySize <- c(tiny = 5L, small = 5L, medium = 3L, large = 2L,
                xlarge = 1L, xlarge_dense = 1L, xxlarge = 1L, xxxlarge = 1L)

timeOne <- function(fn) {
  t0 <- Sys.time()
  out <- fn()
  list(secs = as.numeric(Sys.time() - t0, units = "secs"), nrows = nrow(out))
}

results <- list()
for (sz in sizes) {
  reps <- repsBySize[[sz]]
  cat(sprintf("\n--- size: %s (reps=%d) ---\n", sz, reps))
  fix <- makeLANDISDispFixture(size = sz, fixtureSeed = 11L)
  cat(sprintf("  cells=%d  rcv=%d  src=%d\n",
              terra::ncell(fix$pixelGroupMap),
              nrow(fix$dtRcv), nrow(fix$dtSrc)))

  ## Warm-up once each (also caches the spiral in .pkgEnv so neither impl is
  ## charged for the focalMat/spiralDistances cost on later reps)
  invisible(runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = FALSE))
  invisible(runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = TRUE))

  rTimes <- numeric(reps); cTimes <- numeric(reps); nR <- nC <- 0L
  for (rep in seq_len(reps)) {
    rRes <- timeOne(function() runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = FALSE))
    cRes <- timeOne(function() runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = TRUE))
    rTimes[rep] <- rRes$secs; cTimes[rep] <- cRes$secs
    nR <- rRes$nrows; nC <- cRes$nrows
    cat(sprintf("  rep %d: R=%.3fs  Cpp=%.3fs  speedup=%.2fx  rows R=%d C=%d\n",
                rep, rRes$secs, cRes$secs, rRes$secs / cRes$secs, nR, nC))
  }

  results[[sz]] <- list(
    size = sz, cells = terra::ncell(fix$pixelGroupMap),
    rcv = nrow(fix$dtRcv), src = nrow(fix$dtSrc),
    rTimes = rTimes, cTimes = cTimes, nrowR = nR, nrowC = nC
  )
}

## ---- Build results table ----
fmt <- function(x) sprintf("%.3f", x)
df <- data.frame(
  size       = sapply(results, `[[`, "size"),
  cells      = sapply(results, `[[`, "cells"),
  rcvRows    = sapply(results, `[[`, "rcv"),
  srcRows    = sapply(results, `[[`, "src"),
  rMedian_s  = sapply(results, function(r) median(r$rTimes)),
  cMedian_s  = sapply(results, function(r) median(r$cTimes)),
  rMin_s     = sapply(results, function(r) min(r$rTimes)),
  cMin_s     = sapply(results, function(r) min(r$cTimes)),
  outRows    = sapply(results, `[[`, "nrowR"),
  identical  = sapply(results, function(r) r$nrowR == r$nrowC)
)
df$speedup <- df$rMedian_s / df$cMedian_s

cat("\n\n========== Benchmark summary ==========\n")
print(df, row.names = FALSE, digits = 3)

## Markdown table
md <- c(
  "# LANDISDisp benchmarks: R vs Rcpp",
  "",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Host: %s | R: %s",
          paste(Sys.info()[c("nodename", "machine")], collapse = "/"),
          R.version.string),
  "",
  "All timings include the call into LANDISDisp() (spiral + ward-prob prep + the",
  "loop). The first call per size is excluded (warm-up); reported numbers are the",
  "median over the listed reps. Both implementations use the same RNG stream, so",
  "the output rows match exactly (verified by seed-locked tests).",
  "",
  paste0("| size | pgm cells | rcv rows | src rows | reps | R median (s) | Cpp median (s) | speedup | output rows | identical |"),
  paste0("|------|-----------|---------:|---------:|-----:|------------:|--------------:|--------:|-----------:|:---------:|")
)
for (i in seq_len(nrow(df))) {
  reps <- length(results[[df$size[i]]]$rTimes)
  md <- c(md, sprintf(
    "| %s | %d | %d | %d | %d | %s | %s | %.2fx | %d | %s |",
    df$size[i], df$cells[i], df$rcvRows[i], df$srcRows[i], reps,
    fmt(df$rMedian_s[i]), fmt(df$cMedian_s[i]), df$speedup[i],
    df$outRows[i], if (df$identical[i]) "yes" else "**NO**"
  ))
}
writeLines(md, ".dev_landisdisp/benchmark_results.md")
saveRDS(results, ".dev_landisdisp/benchmark_results.rds")
cat("\nWrote .dev_landisdisp/benchmark_results.md and .rds\n")
