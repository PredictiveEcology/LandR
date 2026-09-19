#!/usr/bin/env Rscript
## Profile LANDISDisp on the xxxlarge fixture (9M cells, ~4% rcv) to confirm
## where wall time goes. Two angles:
##
##  1) Manual stage timing inside LANDISDisp by wrapping its key sections via
##     a temporary trace, so we can attribute time to:
##        - factor compaction
##        - srcPixelMatrix construction (dtSrcLong join + speciesSrcRasterVecList + cbind)
##        - spiralSeedDispersal{R,Cpp} prep (spiral, distsBySpCode/wardProb)
##        - the inner spiral loop itself (R loop vs spiralLoopCpp call)
##        - post-loop join + filtering
##
##  2) Rprof on each path, with line.profiling = TRUE, to corroborate (1).
##
## Run from package root:
##   Rscript .dev_landisdisp/profile-xxxlarge.R

suppressPackageStartupMessages({
  library(data.table); library(terra); library(SpaDES.tools); library(reproducible); library(Rcpp)
})

.pkgEnv <- new.env(parent = emptyenv())
sourceCpp("src/spiralSeedDispersal.cpp", verbose = FALSE, showOutput = FALSE)
source("R/seedDispersalLANDIS.R")
source("tests/testthat/helper-LANDISDisp-fixtures.R")

cat("Building xxxlarge fixture...\n")
t0 <- Sys.time()
fix <- makeLANDISDispFixture(size = "xxxlarge", fixtureSeed = 11L,
                             successionTimestep = 1L)
cat(sprintf("  fixture build: %.2fs  (cells=%d, rcv pgs=%d, src pgs=%d)\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            terra::ncell(fix$pixelGroupMap),
            uniqueN(fix$dtRcv$pixelGroup), uniqueN(fix$dtSrc$pixelGroup)))

## --- Stage timer ----------------------------------------------------------
## We instrument LANDISDisp by injecting Sys.time() probes into a copy of its
## body. Returns a named numeric of seconds for each stage.
stageTimes <- function(useCpp) {
  ## Build a stand-alone version of LANDISDisp that records per-stage timing
  ## by wrapping each block with Sys.time(). Easier than tracing because we
  ## want stable, well-named segments.
  timed <- function(dtSrc, dtRcv, pixelGroupMap, speciesTable,
                    successionTimestep, verbose = 1, useCpp) {
    times <- list()
    mark <- function(name) {
      times[[name]] <<- Sys.time(); invisible(NULL)
    }
    delta <- function(a, b) as.numeric(times[[b]] - times[[a]], units = "secs")

    mark("t0")

    ## ---- Factor / coercion block (mirrors LANDISDisp top) ----
    dtSrc <- data.table::copy(dtSrc)
    dtRcv <- data.table::copy(dtRcv)
    speciesTable <- data.table::copy(speciesTable)

    origClassWasNumeric <- is.numeric(speciesTable[["speciesCode"]])
    if (is(dtSrc$speciesCode, "numeric")) {
      data.table::set(dtSrc, NULL, c("speciesCode"),
                      factor(dtSrc[["speciesCode"]]))
      origLevels <- levels(dtSrc[["speciesCode"]])
      speciesTable <- speciesTable[as.numeric(origLevels), ]
      data.table::set(speciesTable, NULL, c("speciesCode"),
                      factor(speciesTable[["speciesCode"]], levels = origLevels))
      data.table::set(dtRcv, NULL, c("speciesCode"),
                      factor(dtRcv[["speciesCode"]], levels = origLevels))
    }
    if (is.factor(dtSrc$speciesCode)) {
      origLevels <- levels(dtSrc$speciesCode)
      dtSrc[, speciesCode2 := as.integer(speciesCode)]
      dtRcv[, speciesCode2 := as.integer(speciesCode)]
      speciesTable[, speciesCode2 := as.integer(speciesCode)]
      data.table::set(dtSrc, NULL, "speciesCode", NULL)
      data.table::set(dtRcv, NULL, "speciesCode", NULL)
      data.table::set(speciesTable, NULL, "speciesCode", NULL)
      data.table::setnames(dtSrc, "speciesCode2", "speciesCode")
      data.table::setnames(dtRcv, "speciesCode2", "speciesCode")
      data.table::setnames(speciesTable, "speciesCode2", "speciesCode")
      if (!"species" %in% colnames(speciesTable)) {
        data.table::set(speciesTable, NULL, "species",
                        paste0("Spp_", speciesTable[["speciesCode"]]))
      }
      data.table::setorderv(speciesTable, "speciesCode")
      data.table::setorderv(dtSrc, "speciesCode")
      data.table::setorderv(dtRcv, "speciesCode")
    }

    mark("t_factor")

    ## ---- srcPixelMatrix construction ----
    pgv <- as.vector(pixelGroupMap[])
    rasVectorTemplate <- rep(NA_integer_, terra::ncell(pixelGroupMap))
    srcSpeciesCodes <- sort(unique(dtSrc$speciesCode))
    names(srcSpeciesCodes) <- as.character(srcSpeciesCodes)
    cellsCanSrc <- which(pgv %in% dtSrc$pixelGroup)
    dtSrcLong <- data.table::data.table(pixelGroup = pgv[cellsCanSrc],
                                        pixelIndex = cellsCanSrc)
    dtSrcLong <- dtSrc[, c("pixelGroup", "speciesCode")][dtSrcLong, on = "pixelGroup",
                                                          allow.cartesian = TRUE]
    data.table::set(dtSrcLong, NULL, "pixelGroup", NULL)
    data.table::setkeyv(dtSrcLong, "speciesCode")

    srcSpeciesByIndex <- split(dtSrcLong, by = "speciesCode")
    speciesSrcRasterVecList <- lapply(srcSpeciesCodes, function(sc) {
      rasVectorTemplate[srcSpeciesByIndex[[as.character(sc)]][["pixelIndex"]]] <- sc
      rasVectorTemplate
    })
    maxSpCode <- max(as.integer(srcSpeciesCodes))
    speciesSrcRasterVecList <- lapply(seq_len(maxSpCode), function(ind) {
      if (as.character(ind) %in% names(speciesSrcRasterVecList)) {
        speciesSrcRasterVecList[[as.character(ind)]]
      }
    })
    srcPixelMatrix <- do.call(cbind, speciesSrcRasterVecList)

    mark("t_srcMatrix")

    ## ---- dtRcvLong setup ----
    cellSize <- unique(terra::res(pixelGroupMap))[1]
    dtRcvNew <- dtRcv[unique(dtSrc[, "speciesCode"], by = "speciesCode"),
                     on = "speciesCode", nomatch = NULL]
    cellsCanRcv <- which(pgv %in% dtRcvNew$pixelGroup)
    dtRcvLong <- data.table::data.table(pixelGroup = pgv[cellsCanRcv],
                                        pixelIndex = cellsCanRcv)
    dtRcvSmall <- dtRcvNew[, c("pixelGroup", "speciesCode")]
    dtRcvLong <- dtRcvLong[dtRcvSmall, on = "pixelGroup",
                           allow.cartesian = TRUE, nomatch = NULL]
    data.table::setorderv(dtRcvLong, c("pixelIndex", "speciesCode"))

    mark("t_rcvLong")

    ## ---- spiralSeedDispersal(R or Cpp) ----
    spiralFn <- if (useCpp) spiralSeedDispersalCpp else spiralSeedDispersalR
    out <- spiralFn(speciesTable, pixelGroupMap, dtRcvLong,
                    srcPixelMatrix, cellSize, k = 0.95, b = 0.01,
                    successionTimestep, verbose, dispersalFn = Ward)

    mark("t_dispersal")

    list(
      output = out,
      times = c(
        factor       = delta("t0", "t_factor"),
        srcMatrix    = delta("t_factor", "t_srcMatrix"),
        rcvLong      = delta("t_srcMatrix", "t_rcvLong"),
        dispersalAll = delta("t_rcvLong", "t_dispersal"),
        TOTAL        = delta("t0", "t_dispersal")
      )
    )
  }

  set.seed(42L)
  ## warm-up (cache spiral)
  invisible(timed(fix$dtSrc, fix$dtRcvFull, fix$pixelGroupMap,
                  fix$speciesTable, fix$successionTimestep, useCpp = useCpp))
  set.seed(42L)
  res <- timed(fix$dtSrc, fix$dtRcvFull, fix$pixelGroupMap,
               fix$speciesTable, fix$successionTimestep, useCpp = useCpp)
  res$times
}

cat("\n--- Stage timing: useCpp = TRUE ---\n")
ttCpp <- stageTimes(useCpp = TRUE)
print(round(ttCpp, 3))

cat("\n--- Stage timing: useCpp = FALSE ---\n")
ttR <- stageTimes(useCpp = FALSE)
print(round(ttR, 3))

cat("\n--- Comparison (Cpp vs R) ---\n")
df <- data.frame(
  stage = names(ttCpp),
  cpp_s = round(ttCpp, 3),
  r_s   = round(ttR, 3),
  delta = round(ttR - ttCpp, 3)
)
df$pct_cpp <- round(100 * df$cpp_s / df$cpp_s[df$stage == "TOTAL"], 1)
df$pct_r   <- round(100 * df$r_s   / df$r_s[df$stage == "TOTAL"], 1)
print(df, row.names = FALSE)

## --------------------------------------------------------------------------
## (2) Rprof on each path. line.profiling captures intra-function lines so we
## can see where inside spiralSeedDispersalR the R impl spends its time.
## --------------------------------------------------------------------------
profPath <- function(useCpp, file) {
  set.seed(42L)
  invisible(LANDISDisp(  ## warm-up
    dtSrc = fix$dtSrc, dtRcv = fix$dtRcvFull, pixelGroupMap = fix$pixelGroupMap,
    speciesTable = fix$speciesTable, successionTimestep = fix$successionTimestep,
    verbose = 1, useCpp = useCpp))
  set.seed(42L)
  Rprof(filename = file, line.profiling = TRUE, interval = 0.02)
  invisible(LANDISDisp(
    dtSrc = fix$dtSrc, dtRcv = fix$dtRcvFull, pixelGroupMap = fix$pixelGroupMap,
    speciesTable = fix$speciesTable, successionTimestep = fix$successionTimestep,
    verbose = 1, useCpp = useCpp))
  Rprof(NULL)
}

cat("\n=== Rprof: useCpp = TRUE ===\n")
profPath(TRUE, ".dev_landisdisp/prof_cpp.out")
print(summaryRprof(".dev_landisdisp/prof_cpp.out", lines = "show")$by.self[1:15, ])

cat("\n=== Rprof: useCpp = FALSE ===\n")
profPath(FALSE, ".dev_landisdisp/prof_r.out")
print(summaryRprof(".dev_landisdisp/prof_r.out", lines = "show")$by.self[1:20, ])

cat("\n[done] traces in .dev_landisdisp/prof_cpp.out and prof_r.out\n")
