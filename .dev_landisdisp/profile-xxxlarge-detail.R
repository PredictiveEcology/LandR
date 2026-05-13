#!/usr/bin/env Rscript
## Finer-grained profile of LANDISDisp on xxxlarge: break the two big stages
## (srcMatrix, dispersalAll) into their internal sub-steps.

suppressPackageStartupMessages({
  library(data.table); library(terra); library(SpaDES.tools); library(reproducible); library(Rcpp)
})

.pkgEnv <- new.env(parent = emptyenv())
sourceCpp("src/spiralSeedDispersal.cpp", verbose = FALSE, showOutput = FALSE)
source("R/seedDispersalLANDIS.R")
source("tests/testthat/helper-LANDISDisp-fixtures.R")

cat("Building xxxlarge fixture...\n")
fix <- makeLANDISDispFixture(size = "xxxlarge", fixtureSeed = 11L,
                             successionTimestep = 1L)
cat(sprintf("  cells=%d, rcv pgs=%d, src pgs=%d, srcRows=%d\n",
            terra::ncell(fix$pixelGroupMap),
            uniqueN(fix$dtRcv$pixelGroup), uniqueN(fix$dtSrc$pixelGroup),
            nrow(fix$dtSrc)))

## ---- replicate srcPixelMatrix construction with sub-stage timing --------
buildSrcMatrix <- function(dtSrc, dtRcv, speciesTable, pixelGroupMap) {
  ## (post-factor-compaction; same logic as inside LANDISDisp)
  times <- list()
  mark <- function(name) { times[[name]] <<- Sys.time(); invisible(NULL) }
  delta <- function(a, b) as.numeric(times[[b]] - times[[a]], units = "secs")
  mark("t0")

  pgv <- as.vector(pixelGroupMap[])
  mark("t_pgv")

  rasVectorTemplate <- rep(NA_integer_, terra::ncell(pixelGroupMap))
  mark("t_template")

  srcSpeciesCodes <- sort(unique(dtSrc$speciesCode))
  names(srcSpeciesCodes) <- as.character(srcSpeciesCodes)
  cellsCanSrc <- which(pgv %in% dtSrc$pixelGroup)
  mark("t_which")

  dtSrcLong <- data.table(pixelGroup = pgv[cellsCanSrc], pixelIndex = cellsCanSrc)
  dtSrcLong <- dtSrc[, c("pixelGroup", "speciesCode")][dtSrcLong, on = "pixelGroup",
                                                         allow.cartesian = TRUE]
  set(dtSrcLong, NULL, "pixelGroup", NULL)
  setkeyv(dtSrcLong, "speciesCode")
  mark("t_srcLong")

  srcSpeciesByIndex <- split(dtSrcLong, by = "speciesCode")
  mark("t_split")

  speciesSrcRasterVecList <- lapply(srcSpeciesCodes, function(sc) {
    rasVectorTemplate[srcSpeciesByIndex[[as.character(sc)]][["pixelIndex"]]] <- sc
    rasVectorTemplate
  })
  mark("t_lapply")

  maxSpCode <- max(as.integer(srcSpeciesCodes))
  speciesSrcRasterVecList <- lapply(seq_len(maxSpCode), function(ind) {
    if (as.character(ind) %in% names(speciesSrcRasterVecList)) {
      speciesSrcRasterVecList[[as.character(ind)]]
    }
  })
  mark("t_compact")

  srcPixelMatrix <- do.call(cbind, speciesSrcRasterVecList)
  mark("t_cbind")

  list(
    matrix = srcPixelMatrix,
    times = c(
      pgv         = delta("t0", "t_pgv"),
      template    = delta("t_pgv", "t_template"),
      cellsCanSrc = delta("t_template", "t_which"),
      srcLongJoin = delta("t_which", "t_srcLong"),
      split       = delta("t_srcLong", "t_split"),
      writeMatrix = delta("t_split", "t_lapply"),
      compactList = delta("t_lapply", "t_compact"),
      cbind       = delta("t_compact", "t_cbind"),
      TOTAL       = delta("t0", "t_cbind")
    )
  )
}

## ---- replicate spiralSeedDispersalCpp prep with sub-stage timing -------
profDispersalCpp <- function(speciesTable, pixelGroupMap, dtRcvLong,
                             pgv_, dtSrc_, cellSize, k, b,
                             successionTimestep, verbose, dispersalFn) {
  times <- list()
  mark <- function(name) { times[[name]] <<- Sys.time(); invisible(NULL) }
  delta <- function(a, b) as.numeric(times[[b]] - times[[a]], units = "secs")
  mark("t0")

  speciesTable <- copy(speciesTable)
  set(speciesTable, NULL, "seeddistance_maxMinCellSize",
      pmax(cellSize, speciesTable[["seeddistance_max"]]))
  maxDis <- max(speciesTable[, "seeddistance_maxMinCellSize"])

  preExistingSpiral <- paste0("spirals_max", round(maxDis, 6),
                              "_cell", round(cellSize, 6))
  if (!exists(preExistingSpiral, envir = .pkgEnv)) {
    .pkgEnv[[preExistingSpiral]] <- spiralDistances(pixelGroupMap, maxDis, cellSize)
  }
  spiral <- .pkgEnv[[preExistingSpiral]]
  mark("t_spiral")

  speciesTableSmall <- speciesTable[, c("speciesCode", "seeddistance_eff", "seeddistance_max")]
  uniqueDists <- unique(spiral[, "dists", drop = FALSE]) * cellSize
  numSp <- NROW(speciesTable)
  spSeq <- seq(numSp)
  distsBySpCode <- as.data.table(expand.grid(dists = uniqueDists,
                                              speciesCode = speciesTable[["speciesCode"]]))
  set(distsBySpCode, NULL, "seeddistance_max", speciesTableSmall[
    distsBySpCode[["speciesCode"]], "seeddistance_max"])
  set(distsBySpCode, NULL, "seeddistance_eff", speciesTableSmall[
    distsBySpCode[["speciesCode"]], "seeddistance_eff"])
  set(distsBySpCode, NULL, "wardProb",
      pmin(1, dispersalFn(dist = distsBySpCode$dists, cellSize = cellSize,
                          effDist = distsBySpCode$seeddistance_eff,
                          maxDist = distsBySpCode$seeddistance_max, k = k, b = b)))
  set(distsBySpCode, NULL, c("seeddistance_max", "seeddistance_eff"), NULL)
  setorderv(distsBySpCode, c("dists", "speciesCode"))
  if (successionTimestep > 1) {
    set(distsBySpCode, NULL, "wardProb",
        1 - (1 - distsBySpCode[["wardProb"]])^successionTimestep)
  }
  numUniqueDists <- length(uniqueDists)
  wardProbByDist <- matrix(distsBySpCode[["wardProb"]],
                           nrow = numUniqueDists, ncol = numSp, byrow = TRUE)
  mark("t_wardprob")

  rcvFull <- dtRcvLong[, c("pixelIndex", "speciesCode")]
  rcvFull <- rcvFull[speciesTable[, c("seeddistance_max", "speciesCode")],
                     on = "speciesCode", nomatch = NULL]
  rc1 <- rowColFromCell(pixelGroupMap, rcvFull[["pixelIndex"]])
  colnames(rc1) <- c("row", "col")
  rowOrig <- as.integer(rc1[, "row"])
  colOrig <- as.integer(rc1[, "col"])
  curDists <- drop(spiral[, 3]) * cellSize
  spiralRow <- as.integer(spiral[, "row"])
  spiralCol <- as.integer(spiral[, "col"])
  activeSpMaxDist <- numeric(numSp + 1L); activeSpMax <- numeric(numSp + 1L)
  for (rowIdx in seq_len(NROW(speciesTable))) {
    sc <- speciesTable[["speciesCode"]][rowIdx]
    activeSpMaxDist[sc + 1L] <- speciesTable[["seeddistance_maxMinCellSize"]][rowIdx]
    activeSpMax[sc + 1L]     <- speciesTable[["seeddistance_max"]][rowIdx]
  }
  mark("t_rcvprep")

  res <- spiralLoopCpp(
    pixelIndex_in       = as.integer(rcvFull[["pixelIndex"]]),
    speciesCode_in      = as.integer(rcvFull[["speciesCode"]]),
    rowOrig_in          = rowOrig, colOrig_in = colOrig,
    seeddist_max_perRow = as.integer(rcvFull[["seeddistance_max"]]),
    spiralRow = spiralRow, spiralCol = spiralCol,
    spiralCurDist = curDists,
    pgmRows = nrow(pixelGroupMap), pgmCols = ncol(pixelGroupMap),
    pgv = as.integer(pgv_),
    srcPg = as.integer(dtSrc_[["pixelGroup"]]),
    srcSpeciesCode = as.integer(dtSrc_[["speciesCode"]]),
    numSp = as.integer(numSp),
    wardProbByDist = wardProbByDist,
    activeSpMaxDist = activeSpMaxDist, activeSpMax = activeSpMax,
    cellSize = cellSize,
    successionTimestep = as.integer(successionTimestep),
    verbose = as.integer(verbose), wardAlreadyExp = TRUE, debug = FALSE
  )
  mark("t_loop")

  whSuccess <- which(res$Success)
  if (verbose >= 1) {
    set(rcvFull, NULL, "DistOfSuccess", res$DistOfSuccess)
    fails <- which(is.na(rcvFull[["DistOfSuccess"]]))
    if (length(fails)) {
      set(rcvFull, NULL, "ReasonForStop", NA_character_)
      set(rcvFull, fails, "ReasonForStop", "RanOutOfDistance")
    }
  }
  set(rcvFull, NULL, "seeddistance_max", NULL)
  if (length(whSuccess) == 0L) rcvFull <- rcvFull[0] else rcvFull <- rcvFull[whSuccess]
  speciesCodeCols <- intersect(c("species", "speciesCode"), colnames(speciesTable))
  rcvFull <- rcvFull[speciesTable[, ..speciesCodeCols], on = "speciesCode", nomatch = NULL]
  mark("t_post")

  list(
    output = rcvFull,
    times = c(
      spiralBuild = delta("t0", "t_spiral"),
      wardProbTab = delta("t_spiral", "t_wardprob"),
      rcvPrep     = delta("t_wardprob", "t_rcvprep"),
      cppLoop     = delta("t_rcvprep", "t_loop"),
      postProc    = delta("t_loop", "t_post"),
      TOTAL       = delta("t0", "t_post")
    )
  )
}

## ---- run ----------------------------------------------------------------
## Replicate the LANDISDisp pre-loop work to get the inputs we need
dtSrc <- copy(fix$dtSrc); dtRcv <- copy(fix$dtRcvFull); speciesTable <- copy(fix$speciesTable)

set(dtSrc, NULL, "speciesCode", factor(dtSrc[["speciesCode"]]))
origLevels <- levels(dtSrc[["speciesCode"]])
speciesTable <- speciesTable[as.numeric(origLevels), ]
set(speciesTable, NULL, "speciesCode", factor(speciesTable[["speciesCode"]], levels = origLevels))
set(dtRcv, NULL, "speciesCode", factor(dtRcv[["speciesCode"]], levels = origLevels))
dtSrc[, speciesCode2 := as.integer(speciesCode)]
dtRcv[, speciesCode2 := as.integer(speciesCode)]
speciesTable[, speciesCode2 := as.integer(speciesCode)]
set(dtSrc, NULL, "speciesCode", NULL); set(dtRcv, NULL, "speciesCode", NULL)
set(speciesTable, NULL, "speciesCode", NULL)
setnames(dtSrc, "speciesCode2", "speciesCode")
setnames(dtRcv, "speciesCode2", "speciesCode")
setnames(speciesTable, "speciesCode2", "speciesCode")

## warm-up
invisible(buildSrcMatrix(dtSrc, dtRcv, speciesTable, fix$pixelGroupMap))
res1 <- buildSrcMatrix(dtSrc, dtRcv, speciesTable, fix$pixelGroupMap)
cat("\n--- srcPixelMatrix sub-stages ---\n")
print(round(res1$times, 3))
srcPixelMatrix <- res1$matrix

## build dtRcvLong
pgv <- as.vector(fix$pixelGroupMap[])
dtRcvNew <- dtRcv[unique(dtSrc[, "speciesCode"], by = "speciesCode"),
                  on = "speciesCode", nomatch = NULL]
cellsCanRcv <- which(pgv %in% dtRcvNew$pixelGroup)
dtRcvLong <- data.table(pixelGroup = pgv[cellsCanRcv], pixelIndex = cellsCanRcv)
dtRcvLong <- dtRcvLong[dtRcvNew[, c("pixelGroup", "speciesCode")],
                       on = "pixelGroup", allow.cartesian = TRUE, nomatch = NULL]
setorderv(dtRcvLong, c("pixelIndex", "speciesCode"))
cat(sprintf("dtRcvLong rows: %d\n", nrow(dtRcvLong)))

cellSize <- terra::res(fix$pixelGroupMap)[1]
## warm-up
invisible(profDispersalCpp(speciesTable, fix$pixelGroupMap, dtRcvLong,
                           pgv_ = pgv, dtSrc_ = dtSrc,
                           cellSize, k = 0.95, b = 0.01,
                           successionTimestep = 1, verbose = 1, dispersalFn = Ward))
set.seed(42L)
res2 <- profDispersalCpp(speciesTable, fix$pixelGroupMap, dtRcvLong,
                         pgv_ = pgv, dtSrc_ = dtSrc,
                         cellSize, k = 0.95, b = 0.01,
                         successionTimestep = 1, verbose = 1, dispersalFn = Ward)
cat("\n--- spiralSeedDispersalCpp sub-stages ---\n")
print(round(res2$times, 3))
