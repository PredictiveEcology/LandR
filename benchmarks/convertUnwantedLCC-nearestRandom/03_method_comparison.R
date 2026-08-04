## Real-landscape head-to-head for convertUnwantedLCC(method=): the pre-1.2.0.9004 spiral
## (iterative spread2, stochastic) vs method = "nearest" (deterministic) vs the new
## method = "nearestRandom" (stochastic, weighted by local abundance).
##
## The question a reviewer needs answered is NOT "does nearestRandom match the old output
## pixel-for-pixel" -- it cannot, and neither could two runs of the old algorithm, which was
## itself random. It is "does nearestRandom put the *same mix* of land-cover classes on the
## ground as the old algorithm did, where 'nearest' does not". So alongside per-pixel
## agreement we report the composition of the assigned classes, and the total-variation
## distance between each method's composition and the old algorithm's. The old algorithm's
## own seed-to-seed values are the noise floor every other column must be read against.
##
## Inputs are the same four real SCANFI + FAO landscapes used by the previous PR's
## benchmark bundle (class 240 = FAO-forest pixels that are not a forest LCC class), plus
## the real Western-Alberta-Upland study area for the pathological case. They are private
## LandWeb data, so point at them with:
##   LCC_BENCH_DIR=<dir with v2_input_*.tif> WAU_LCC=<path to WAU rstLCC.tif> \
##     LANDR_SRC=. Rscript benchmarks/convertUnwantedLCC-nearestRandom/03_method_comparison.R

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  ## load_all(), not library(LandR): "nearestRandom" needs this branch's compiled
  ## windowCountsByClassCpp(), which no installed LandR has
  pkgload::load_all(Sys.getenv("LANDR_SRC", "."), quiet = TRUE)
})
terraOptions(memfrac = 0.2)

BENCH <- Sys.getenv("LCC_BENCH_DIR", "")
WAU <- Sys.getenv("WAU_LCC", "")
OUT <- file.path(Sys.getenv("LANDR_SRC", "."), "benchmarks", "convertUnwantedLCC-nearestRandom")
UNW <- 240L
EXCL <- c(0L, 20L, 30L) # not available as replacements (remapDT sends these to NA)
SEEDS <- c(123L, 987L, 5150L)
paddedFloatToChar <- reproducible::paddedFloatToChar
.resample <- function(x, ...) x[sample.int(length(x), ...)] # LandR-internal helper, inlined

cuNew <- convertUnwantedLCC

## the previous (pre-1.2.0.9004) implementation: iterative spread2, random pick among all
## valid cells found at the radius where the first one appeared
cuSpiral <- function(classesToReplace, rstLCC, availableERC_by_Sp, cap = Inf) {
  a <- data.table::copy(availableERC_by_Sp)
  if (!"speciesCode" %in% names(a)) a[, speciesCode := "allSpecies"]
  if (!"pixelIndex" %in% names(a)) a[, pixelIndex := seq(ncell(rstLCC))]
  unw <- sort(unique(a[
    gsub(".*_", "", initialEcoregionCode) %in% as.character(classesToReplace), pixelIndex
  ]))
  ERG2 <- unique(
    a[!gsub(".*_", "", initialEcoregionCode) %in% classesToReplace],
    by = c("speciesCode", "initialEcoregionCode")
  )
  it <- 1L
  cur <- length(unw)
  rep0 <- 0L
  out3 <- NULL
  t0 <- Sys.time()
  while (length(unw) > 0) {
    if (as.numeric(Sys.time() - t0, units = "secs") > cap) {
      return(list(out = unique(out3), finished = FALSE, remaining = length(unw)))
    }
    o <- SpaDES.tools::spread2(
      rstLCC, start = unw, asRaster = FALSE, iterations = it, allowOverlap = TRUE, spreadProb = 1
    )
    o <- o[initialPixels != pixels]
    it <- it + 1L
    o[, lcc := as.vector(rstLCC[])[pixels]][lcc %in% classesToReplace, lcc := NA]
    o <- na.omit(o)
    o5 <- a[o[, state := NULL],
      allow.cartesian = TRUE, on = c("pixelIndex" = "initialPixels"), nomatch = NA
    ]
    o5[, possERC := lcc]
    o6 <- na.omit(o5[ERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode"), nomatch = NA])
    rm0 <- o5[!ERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode")]
    o6 <- o6[!possERC %in% unique(rm0$possERC)]
    if (cur == length(unw)) rep0 <- rep0 + 1L else { cur <- length(unw); rep0 <- 0L }
    if (rep0 > 5) unw <- integer()
    if (nrow(o6) > 0) {
      keep <- o6[, list(k = .resample(.I, 1)), by = pixelIndex]
      o2 <- o6[keep$k][, list(pixelIndex, ecoregionGroup = as.integer(lcc))]
      unw <- unw[!unw %in% o2$pixelIndex]
      out3 <- if (is.null(out3)) o2 else rbindlist(list(o2, out3))
    }
  }
  list(out = unique(out3), finished = TRUE, remaining = 0L)
}

availDT <- function(r) {
  v <- values(r)[, 1]
  keep <- which(!is.na(v) & !v %in% EXCL)
  data.table(pixelIndex = keep, initialEcoregionCode = as.integer(v[keep]))
}

## assigned-class composition over the unwanted pixels, as proportions
composition <- function(out, classes) {
  tab <- table(factor(as.integer(out$ecoregionGroup), levels = classes))
  as.numeric(tab) / max(1, sum(tab))
}
## total variation distance between two compositions: 0 = identical mix, 1 = disjoint
tvd <- function(p, q) 0.5 * sum(abs(p - q))

agree <- function(a, b) {
  m <- merge(a[!is.na(ecoregionGroup)], b[!is.na(ecoregionGroup)], by = "pixelIndex")
  100 * mean(as.integer(m$ecoregionGroup.x) == as.integer(m$ecoregionGroup.y))
}

timeIt <- function(expr) {
  t <- system.time(val <- force(expr))[["elapsed"]]
  list(t = t, val = val)
}

rows <- list()
comps <- list()
labs <- c("small", "medium", "large", "bigblob")

for (lab in labs) {
  f <- file.path(BENCH, paste0("v2_input_", lab, ".tif"))
  if (!file.exists(f)) {
    message("skipping ", lab, ": ", f, " not found")
    next
  }
  L <- rast(f)
  v <- values(L)[, 1]
  aDT <- availDT(L)
  nUnw <- sum(v == UNW, na.rm = TRUE)
  classes <- sort(unique(v[!is.na(v) & v != UNW & !v %in% EXCL]))

  ## old algorithm, once per seed
  oldRuns <- lapply(SEEDS, function(s) {
    set.seed(s)
    r <- timeIt(cuSpiral(UNW, L, aDT))
    list(t = r$t, out = r$val$out, finished = r$val$finished)
  })
  oldOK <- vapply(oldRuns, `[[`, logical(1), "finished")

  ## new: deterministic, then stochastic once per seed
  det <- timeIt(suppressMessages(cuNew(UNW, L, copy(aDT), doAssertion = FALSE)))
  rnd <- lapply(SEEDS, function(s) {
    set.seed(s)
    timeIt(suppressMessages(cuNew(UNW, L, copy(aDT), doAssertion = FALSE, method = "nearestRandom")))
  })

  ## compositions (old averaged over seeds; it is the reference mix)
  cOldEach <- lapply(oldRuns[oldOK], function(o) composition(o$out, classes))
  cOld <- Reduce(`+`, cOldEach) / length(cOldEach)
  cDet <- composition(det$val, classes)
  cRndEach <- lapply(rnd, function(r) composition(r$val, classes))
  cRnd <- Reduce(`+`, cRndEach) / length(cRndEach)

  ## noise floors: what the old algorithm's own seed-to-seed variation looks like
  oldPairs <- utils::combn(which(oldOK), 2, simplify = FALSE)
  selfAgree <- mean(vapply(oldPairs, function(ij) {
    agree(oldRuns[[ij[1]]]$out, oldRuns[[ij[2]]]$out)
  }, numeric(1)))
  selfTVD <- mean(vapply(oldPairs, function(ij) {
    tvd(cOldEach[[which(which(oldOK) == ij[1])]], cOldEach[[which(which(oldOK) == ij[2])]])
  }, numeric(1)))

  rows[[lab]] <- data.table(
    landscape = lab,
    ncell = ncell(L),
    unwanted = nUnw,
    t_old = round(mean(vapply(oldRuns, `[[`, numeric(1), "t")), 2),
    t_nearest = round(det$t, 2),
    t_nearestRandom = round(mean(vapply(rnd, `[[`, numeric(1), "t")), 2),
    agree_old_vs_old = round(selfAgree, 1),
    agree_nearest_vs_old = round(mean(vapply(oldRuns[oldOK], function(o) agree(det$val, o$out), numeric(1))), 1),
    agree_rand_vs_old = round(mean(mapply(function(r, o) agree(r$val, o$out), rnd[seq_len(sum(oldOK))], oldRuns[oldOK])), 1),
    tvd_old_vs_old = round(selfTVD, 4),
    tvd_nearest_vs_old = round(tvd(cDet, cOld), 4),
    tvd_rand_vs_old = round(tvd(cRnd, cOld), 4)
  )
  comps[[lab]] <- data.table(
    landscape = lab, class = classes, old = round(cOld, 4),
    nearest = round(cDet, 4), nearestRandom = round(cRnd, 4)
  )
  cat(sprintf(
    ">>> %-8s ncell=%8d unw=%5d | t: old=%6.2fs near=%5.2fs rand=%5.2fs | agree vs old: self=%.1f%% near=%.1f%% rand=%.1f%% | TVD vs old: self=%.4f near=%.4f rand=%.4f\n",
    lab, ncell(L), nUnw, rows[[lab]]$t_old, rows[[lab]]$t_nearest, rows[[lab]]$t_nearestRandom,
    rows[[lab]]$agree_old_vs_old, rows[[lab]]$agree_nearest_vs_old, rows[[lab]]$agree_rand_vs_old,
    rows[[lab]]$tvd_old_vs_old, rows[[lab]]$tvd_nearest_vs_old, rows[[lab]]$tvd_rand_vs_old
  ))
  flush.console()
  rm(L)
  gc(verbose = FALSE)
}

tab <- rbindlist(rows)
print(as.data.frame(tab), row.names = FALSE)
fwrite(tab, file.path(OUT, "real_landscapes_methods.csv"))
fwrite(rbindlist(comps), file.path(OUT, "assigned_class_composition.csv"))

## ---- the pathological case: does nearestRandom stay fast where the old one never finished?
if (nzchar(WAU) && file.exists(WAU)) {
  L <- rast(WAU)
  v <- values(L)[, 1]
  if (sum(v == UNW, na.rm = TRUE) > 0) {
    aDT <- availDT(L)
    det <- timeIt(suppressMessages(cuNew(UNW, L, copy(aDT), doAssertion = FALSE)))
    set.seed(123)
    rnd <- timeIt(suppressMessages(cuNew(UNW, L, copy(aDT), doAssertion = FALSE, method = "nearestRandom")))
    wau <- data.table(
      study_area = "Western Alberta Upland", ncell = ncell(L),
      unwanted = sum(v == UNW, na.rm = TRUE),
      t_nearest = round(det$t, 2), t_nearestRandom = round(rnd$t, 2),
      agree_rand_vs_nearest = round(agree(rnd$val, det$val), 1)
    )
    print(as.data.frame(wau), row.names = FALSE)
    fwrite(wau, file.path(OUT, "wau_methods.csv"))
  } else {
    message("WAU raster has no class ", UNW, " pixels; skipping")
  }
}
