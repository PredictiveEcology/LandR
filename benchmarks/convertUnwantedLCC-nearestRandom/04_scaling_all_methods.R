## Does restoring the stochastic behaviour restore the blow-up? Self-contained scaling
## sweep: a single unwanted blob of increasing radius, timed under the pre-1.2.0.9004
## spiral (iterative spread2), method = "nearest", and method = "nearestRandom".
##
## The old implementation's cost grew with the square of the blob radius because it
## re-spread outward from every remaining unwanted pixel once per pass. "nearestRandom"
## must not reintroduce that: it adds a second set of distance transforms and one
## summed-area table per candidate class, all of which are O(ncell) and independent of
## how deep the blob is.
##
## Run: LANDR_SRC=. Rscript benchmarks/convertUnwantedLCC-nearestRandom/04_scaling_all_methods.R

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  pkgload::load_all(Sys.getenv("LANDR_SRC", "."), quiet = TRUE)
})
OUT <- file.path(Sys.getenv("LANDR_SRC", "."), "benchmarks", "convertUnwantedLCC-nearestRandom")
CAP <- as.numeric(Sys.getenv("BENCH_CAP", "120")) # wall-clock cap for the old implementation
paddedFloatToChar <- reproducible::paddedFloatToChar
.resample <- function(x, ...) x[sample.int(length(x), ...)]

cuSpiral <- function(classesToReplace, rstLCC, availableERC_by_Sp, cap = Inf) {
  a <- data.table::copy(availableERC_by_Sp)
  if (!"speciesCode" %in% names(a)) a[, speciesCode := "allSpecies"]
  if (!"pixelIndex" %in% names(a)) a[, pixelIndex := seq(ncell(rstLCC))]
  unw <- sort(unique(a[initialEcoregionCode %in% classesToReplace, pixelIndex]))
  ERG2 <- unique(
    a[!initialEcoregionCode %in% classesToReplace],
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

blob <- function(radius) {
  n <- as.integer(radius * 2.6)
  r <- rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "EPSG:3978")
  xy <- xyFromCell(r, 1:ncell(r))
  v <- ifelse(xy[, 1] < n / 2, 210L, 220L)
  d <- sqrt((xy[, 1] - n / 2)^2 + (xy[, 2] - n / 2)^2)
  v[d < radius] <- 240L
  setValues(r, v)
}

radii <- as.integer(strsplit(Sys.getenv("BENCH_RADII", "10,20,40,80,160"), ",")[[1]])
res <- rbindlist(lapply(radii, function(rad) {
  r <- blob(rad)
  v <- values(r)[, 1]
  aDT <- data.table(pixelIndex = seq_len(ncell(r)), initialEcoregionCode = as.integer(v))

  set.seed(1)
  tOld <- system.time(oOld <- suppressMessages(cuSpiral(240L, r, aDT, cap = CAP)))[["elapsed"]]
  tDet <- system.time(suppressMessages(
    convertUnwantedLCC(240L, r, copy(aDT), doAssertion = FALSE)
  ))[["elapsed"]]
  set.seed(1)
  tRnd <- system.time(suppressMessages(
    convertUnwantedLCC(240L, r, copy(aDT), doAssertion = FALSE, method = "nearestRandom")
  ))[["elapsed"]]

  out <- data.table(
    radius = rad, ncell = ncell(r), unwanted = sum(v == 240L),
    spiral_s = if (oOld$finished) round(tOld, 2) else NA_real_,
    spiral_done = oOld$finished,
    nearest_s = round(tDet, 2),
    nearestRandom_s = round(tRnd, 2),
    speedup_vs_spiral = if (oOld$finished) round(tOld / max(tRnd, 1e-6)) else NA_integer_
  )
  cat(sprintf(
    ">>> radius=%3d ncell=%8d unw=%7d | spiral=%s nearest=%5.2fs nearestRandom=%5.2fs\n",
    rad, ncell(r), sum(v == 240L),
    if (oOld$finished) sprintf("%7.2fs", tOld) else sprintf("DNF(>%gs, %d left)", CAP, oOld$remaining),
    tDet, tRnd
  ))
  flush.console()
  out
}))
print(as.data.frame(res), row.names = FALSE)
fwrite(res, file.path(OUT, "scaling_all_methods.csv"))

png(file.path(OUT, "fig1_scaling_all_methods.png"), width = 1100, height = 470, res = 110)
par(mar = c(4, 4, 3, 1))
ylim <- range(c(res$spiral_s, res$nearest_s, res$nearestRandom_s), na.rm = TRUE)
plot(res$radius, res$spiral_s,
  type = "b", log = "xy", pch = 19, col = "firebrick", ylim = ylim,
  xlab = "unwanted blob radius (cells)", ylab = "time (s, log)",
  main = "convertUnwantedLCC(): cost vs depth of the unwanted blob"
)
lines(res$radius, res$nearest_s, type = "b", pch = 19, col = "steelblue")
lines(res$radius, res$nearestRandom_s, type = "b", pch = 17, col = "darkgreen", lty = 2)
dnf <- res[spiral_done == FALSE]
if (nrow(dnf)) {
  points(dnf$radius, rep(ylim[2], nrow(dnf)), pch = 4, col = "firebrick", cex = 1.4)
  text(dnf$radius, rep(ylim[2], nrow(dnf)), "spiral DNF", pos = 1, cex = 0.7, col = "firebrick")
}
legend("topleft",
  c("spiral (pre-1.2.0.9004)", "nearest", "nearestRandom"),
  col = c("firebrick", "steelblue", "darkgreen"), pch = c(19, 19, 17), lty = c(1, 1, 2), bty = "n"
)
dev.off()
cat("wrote fig1_scaling_all_methods.png\n")
