## Self-contained scaling benchmark for convertUnwantedLCC.
## Compares the previous iterative-spread2 implementation ("spiral") against the new
## terra::distance() nearest-allocation, on a synthetic landscape with a single unwanted
## blob of increasing radius. Demonstrates the O(radius^2) -> ~O(1) change and that the
## two agree within the old algorithm's own (random tie-break) run-to-run variance.
##
## Run against this branch:  LANDR_SRC=. Rscript benchmarks/convertUnwantedLCC-nearest-alloc/01_scaling_benchmark.R
suppressPackageStartupMessages({ library(terra); library(data.table); library(LandR) })
paddedFloatToChar <- reproducible::paddedFloatToChar
.resample <- function(x, ...) x[sample.int(length(x), ...)]   # LandR-internal helper, inlined

## new implementation (this branch): eval from source if provided, else installed LandR
LANDR_SRC <- Sys.getenv("LANDR_SRC", "")
cuNew <- LandR::convertUnwantedLCC
if (nzchar(LANDR_SRC) && file.exists(file.path(LANDR_SRC, "R/cohorts.R"))) {
  for (e in parse(file.path(LANDR_SRC, "R/cohorts.R")))
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("convertUnwantedLCC"))) { cuNew <- eval(e); break }
}

## previous implementation (iterative spread2 "spiral"), for a same-process comparison
cuSpiral <- function(classesToReplace, rstLCC, availableERC_by_Sp) {
  a <- data.table::copy(availableERC_by_Sp)
  hasPreDash <- all(grepl("_", a$initialEcoregionCode))
  if (!"speciesCode" %in% names(a)) a[, speciesCode := "allSpecies"]
  if (!"pixelIndex" %in% names(a)) a[, pixelIndex := seq(ncell(rstLCC))]
  unw <- sort(unique(a[gsub(".*_", "", initialEcoregionCode) %in% as.character(classesToReplace), pixelIndex]))
  ERG2 <- unique(a[!gsub(".*_", "", initialEcoregionCode) %in% classesToReplace], by = c("speciesCode", "initialEcoregionCode"))
  it <- 1L; cur <- length(unw); rep0 <- 0L; out3 <- NULL
  while (length(unw) > 0) {
    o <- SpaDES.tools::spread2(rstLCC, start = unw, asRaster = FALSE, iterations = it, allowOverlap = TRUE, spreadProb = 1)
    o <- o[initialPixels != pixels]; it <- it + 1L
    o[, lcc := as.vector(rstLCC[])[pixels]][lcc %in% classesToReplace, lcc := NA]; o <- na.omit(o)
    o5 <- a[o[, state := NULL], allow.cartesian = TRUE, on = c("pixelIndex" = "initialPixels"), nomatch = NA]
    o5[, possERC := lcc]
    o6 <- na.omit(o5[ERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode"), nomatch = NA])
    rm <- o5[!ERG2, on = c("speciesCode", "possERC" = "initialEcoregionCode")]
    o6 <- o6[!possERC %in% unique(rm$possERC)]
    if (cur == length(unw)) rep0 <- rep0 + 1L else { cur <- length(unw); rep0 <- 0L }
    if (rep0 > 5) unw <- integer()
    if (nrow(o6) > 0) {
      keep <- o6[, list(k = .resample(.I, 1)), by = pixelIndex]
      o2 <- o6[keep$k][, list(pixelIndex, ecoregionGroup = as.integer(lcc))]
      unw <- unw[!unw %in% o2$pixelIndex]
      out3 <- if (is.null(out3)) o2 else rbindlist(list(o2, out3))
    }
  }
  unique(out3)
}

blob <- function(radius) {
  n <- as.integer(radius * 2.6)
  r <- rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  xy <- xyFromCell(r, 1:ncell(r)); v <- ifelse(xy[, 1] < n / 2, 210L, 220L)
  d <- sqrt((xy[, 1] - n / 2)^2 + (xy[, 2] - n / 2)^2); v[d < radius] <- 240L
  setValues(r, v)
}

radii <- as.integer(strsplit(Sys.getenv("BENCH_RADII", "10,20,40,80"), ",")[[1]])
res <- rbindlist(lapply(radii, function(rad) {
  r <- blob(rad); v <- values(r)[, 1]
  aDT <- data.table(pixelIndex = seq_len(ncell(r)), initialEcoregionCode = as.integer(v))
  tOld <- system.time(oOld <- suppressMessages(cuSpiral(240L, r, aDT)))[["elapsed"]]
  tNew <- system.time(oNew <- suppressMessages(cuNew(240L, r, data.table::copy(aDT), doAssertion = FALSE)))[["elapsed"]]
  uw <- which(v == 240L)
  ag <- 100 * mean(oOld[match(uw, pixelIndex), ecoregionGroup] == oNew[match(uw, pixelIndex), ecoregionGroup], na.rm = TRUE)
  data.table(radius = rad, unwanted = length(uw), spiral_s = round(tOld, 2), nearest_s = round(tNew, 2),
             speedup = round(tOld / max(tNew, 1e-6)), agree_pct = round(ag, 1))
}))
print(res)

png(file.path(Sys.getenv("BENCH_OUT", "."), "fig1_scaling.png"), width = 1100, height = 460, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(res$radius, res$spiral_s, type = "b", log = "y", pch = 19, col = "firebrick",
     xlab = "unwanted blob radius (cells)", ylab = "time (s, log)", main = "run time vs blob radius",
     ylim = range(c(res$spiral_s, res$nearest_s)))
lines(res$radius, res$nearest_s, type = "b", pch = 19, col = "steelblue")
legend("topleft", c("spiral (old)", "nearest-alloc (new)"), col = c("firebrick", "steelblue"), pch = 19, bty = "n")
plot(res$radius, res$agree_pct, type = "b", pch = 19, ylim = c(0, 100),
     xlab = "unwanted blob radius (cells)", ylab = "new vs old agreement (%)", main = "agreement (imputation)")
dev.off()
cat("wrote fig1_scaling.png\n")
