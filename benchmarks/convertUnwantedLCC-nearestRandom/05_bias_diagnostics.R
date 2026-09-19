## Where the deterministic allocation's bias actually lives.
##
## PR #196 raised the concern that a deterministic pick creates a visible artifact -- the
## recollection being that it "always chose the north east (or whatever) replacement". In
## *this* implementation the bias is not directional: `method = "nearest"` takes the class
## whose nearest cell is closest, and only when two or more classes tie at that distance
## does the tie-break decide -- always in favour of the lowest class code. So the artifact
## is a bias toward low-numbered classes, with no spatial signature.
##
## This script measures both, to show which one is real:
##
##  1. TIE BIAS -- of the unwanted pixels where >= 2 classes tie at the minimum distance,
##     the share assigned the lowest-coded of the tied classes. "nearest" is 100% by
##     construction; the old algorithm and "nearestRandom" should sit near chance.
##  2. SPATIAL STRUCTURE -- of all rook-adjacent pairs of unwanted pixels, the share
##     assigned the same class. A directional/patch artifact would push this up.
##
## Run: LCC_BENCH_DIR=<dir with v2_input_*.tif> LANDR_SRC=. \
##        Rscript benchmarks/convertUnwantedLCC-nearestRandom/05_bias_diagnostics.R

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  pkgload::load_all(Sys.getenv("LANDR_SRC", "."), quiet = TRUE)
})
terraOptions(memfrac = 0.2)

BENCH <- Sys.getenv("LCC_BENCH_DIR", "")
OUT <- file.path(Sys.getenv("LANDR_SRC", "."), "benchmarks", "convertUnwantedLCC-nearestRandom")
UNW <- 240L
EXCL <- c(0L, 20L, 30L)
paddedFloatToChar <- reproducible::paddedFloatToChar
.resample <- function(x, ...) x[sample.int(length(x), ...)]

cuSpiral <- function(classesToReplace, rstLCC, availableERC_by_Sp) {
  a <- data.table::copy(availableERC_by_Sp)
  if (!"speciesCode" %in% names(a)) a[, speciesCode := "allSpecies"]
  if (!"pixelIndex" %in% names(a)) a[, pixelIndex := seq(ncell(rstLCC))]
  unw <- sort(unique(a[initialEcoregionCode %in% classesToReplace, pixelIndex]))
  ERG2 <- unique(a[!initialEcoregionCode %in% classesToReplace],
    by = c("speciesCode", "initialEcoregionCode")
  )
  it <- 1L
  cur <- length(unw)
  rep0 <- 0L
  out3 <- NULL
  while (length(unw) > 0) {
    o <- SpaDES.tools::spread2(rstLCC,
      start = unw, asRaster = FALSE,
      iterations = it, allowOverlap = TRUE, spreadProb = 1
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
  unique(out3)
}

availDT <- function(r) {
  v <- values(r)[, 1]
  keep <- which(!is.na(v) & !v %in% EXCL)
  data.table(pixelIndex = keep, initialEcoregionCode = as.integer(v[keep]))
}

## raster of the assigned class, NA everywhere except the unwanted pixels
assignedOnly <- function(r, out) {
  v <- rep(NA_integer_, ncell(r))
  o <- out[!is.na(pixelIndex) & !is.na(ecoregionGroup)]
  v[o$pixelIndex] <- as.integer(o$ecoregionGroup)
  setValues(rast(r), v)
}

## of all rook-adjacent pairs where BOTH cells are unwanted, the fraction assigned the same
## class. 1 = the blob resolved into solid single-class patches.
adjAgreement <- function(a) {
  m <- matrix(values(a)[, 1], nrow = nrow(a), ncol = ncol(a), byrow = TRUE)
  h <- cbind(m[, -ncol(m)], NA)[, seq_len(ncol(m) - 1), drop = FALSE]
  pairsH <- cbind(as.vector(m[, -ncol(m)]), as.vector(m[, -1]))
  pairsV <- cbind(as.vector(m[-nrow(m), ]), as.vector(m[-1, ]))
  p <- rbind(pairsH, pairsV)
  p <- p[!is.na(p[, 1]) & !is.na(p[, 2]), , drop = FALSE]
  if (nrow(p) == 0) return(NA_real_)
  100 * mean(p[, 1] == p[, 2])
}

## for each unwanted pixel: which candidate classes tie at the minimum distance, and what is
## the lowest-coded of them. This is the only place `method = "nearest"` exercises a choice.
tieInfo <- function(r, unw) {
  v <- values(r)[, 1]
  cls <- sort(unique(v[!is.na(v) & v != UNW & !v %in% EXCL]))
  D <- vapply(cls, function(cc) {
    m <- setValues(rast(r), ifelse(v == cc, 1L, NA_integer_))
    values(distance(m))[unw, 1]
  }, numeric(length(unw)))
  minD <- do.call(pmin, as.data.frame(D))
  tied <- D <= minD + 1e-9
  nTied <- rowSums(tied)
  lowestTied <- cls[apply(tied, 1, function(z) which(z)[1L])] # cls is sorted ascending
  data.table(pixelIndex = unw, nTied = nTied, lowestTied = lowestTied)
}

## of the pixels with a genuine tie, the share given the lowest-coded tied class
tieBias <- function(out, ti) {
  m <- merge(out[!is.na(ecoregionGroup)], ti[nTied >= 2L], by = "pixelIndex")
  if (nrow(m) == 0) return(NA_real_)
  100 * mean(as.integer(m$ecoregionGroup) == m$lowestTied)
}

labs <- c("medium", "bigblob")
cols <- c("40" = "#fdbf6f", "50" = "#b2df8a", "100" = "#cab2d6",
          "210" = "#1b9e77", "220" = "#d95f02", "230" = "#7570b3", "240" = "grey55")

rows <- list()
panels <- list()
for (lab in labs) {
  f <- file.path(BENCH, paste0("v2_input_", lab, ".tif"))
  if (!file.exists(f)) {
    message("skipping ", lab, ": not found")
    next
  }
  L <- rast(f)
  aDT <- availDT(L)
  unw <- which(values(L)[, 1] == UNW)
  ti <- tieInfo(L, unw)

  set.seed(123)
  oldOut <- cuSpiral(UNW, L, aDT)
  ## the removed 1.2.0.9004 rule, reconstructed locally: nearest class, ties -> lowest code.
  ## `tieInfo()` already resolves exactly that, so no package call is needed (nor possible).
  lowOut <- data.table(pixelIndex = ti$pixelIndex, ecoregionGroup = ti$lowestTied)
  detOut <- suppressMessages(
    convertUnwantedLCC(UNW, L, copy(aDT), doAssertion = FALSE, method = "nearestWeighted")
  )
  set.seed(123)
  rndOut <- suppressMessages(
    convertUnwantedLCC(UNW, L, copy(aDT), doAssertion = FALSE, method = "nearestRandom")
  )
  old <- assignedOnly(L, oldOut)
  low <- assignedOnly(L, lowOut)
  det <- assignedOnly(L, detOut)
  rnd <- assignedOnly(L, rndOut)

  rows[[lab]] <- data.table(
    landscape = lab,
    unwanted = length(unw),
    pct_tied = round(100 * mean(ti$nTied >= 2L), 1),
    tie_spiral = round(tieBias(oldOut, ti), 1),
    tie_lowestCode = round(tieBias(lowOut, ti), 1),
    tie_nearestWeighted = round(tieBias(detOut, ti), 1),
    tie_nearestRandom = round(tieBias(rndOut, ti), 1),
    adj_spiral = round(adjAgreement(old), 1),
    adj_lowestCode = round(adjAgreement(low), 1),
    adj_nearestWeighted = round(adjAgreement(det), 1),
    adj_nearestRandom = round(adjAgreement(rnd), 1)
  )
  panels[[lab]] <- list(L = L, old = old, low = low, det = det, rnd = rnd)
  cat(sprintf(
    ">>> %-8s tied=%.1f%% | lowest-tied-class share: spiral=%.1f%% lowestCode=%.1f%% wtd=%.1f%% rand=%.1f%% | adjacency: spiral=%.1f%% lowestCode=%.1f%% wtd=%.1f%% rand=%.1f%%\n",
    lab, rows[[lab]]$pct_tied, rows[[lab]]$tie_spiral, rows[[lab]]$tie_lowestCode,
    rows[[lab]]$tie_nearestWeighted, rows[[lab]]$tie_nearestRandom, rows[[lab]]$adj_spiral,
    rows[[lab]]$adj_lowestCode, rows[[lab]]$adj_nearestWeighted, rows[[lab]]$adj_nearestRandom
  ))
  flush.console()
}

tab <- rbindlist(rows)
print(as.data.frame(tab), row.names = FALSE)
fwrite(tab, file.path(OUT, "bias_diagnostics.csv"))

png(file.path(OUT, "fig2_bias_diagnostics.png"),
  width = 1500, height = 340 * length(panels), res = 115
)
par(mfrow = c(length(panels), 5), mar = c(1.5, 1.5, 3, 1))
for (lab in names(panels)) {
  p <- panels[[lab]]
  plot(p$L, main = sprintf("%s: input (grey = unwanted)", lab),
    col = cols, type = "classes", legend = FALSE, axes = FALSE
  )
  plot(p$old, main = sprintf("spiral (old, stochastic) adj %.1f%%", rows[[lab]]$adj_spiral),
    col = cols, type = "classes", legend = FALSE, axes = FALSE
  )
  plot(p$low, main = sprintf("lowest-code (REMOVED) adj %.1f%%", rows[[lab]]$adj_lowestCode),
    col = cols, type = "classes", legend = FALSE, axes = FALSE
  )
  plot(p$det, main = sprintf("nearestWeighted adj %.1f%%", rows[[lab]]$adj_nearestWeighted),
    col = cols, type = "classes", legend = FALSE, axes = FALSE
  )
  plot(p$rnd, main = sprintf("nearestRandom adj %.1f%%", rows[[lab]]$adj_nearestRandom),
    col = cols, type = "classes", legend = FALSE, axes = FALSE
  )
}
dev.off()
cat("wrote fig2_spatial_artifacts.png\n")
