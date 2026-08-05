## Candidate replacements for the lowest-class tie-break.
##
## `method = "nearest"` only makes a choice where >= 2 classes tie at the minimum distance
## (35-41% of unwanted pixels). Resolving those to the lowest class code biases the result
## toward low-numbered classes, i.e. toward non-forest cover. This prototypes deterministic
## alternatives, scored the same way as everything else in this bundle: total-variation
## distance between the assigned-class composition and the old algorithm's, read against the
## old algorithm's own seed-to-seed noise floor.
##
##   lowestCode  -- current "nearest": ties -> lowest class code
##   modalTie    -- ties -> whichever tied class has most cells in the pixel's window
##                  (deterministic; falls back to lowest code only if the counts also tie)
##   modalAll    -- ignore ties; take the class with most cells in the window, full stop
##                  (the deterministic argmax of the same weights nearestRandom samples)
##   hashWeighted-- nearestRandom's weighted draw, but with the uniform taken from a hash of
##                  the pixel index instead of the RNG: deterministic and seed-free, while
##                  keeping the abundance weighting exactly
##
## Run: LCC_BENCH_DIR=<dir with v2_input_*.tif> LANDR_SRC=. \
##        Rscript benchmarks/convertUnwantedLCC-nearestRandom/07_tiebreak_candidates.R

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
SEEDS <- c(123L, 987L, 5150L)
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
  it <- 1L; cur <- length(unw); rep0 <- 0L; out3 <- NULL
  while (length(unw) > 0) {
    o <- SpaDES.tools::spread2(rstLCC, start = unw, asRaster = FALSE,
      iterations = it, allowOverlap = TRUE, spreadProb = 1)
    o <- o[initialPixels != pixels]; it <- it + 1L
    o[, lcc := as.vector(rstLCC[])[pixels]][lcc %in% classesToReplace, lcc := NA]
    o <- na.omit(o)
    o5 <- a[o[, state := NULL], allow.cartesian = TRUE,
      on = c("pixelIndex" = "initialPixels"), nomatch = NA]
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

## distance (in cells) and window counts for every (unwanted pixel, candidate class)
buildCand <- function(r, unw, cls) {
  v <- values(r)[, 1]
  nR <- nrow(r); nC <- ncol(r)
  unit <- rast(nrows = nR, ncols = nC, xmin = 0, xmax = nC, ymin = 0, ymax = nR, crs = "EPSG:3978")
  D <- vapply(cls, function(cc) {
    values(distance(setValues(unit, ifelse(v == cc, 1L, NA_integer_))))[unw, 1]
  }, numeric(length(unw)))
  minD <- do.call(pmin, as.data.frame(D))
  k <- pmin(nC, pmax(1L, as.integer(ceiling(minD))))
  cnt <- LandR:::windowCountsByClassCpp(
    lccVals = as.integer(v), candClasses = as.integer(cls),
    nrow = nR, ncol = nC, cells0 = as.integer(unw - 1L),
    kx = k, ky = pmin(nR, k)
  )
  self <- as.integer(v[unw])
  for (j in seq_along(cls)) {
    hit <- !is.na(self) & self == cls[j]
    cnt[hit, j] <- cnt[hit, j] - 1L
  }
  list(D = D, minD = minD, cnt = cnt, cls = cls)
}

## deterministic, seed-free uniform in [0,1) from the pixel index (xorshift-style mix)
hashUnif <- function(i) {
  x <- bitwXor(i, bitwShiftL(i, 13L))
  x <- bitwXor(x, bitwShiftR(x, 17L))
  x <- bitwXor(x, bitwShiftL(x, 5L))
  (abs(x) %% 100000L) / 100000
}

pickers <- list(
  lowestCode = function(cd) {
    tied <- cd$D <= cd$minD + 1e-9
    cd$cls[apply(tied, 1, function(z) which(z)[1L])]
  },
  modalTie = function(cd) {
    tied <- cd$D <= cd$minD + 1e-9
    w <- cd$cnt; w[!tied] <- -1L
    cd$cls[max.col(w, ties.method = "first")]
  },
  modalAll = function(cd) cd$cls[max.col(cd$cnt, ties.method = "first")],
  hashWeighted = function(cd) {
    u <- hashUnif(seq_len(nrow(cd$cnt)))
    w <- cd$cnt; w[w < 0L] <- 0L
    cw <- t(apply(w, 1, cumsum))
    tot <- cw[, ncol(cw)]
    target <- u * tot
    idx <- max.col(cw >= target & !cbind(FALSE, cw[, -ncol(cw), drop = FALSE] >= target),
                   ties.method = "first")
    cd$cls[idx]
  }
)

comp <- function(x, cls) {
  tab <- table(factor(x, levels = cls))
  as.numeric(tab) / max(1, sum(tab))
}
tvd <- function(p, q) 0.5 * sum(abs(p - q))

rows <- list()
for (lab in c("small", "medium", "large", "bigblob")) {
  f <- file.path(BENCH, paste0("v2_input_", lab, ".tif"))
  if (!file.exists(f)) next
  L <- rast(f)
  v <- values(L)[, 1]
  unw <- which(v == UNW)
  cls <- sort(unique(v[!is.na(v) & v != UNW & !v %in% EXCL]))
  aDT <- availDT(L)
  cd <- buildCand(L, unw, cls)

  cOldEach <- lapply(SEEDS, function(s) {
    set.seed(s)
    comp(as.integer(cuSpiral(UNW, L, aDT)$ecoregionGroup), cls)
  })
  cOld <- Reduce(`+`, cOldEach) / length(cOldEach)
  floorTVD <- mean(combn(length(cOldEach), 2, function(ij) tvd(cOldEach[[ij[1]]], cOldEach[[ij[2]]])))

  cRnd <- {
    set.seed(123)
    comp(as.integer(suppressMessages(convertUnwantedLCC(
      UNW, L, copy(aDT), doAssertion = FALSE, method = "nearestRandom"
    ))$ecoregionGroup), cls)
  }
  ## the shipped deterministic method, so this table scores what actually landed rather than
  ## only the R prototype of it
  cWtd <- comp(as.integer(suppressMessages(convertUnwantedLCC(
    UNW, L, copy(aDT), doAssertion = FALSE, method = "nearestWeighted"
  ))$ecoregionGroup), cls)

  res <- vapply(names(pickers), function(nm) tvd(comp(pickers[[nm]](cd), cls), cOld), numeric(1))
  rows[[lab]] <- data.table(
    landscape = lab, unwanted = length(unw),
    tvd_floor_oldSelf = round(floorTVD, 4),
    tvd_lowestCode = round(res[["lowestCode"]], 4),
    tvd_modalTie = round(res[["modalTie"]], 4),
    tvd_modalAll = round(res[["modalAll"]], 4),
    tvd_nearestWeighted = round(tvd(cWtd, cOld), 4),
    tvd_nearestRandom = round(tvd(cRnd, cOld), 4)
  )
  print(rows[[lab]])
  flush.console()
  rm(L, cd); gc(verbose = FALSE)
}

tab <- rbindlist(rows)
cat("\n== TVD from the old algorithm's assigned-class mix (lower = closer; floor = old vs old) ==\n")
print(as.data.frame(tab), row.names = FALSE)
fwrite(tab, file.path(OUT, "tiebreak_candidates.csv"))
