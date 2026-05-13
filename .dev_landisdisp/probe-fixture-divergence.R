## Narrow probe: which part of the fixture is platform-divergent?
##
## We already know R == Cpp on each OS, but both differ Linux<->Windows for
## the same seed+fixture. The fixture digest itself differs across OSes — so
## one of: pgm values, dtSrc, dtRcvFull, OR the RNG state advances differently
## somewhere upstream of these.
##
## This script hashes each component separately and runs two variants:
##   (A) "stock"     — uses SpaDES.tools::randomPolygons (the real fixture)
##   (B) "handMade"  — uses a deterministic hand-built pixelGroupMap
##
## If (B) matches across OSes but (A) doesn't, randomPolygons/terra is the
## divergence point. If (B) also differs, it's the LANDISDisp call itself
## (terra cell ordering, spiral, etc.).
##
## Run on Linux + Windows; paste both outputs.

suppressPackageStartupMessages({
  library(LandR)
  library(data.table)
  library(terra)
  library(SpaDES.tools)
  library(digest)
})

cat("OS:", Sys.info()[["sysname"]], " R:", R.version.string,
    " terra:", as.character(packageVersion("terra")),
    " SpaDES.tools:", as.character(packageVersion("SpaDES.tools")), "\n\n")

speciesTable <- data.table(
  species          = c("Abie_bal","Betu_pap","Pice_gla","Pice_mar",
                       "Pinu_ban","Popu_tre","Lari_lar"),
  speciesCode      = 1:7,
  seeddistance_eff = c(  25, 200,  30,  80,  30, 120,  50),
  seeddistance_max = c( 160, 500, 200, 200, 100,1500, 200)
)
spCodes <- speciesTable$speciesCode
sha <- function(x) substr(digest::digest(x, algo = "sha256"), 1, 16)

## ---- (A) Stock fixture: randomPolygons + sample() ----
cat("--- (A) randomPolygons-based fixture ---\n")
set.seed(11L)
pgmEmpty <- terra::rast(xmin = 0, xmax = 2000, ymin = 0, ymax = 2000,
                        resolution = c(100, 100), vals = 1L)

## State BEFORE randomPolygons
rngBefore <- .Random.seed
pgmA <- SpaDES.tools::randomPolygons(pgmEmpty, numTypes = 8L)
rngAfter <- .Random.seed
pgmA[] <- as.integer(terra::values(pgmA))

cat("  pgmA values  digest :", sha(as.integer(terra::values(pgmA))), "\n")
cat("  pgmA dim/range      :", paste(dim(pgmA), collapse = "x"),
    " range=", paste(range(terra::values(pgmA)), collapse = "-"), "\n")
cat("  pgmA tabulation     :",
    paste(as.integer(table(terra::values(pgmA))), collapse = ","), "\n")
cat("  RNG state delta     : before=", sha(rngBefore),
    " after=", sha(rngAfter), "\n")

## Continue with the stock sampling
rcvPGs <- 1:4; srcPGs <- 5:8
rcvList <- lapply(rcvPGs, function(pg) {
  n <- sample.int(4L, 1L)
  data.table(pixelGroup = pg, speciesCode = sort(sample(spCodes, size = n)))
})
srcList <- lapply(srcPGs, function(pg) {
  n <- sample.int(4L, 1L)
  data.table(pixelGroup = pg, speciesCode = sort(sample(spCodes, size = n)))
})
dtRcvA <- rbindlist(rcvList); dtSrcA <- rbindlist(srcList)
dtRcvFullA <- speciesTable[dtRcvA, on = "speciesCode"]

cat("  dtSrcA digest       :", sha(dtSrcA), "  nrow=", nrow(dtSrcA), "\n")
cat("  dtRcvFullA digest   :", sha(dtRcvFullA), "  nrow=", nrow(dtRcvFullA), "\n")

set.seed(42L)
outA <- LANDISDisp(dtSrc = dtSrcA, dtRcv = dtRcvFullA,
                   pixelGroupMap = pgmA, speciesTable = speciesTable,
                   successionTimestep = 10L, verbose = 1, useCpp = FALSE)
outA <- as.data.table(outA); setattr(outA, "ReasonForStop", NULL)
outA <- outA[, intersect(c("pixelIndex","speciesCode","DistOfSuccess","species"),
                         colnames(outA)), with = FALSE]
setorderv(outA, intersect(c("pixelIndex","speciesCode"), colnames(outA)))
cat("  LANDISDisp(A) out   :", sha(outA), "  nrow=", nrow(outA), "\n\n")


## ---- (B) Hand-built pgm: skip randomPolygons entirely ----
cat("--- (B) hand-built fixture (no randomPolygons) ---\n")
## Deterministic checkerboard-ish pgm: 8 groups, fixed assignment by
## (row %/% 5, col %/% 5). 20x20 = 400 cells; 16 blocks; mod 8 -> 8 groups.
nx <- 20; ny <- 20
rowIdx <- rep(seq_len(ny), each = nx)
colIdx <- rep(seq_len(nx), times = ny)
vals <- as.integer(((rowIdx %/% 5L) * 4L + (colIdx %/% 5L)) %% 8L + 1L)
pgmB <- terra::rast(xmin = 0, xmax = 2000, ymin = 0, ymax = 2000,
                    resolution = c(100, 100), vals = vals)
pgmB[] <- as.integer(terra::values(pgmB))

cat("  pgmB values digest  :", sha(as.integer(terra::values(pgmB))), "\n")
cat("  pgmB tabulation     :",
    paste(as.integer(table(terra::values(pgmB))), collapse = ","), "\n")

## Hand-coded src/rcv (no RNG at all in fixture)
dtRcvB <- data.table(
  pixelGroup  = c(1,1,2,2,3,3,4,4),
  speciesCode = c(1,2,3,4,5,6,7,1)
)
dtSrcB <- data.table(
  pixelGroup  = c(5,5,6,6,7,7,8,8),
  speciesCode = c(2,3,4,5,6,7,1,2)
)
dtRcvFullB <- speciesTable[dtRcvB, on = "speciesCode"]

cat("  dtSrcB digest       :", sha(dtSrcB), "\n")
cat("  dtRcvFullB digest   :", sha(dtRcvFullB), "\n")

set.seed(42L)
outB <- LANDISDisp(dtSrc = dtSrcB, dtRcv = dtRcvFullB,
                   pixelGroupMap = pgmB, speciesTable = speciesTable,
                   successionTimestep = 10L, verbose = 1, useCpp = FALSE)
outB <- as.data.table(outB); setattr(outB, "ReasonForStop", NULL)
outB <- outB[, intersect(c("pixelIndex","speciesCode","DistOfSuccess","species"),
                         colnames(outB)), with = FALSE]
setorderv(outB, intersect(c("pixelIndex","speciesCode"), colnames(outB)))
cat("  LANDISDisp(B) out   :", sha(outB), "  nrow=", nrow(outB), "\n")
