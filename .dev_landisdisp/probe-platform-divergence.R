## Cross-platform divergence probe for LANDISDisp.
##
## Run this script on each OS (Linux, Windows, macOS) and compare the
## printed digests. Fills in this 2x3 table:
##
##   OS \ Impl |  R (useCpp=FALSE)  |  Cpp (useCpp=TRUE)
##   ----------+--------------------+--------------------
##   Linux     |        ?           |         ?
##   Windows   |        ?           |         ?
##   macOS     |        ?           |         ?
##
## Diagnosis:
##  - All four cells equal      => no divergence; bug is elsewhere.
##  - R column equal across OSes, Cpp column differs => Cpp port is platform-
##    sensitive (likely libm pow, FMA, or signed-overflow UB).
##  - R column differs across OSes => R reference itself is platform-sensitive
##    (less likely; would point at upstream pkgs or set.seed semantics).
##  - Cpp == R within an OS but both differ across OSes => fixture itself is
##    platform-sensitive (suspect SpaDES.tools::randomPolygons or terra).
##
## Usage:
##   Rscript .dev_landisdisp/probe-platform-divergence.R
##
## On a fresh machine that doesn't have LandR's LANDISDisp branch installed,
## set INSTALL = TRUE below (or pass --install).

INSTALL <- isTRUE("--install" %in% commandArgs(trailingOnly = TRUE))

if (INSTALL) {
  if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
  pak::pak("PredictiveEcology/LandR@LANDISDisp", ask = FALSE)
}

suppressPackageStartupMessages({
  library(LandR)
  library(data.table)
  library(terra)
  library(SpaDES.tools)
  library(digest)
})

cat("=== Environment ===\n")
cat("OS         :", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")
cat("R version  :", R.version.string, "\n")
cat("LandR ver  :", as.character(packageVersion("LandR")), "\n")
cat("terra ver  :", as.character(packageVersion("terra")), "\n")
cat("SpaDES.tls :", as.character(packageVersion("SpaDES.tools")), "\n")
cat("Endian     :", .Platform$endian, "\n\n")

## ---- deterministic fixture (mirrors helper-LANDISDisp-fixtures.R "tiny") ----
makeFix <- function(fixtureSeed = 11L, successionTimestep = 10L) {
  set.seed(fixtureSeed)
  pgm <- terra::rast(xmin = 0, xmax = 2000, ymin = 0, ymax = 2000,
                     resolution = c(100, 100), vals = 1L)
  pgm <- SpaDES.tools::randomPolygons(pgm, numTypes = 8L)
  pgm[] <- as.integer(terra::values(pgm))

  speciesTable <- data.table(
    species          = c("Abie_bal","Betu_pap","Pice_gla","Pice_mar",
                         "Pinu_ban","Popu_tre","Lari_lar"),
    speciesCode      = 1:7,
    seeddistance_eff = c(  25, 200,  30,  80,  30, 120,  50),
    seeddistance_max = c( 160, 500, 200, 200, 100,1500, 200)
  )
  spCodes <- speciesTable$speciesCode

  rcvPGs <- 1:4   # half receivers, half sources for pgs=8
  srcPGs <- 5:8

  rcvList <- lapply(rcvPGs, function(pg) {
    n <- sample.int(4L, 1L)
    data.table(pixelGroup = pg, speciesCode = sort(sample(spCodes, size = n)))
  })
  srcList <- lapply(srcPGs, function(pg) {
    n <- sample.int(4L, 1L)
    data.table(pixelGroup = pg, speciesCode = sort(sample(spCodes, size = n)))
  })
  dtRcv <- rbindlist(rcvList)
  dtSrc <- rbindlist(srcList)
  dtRcvFull <- speciesTable[dtRcv, on = "speciesCode"]

  list(pgm = pgm, dtSrc = dtSrc, dtRcvFull = dtRcvFull,
       speciesTable = speciesTable, successionTimestep = successionTimestep)
}

runOnce <- function(fix, runSeed, useCpp) {
  set.seed(runSeed)
  out <- LANDISDisp(
    dtSrc              = fix$dtSrc,
    dtRcv              = fix$dtRcvFull,
    pixelGroupMap      = fix$pgm,
    speciesTable       = fix$speciesTable,
    successionTimestep = fix$successionTimestep,
    verbose            = 1,
    useCpp             = useCpp
  )
  out <- as.data.table(out)
  setattr(out, "ReasonForStop", NULL)
  cols <- intersect(c("pixelIndex","speciesCode","DistOfSuccess","species"),
                    colnames(out))
  out <- out[, ..cols]
  setorderv(out, intersect(c("pixelIndex","speciesCode"), cols))
  out[]
}

## ---- fixture digest (sanity-check that the *input* is identical) ----
fix <- makeFix(fixtureSeed = 11L, successionTimestep = 10L)
fixDigest <- digest::digest(list(
  pgmVals  = as.integer(terra::values(fix$pgm)),
  dtSrc    = fix$dtSrc,
  dtRcvF   = fix$dtRcvFull,
  spTable  = fix$speciesTable
), algo = "sha256")
cat("=== Fixture digest (should match across OSes!) ===\n")
cat("  ", fixDigest, "\n\n")

## ---- run all four cells ----
cat("=== LANDISDisp output digests ===\n")
specs <- list(
  list(runSeed =   42L, label = "tiny fix=11 run=42 ts=10"),
  list(runSeed = 1729L, label = "tiny fix=11 run=1729 ts=10")
)
for (s in specs) {
  outR   <- runOnce(fix, runSeed = s$runSeed, useCpp = FALSE)
  outCpp <- runOnce(fix, runSeed = s$runSeed, useCpp = TRUE)
  cat(sprintf("[%s]\n", s$label))
  cat(sprintf("  R   (useCpp=FALSE): nrow=%5d  digest=%s\n",
              nrow(outR),   digest::digest(outR,   algo = "sha256")))
  cat(sprintf("  Cpp (useCpp=TRUE ): nrow=%5d  digest=%s\n",
              nrow(outCpp), digest::digest(outCpp, algo = "sha256")))
  cat(sprintf("  Cpp == R within this OS: %s\n\n",
              isTRUE(all.equal(outCpp, outR))))
}
