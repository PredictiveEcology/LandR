## Helpers shared between seed-locked tests and dev benchmark scripts.
##
## All fixtures are built deterministically from `fixtureSeed`; no network
## downloads (no googledrive, no getSpeciesTable). The species table here is a
## small hand-coded subset roughly modelled after the BSW species set used in
## test-wardDispersalFunction.R but with values frozen so test outputs do not
## drift if the upstream table changes.

.landisDispSpeciesTable <- function(nSpecies = 7L) {
  ## Frozen species tables covering a range of dispersal regimes.
  ## Popu_tre (and its 16-species analogues) is the long-distance disperser but
  ## `eff` is held low enough that per-cell ward probability decays fast enough
  ## that no single species saturates (>60%) the eligible receivers in the
  ## test fixtures even with successionTimestep > 1.
  if (nSpecies == 7L) {
    data.table::data.table(
      species          = c("Abie_bal", "Betu_pap", "Pice_gla", "Pice_mar",
                           "Pinu_ban", "Popu_tre", "Lari_lar"),
      speciesCode      = 1:7,
      seeddistance_eff = c(  25,  200,   30,   80,   30,   120,   50),
      seeddistance_max = c( 160,  500,  200,  200,  100,  1500,  200)
    )
  } else if (nSpecies == 16L) {
    ## 16-species table — mix of short, mid, and long dispersers; designed
    ## so per-species saturation stays moderate even at successionTimestep=10
    data.table::data.table(
      species = sprintf("Sp_%02d", 1:16),
      speciesCode = 1:16,
      seeddistance_eff = c( 20,  40,  60,  80,
                           100, 120, 150, 180,
                            25,  50,  75, 100,
                            30,  90, 110,  60),
      seeddistance_max = c(120, 200, 250, 300,
                           400, 600, 800, 1000,
                           150, 220, 320, 450,
                           180, 350, 700, 240)
    )
  } else {
    stop("Only nSpecies %in% c(7, 16) is supported.")
  }
}

.landisDispSizeSpec <- function(size = c("tiny", "small", "medium", "large",
                                         "xlarge", "xlarge_dense", "xxlarge",
                                         "xxxlarge")) {
  size <- match.arg(size)
  switch(size,
    tiny         = list(nx =  20, ny =  20, res = 100, pgs =     8, propRcv = 0.5,
                        maxRcvSpPerPG = 4, maxSrcSpPerPG = 4),
    small        = list(nx =  50, ny =  50, res = 100, pgs =    30, propRcv = 0.5,
                        maxRcvSpPerPG = 5, maxSrcSpPerPG = 5),
    medium       = list(nx = 120, ny = 120, res = 100, pgs =   100, propRcv = 0.5,
                        maxRcvSpPerPG = 5, maxSrcSpPerPG = 5),
    large        = list(nx = 250, ny = 250, res = 100, pgs =   300, propRcv = 0.5,
                        maxRcvSpPerPG = 5, maxSrcSpPerPG = 5),
    xlarge       = list(nx = 500, ny = 500, res = 100, pgs =   800, propRcv = 0.5,
                        maxRcvSpPerPG = 5, maxSrcSpPerPG = 5),
    ## Same raster as xlarge but most pixelGroups are receivers — exercises a
    ## receiver-heavy regime where rcv cells dominate src cells (closer to a
    ## real fire-aftermath workflow where most of the burned landscape is
    ## eligible to receive seed and only a thin ring of survivors is source).
    xlarge_dense = list(nx = 500, ny = 500, res = 100, pgs =  1200, propRcv = 0.8,
                        maxRcvSpPerPG = 6, maxSrcSpPerPG = 5),
    ## Larger raster, also rcv-heavy. Use sparingly — the R reference takes
    ## ~10s on this size on a workstation; the Cpp port handles it in under 3s.
    xxlarge      = list(nx = 800, ny = 800, res = 100, pgs =  2400, propRcv = 0.7,
                        maxRcvSpPerPG = 6, maxSrcSpPerPG = 5),
    ## A real-landscape-scale stress test: 9M cells with ~4% receiver-eligible
    ## (the "small post-fire patch on a big landscape" regime). srcPixelMatrix
    ## uses ~250 MB; the inner loop iterates over hundreds of thousands of
    ## active receivers. R reference ~10s warm; Cpp ~4s.
    xxxlarge     = list(nx = 3000, ny = 3000, res = 100, pgs = 10000, propRcv = 0.04,
                        maxRcvSpPerPG = 6, maxSrcSpPerPG = 5)
  )
}

## Deterministic variable-size partition of `total` cells into `nBlocks` blocks.
## Pure integer arithmetic, no RNG — output is identical across OSes (unlike
## SpaDES.tools::randomPolygons, which is what this fixture used to call and
## which was the source of cross-platform divergence in seed-locked tests).
.varBlockSizes <- function(total, nBlocks, salt = 0L) {
  if (nBlocks <= 1L) return(as.integer(total))
  i <- seq_len(nBlocks) - 1L
  ## weight pattern in 1..5; (* 7) %% 5 cycles through {0,2,4,1,3,...}
  w <- as.integer(((i + as.integer(salt)) * 7L) %% 5L) + 1L
  cumW <- cumsum(w); totW <- cumW[nBlocks]
  bounds <- as.integer(round(cumW / totW * total))
  sizes <- diff(c(0L, bounds))
  ## guarantee no zero-size block (steal one cell from the largest)
  while (any(sizes == 0L)) {
    sizes[which.max(sizes)] <- sizes[which.max(sizes)] - 1L
    sizes[which(sizes == 0L)[1]] <- 1L
  }
  sizes
}

#' Build a deterministic LANDISDisp fixture without any network access.
#'
#' Reproducibility: pure integer arithmetic — no RNG, no SpaDES.tools, no
#' platform-sensitive raster ops. Identical bytes on Linux/Windows/macOS for
#' the same `(size, fixtureSeed, successionTimestep, nSpecies)`. `fixtureSeed`
#' is used as a salt to perturb the layout (block sizes, pgID rotation, species
#' cycle), so different seeds still produce visibly different fixtures.
makeLANDISDispFixture <- function(size = "small", fixtureSeed = 11L,
                                  successionTimestep = 10L,
                                  nSpecies = 7L) {
  spec         <- .landisDispSizeSpec(size)
  speciesTable <- .landisDispSpeciesTable(nSpecies = nSpecies)
  spCodes      <- speciesTable$speciesCode
  numSp        <- length(spCodes)
  salt         <- as.integer(fixtureSeed)

  ## ---- deterministic block-tiled pixelGroupMap ----
  ## Cells are tiled into rectangular blocks of varied (but deterministic) size.
  ## Block index (row-major) is mapped to a pgID via modulo nPgs, so distant
  ## blocks may share a pgID — giving multi-patch pixelGroups like a real
  ## pixelGroupMap (the property `randomPolygons` used to provide).
  nRowBlocks <- max(2L, as.integer(ceiling(sqrt(spec$pgs))))
  nColBlocks <- max(2L, as.integer(ceiling(spec$pgs / nRowBlocks)))
  rowSizes   <- .varBlockSizes(spec$ny, nRowBlocks, salt)
  colSizes   <- .varBlockSizes(spec$nx, nColBlocks, salt + 1L)
  rowBlock   <- rep(seq_along(rowSizes), times = rowSizes)
  colBlock   <- rep(seq_along(colSizes), times = colSizes)
  ## terra::rast fills `vals` row-major (upper-left → lower-right)
  ri <- rep(rowBlock, each = spec$nx)
  ci <- rep(colBlock, times = spec$ny)
  blockIdx <- (ri - 1L) * nColBlocks + (ci - 1L)
  vals <- as.integer(((blockIdx + salt) %% spec$pgs) + 1L)
  pgm <- terra::rast(
    xmin = 0, xmax = spec$nx * spec$res,
    ymin = 0, ymax = spec$ny * spec$res,
    resolution = c(spec$res, spec$res),
    vals = vals
  )

  ## ---- deterministic per-pg species assignment ----
  rcvPGs <- seq_len(round(spec$pgs * spec$propRcv))
  srcPGs <- setdiff(seq_len(spec$pgs), rcvPGs)

  ## For each pg: pick n in 1..min(maxSp, numSp) via (pg + salt) %% maxSp,
  ## then take n consecutive species codes (mod numSp) starting at a
  ## pg-and-salt-dependent offset. n distinct because n <= numSp.
  pgEntries <- function(pgs, maxSp) {
    if (length(pgs) == 0L) {
      return(data.table::data.table(pixelGroup = integer(0),
                                    speciesCode = integer(0)))
    }
    pgs   <- as.integer(pgs)
    maxSp <- min(as.integer(maxSp), numSp)
    nVec  <- ((pgs - 1L + salt) %% maxSp) + 1L
    rows  <- data.table::rbindlist(lapply(seq_along(pgs), function(k) {
      pg  <- pgs[k]; n <- nVec[k]
      off <- as.integer((pg - 1L) * 3L + salt)
      sps <- ((off + 0:(n - 1L)) %% numSp) + 1L
      data.table::data.table(pixelGroup = pg, speciesCode = sort(sps))
    }))
    rows
  }
  dtRcv <- pgEntries(rcvPGs, spec$maxRcvSpPerPG)
  dtSrc <- pgEntries(srcPGs, spec$maxSrcSpPerPG)

  ## dtRcv often arrives joined to speciesTable in real workflows
  dtRcvFull <- speciesTable[dtRcv, on = "speciesCode"]

  list(
    pixelGroupMap      = pgm,
    dtSrc              = dtSrc,
    dtRcv              = dtRcv,
    dtRcvFull          = dtRcvFull,
    speciesTable       = speciesTable,
    successionTimestep = successionTimestep,
    fixtureSeed        = fixtureSeed,
    size               = size,
    spec               = spec
  )
}

#' Run LANDISDisp on a fixture under a fixed RNG seed and return a normalized
#' result data.table suitable for snapshot comparison.
runLANDISDispOnFixture <- function(fix, runSeed, verbose = 1,
                                   landisDispFn = LANDISDisp, ...) {
  set.seed(runSeed)
  out <- landisDispFn(
    dtSrc              = fix$dtSrc,
    dtRcv              = fix$dtRcvFull,
    pixelGroupMap      = fix$pixelGroupMap,
    speciesTable       = fix$speciesTable,
    successionTimestep = fix$successionTimestep,
    verbose            = verbose,
    ...
  )
  ## strip attributes that depend on verbose path (ReasonForStop attr is huge
  ## and includes timestamps); keep only the user-facing columns
  out <- data.table::as.data.table(out)
  data.table::setattr(out, "ReasonForStop", NULL)
  cols <- intersect(c("pixelIndex", "speciesCode", "DistOfSuccess", "species"),
                    colnames(out))
  out <- out[, ..cols]
  data.table::setorderv(out, intersect(c("pixelIndex", "speciesCode"), cols))
  out[]
}
