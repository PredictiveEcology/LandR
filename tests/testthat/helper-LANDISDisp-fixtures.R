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

#' Build a deterministic LANDISDisp fixture without any network access.
#'
#' Reproducibility note: this uses set.seed() *internally*; the surrounding
#' RNG state is restored on exit so a caller can later set their own seed for
#' the dispersal call itself. That way the fixture ordering is independent of
#' the seed used to drive LANDISDisp.
makeLANDISDispFixture <- function(size = "small", fixtureSeed = 11L,
                                  successionTimestep = 10L,
                                  nSpecies = 7L) {
  spec <- .landisDispSizeSpec(size)

  ## isolate fixture RNG so callers can use their own seeds afterwards
  oldSeed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(oldSeed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(fixtureSeed)

  pgm <- terra::rast(
    xmin = 0, xmax = spec$nx * spec$res,
    ymin = 0, ymax = spec$ny * spec$res,
    resolution = c(spec$res, spec$res),
    vals = 1L
  )
  pgm <- SpaDES.tools::randomPolygons(pgm, numTypes = spec$pgs)
  ## ensure integer storage to mimic real pixelGroupMaps
  pgm[] <- as.integer(terra::values(pgm))

  speciesTable <- .landisDispSpeciesTable(nSpecies = nSpecies)
  spCodes <- speciesTable$speciesCode

  rcvPGs <- seq_len(round(spec$pgs * spec$propRcv))
  srcPGs <- setdiff(seq_len(spec$pgs), rcvPGs)

  rcvList <- lapply(rcvPGs, function(pg) {
    n <- sample.int(spec$maxRcvSpPerPG, 1L)
    data.table::data.table(
      pixelGroup  = pg,
      speciesCode = sort(sample(spCodes, size = n))
    )
  })
  srcList <- lapply(srcPGs, function(pg) {
    n <- sample.int(spec$maxSrcSpPerPG, 1L)
    data.table::data.table(
      pixelGroup  = pg,
      speciesCode = sort(sample(spCodes, size = n))
    )
  })

  dtRcv <- data.table::rbindlist(rcvList)
  dtSrc <- data.table::rbindlist(srcList)

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
