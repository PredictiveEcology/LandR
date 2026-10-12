## Inputs and a runner for the result-identical speed-ups of leadingByPixelGroup() (used by
## vegTypeGenerator() and vegTypeMapGenerator()), .initiateNewCohorts(), calcSeverityB() and
## rmMissingCohorts(). fastCohortOpsRun() uses only functions that already existed before the
## speed-ups; fixtures/fastCohortOps-expected.rds holds what it returned on LandR 1.2.0.9052
## (f3c5a3e1, before them), made with tests/testthat/fixtures/make-fastCohortOps-expected.R.

fastOpsSppEquiv <- function() {
  data.table::data.table(
    LandR = c("Acer_rub", "Betu_pap", "Lari_lar", "Pice_gla", "Popu_tre"),
    Type = c("Deciduous", "Deciduous", "Conifer", "Conifer", "Deciduous")
  )
}

## One row per cohort. Species levels include "Abie_bal" (not in sppEquiv) and an unused level.
## Groups 8 and 9 hold two or more deciduous species whose proportions add up to exactly the
## 0.75 (group 8) and 0.25 (group 9) mixedwood thresholds when summed with base::sum(); added
## in plain double arithmetic (as a grouped sum does) they land one bit below / above them.
fastOpsCohorts <- function(bigPixelGroup = FALSE) {
  g <- function(pg, ...) {
    b <- c(...)
    data.table::data.table(pixelGroup = pg, speciesCode = names(b), B = unname(as.integer(b)))
  }
  cd <- data.table::rbindlist(list(
    g(1L, Pice_gla = 3000),
    g(2L, Pice_gla = 2000, Pice_gla = 1500, Popu_tre = 300),
    g(3L, Popu_tre = 8000, Pice_gla = 1000),
    g(4L, Pice_gla = 1000, Popu_tre = 1000),
    g(5L, Lari_lar = 500, Pice_gla = 500),
    g(6L, Acer_rub = 450, Betu_pap = 450, Pice_gla = 100),
    g(7L, Abie_bal = 400, Popu_tre = 600),
    g(8L, Pice_gla = 5406, Acer_rub = 15972, Betu_pap = 138, Popu_tre = 108),
    g(9L, Pice_gla = 20994, Acer_rub = 5974, Betu_pap = 1001, Popu_tre = 23),
    g(10L, Pice_gla = 5500, Acer_rub = 1500, Betu_pap = 1500, Popu_tre = 1500),
    g(11L, Betu_pap = 1000),
    g(12L, Pice_gla = 0, Popu_tre = 0),
    g(13L, Pice_gla = 8466, Acer_rub = 894, Betu_pap = 1719, Popu_tre = 209)
  ))
  if (bigPixelGroup) cd[pixelGroup == 13L, pixelGroup := 1000000L]
  cd[, age := seq_len(.N) * 10L]
  cd[, speciesCode := factor(speciesCode, levels = c("Abie_bal", "Acer_rub", "Betu_pap", "Lari_lar",
                                                      "Pice_gla", "Pinu_ban", "Popu_tre"))]
  cd[rev(seq_len(.N))]
}

fastOpsPixelGroupMap <- function() {
  r <- terra::rast(terra::ext(0, 4, 0, 4), resolution = 1)
  terra::values(r) <- c(13, 1:12, NA, 3, 8)
  r
}

## the seed after each call is part of the result: the calls must leave the RNG as they found it.
## The RNG kinds are pinned because R's defaults differ between versions (R-devel 2026-10), and
## sample() depends on them. Only the state is kept: R-devel also codes .Random.seed[1] differently
fastOpsSetSeed <- function(seed, envir = parent.frame()) {
  old <- RNGkind()
  withr::defer(suppressWarnings(RNGkind(old[1], old[2], old[3])), envir = envir)
  suppressWarnings(set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
                            sample.kind = "Rejection"))
}

## Where long double is plain double (e.g. macOS arm64), base::sum() lands on the other side of the
## mixedwood threshold for pixelGroups 8 and 13, so the original code itself decides them
## differently there and the fixtures, made on x86_64, do not apply
skipIfNoLongDouble <- function() {
  testthat::skip_if(.Machine$sizeof.longdouble <= 8L,
                    "long double is double here: exact-threshold fixtures were made with extended precision")
}

fastOpsLeading <- function() {
  sppEquiv <- fastOpsSppEquiv()
  cd <- fastOpsCohorts()
  seeded <- function(expr) {
    fastOpsSetSeed(1)
    out <- suppressMessages(expr)
    list(out = out, seed = .Random.seed[-1])
  }
  vtg <- function(x, ..., mixedType = 2) {
    seeded(vegTypeGenerator(x, mixedType = mixedType, sppEquiv = sppEquiv, sppEquivCol = "LandR",
                            doAssertion = FALSE, ...))
  }
  vtm <- function(x, mixedType = 2, doAssertion = FALSE, ...) {
    m <- seeded(vegTypeMapGenerator(x, fastOpsPixelGroupMap(), 0.75, mixedType = mixedType,
                                    sppEquiv = sppEquiv, sppEquivCol = "LandR",
                                    doAssertion = doAssertion, ...))
    list(values = terra::values(m$out, mat = FALSE), levels = terra::levels(m$out),
         coltab = terra::coltab(m$out), seed = m$seed)
  }
  naB <- data.table::copy(cd)[, B := {B[2] <- NA_integer_; B}]
  chr <- data.table::copy(cd)[, speciesCode := as.character(speciesCode)]
  cover <- data.table::copy(cd)[, `:=`(cover = B / 100, B = NULL)]
  noMixed <- cd[pixelGroup %in% c(1:3, 5L, 11L)]
  list(
    mixedType2 = vtg(cd, 0.75),
    mixedType2DefaultProp = vtg(cd),
    mixedType2Prop60 = vtg(cd, 0.6),
    mixedType2BigPixelGroup = vtg(fastOpsCohorts(bigPixelGroup = TRUE), 0.75),
    mixedType2NoMixed = vtg(noMixed, 0.75),
    mixedType1 = vtg(cd, 0.75, mixedType = 1),
    mixedType0 = vtg(cd, 0.75, mixedType = 0),
    characterSpeciesCode = vtg(chr, 0.75),
    naB = vtg(naB, 0.75),
    cover = vtg(cover, 0.75),
    map2 = vtm(cd),
    map2Assertion = vtm(cd, doAssertion = TRUE),
    map1 = vtm(cd, mixedType = 1),
    map1Assertion = vtm(cd, mixedType = 1, doAssertion = TRUE),
    map0Assertion = vtm(cd, mixedType = 0, doAssertion = TRUE)
  )
}

## ---- regeneration after fire: updateCohortData(), calcSeverityB(), rmMissingCohorts() ----------

fastOpsEcoregions <- function() factor(c("eco1", "eco2"))

## Pixel groups 1-4 on the map; pixels 5, 6, 12 and 13 burned (0); pixel 25 is NA.
## Pixel group 1 has its first row younger than successionTimestep (10) and older rows after it.
fastOpsRegenInputs <- function() {
  ecos <- fastOpsEcoregions()
  spp <- factor(c("Pice_gla", "Pinu_ban", "Popu_tre"))
  cd <- data.table::data.table(
    pixelGroup = c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L),
    ecoregionGroup = ecos[c(1, 1, 1, 2, 2, 1, 2, 2)],
    speciesCode = spp[c(1, 1, 3, 1, 2, 3, 1, 2)],
    age = c(5L, 80L, 40L, 30L, 60L, 15L, 12L, 9L),
    B = c(100L, 3000L, 1500L, 2000L, 900L, 800L, 600L, 50L),
    mortality = 0L, aNPPAct = 0L
  )
  cd[, sumB := sum(B), by = pixelGroup]
  pgm <- terra::rast(terra::ext(0, 5, 0, 5), resolution = 1)
  v <- rep(c(1, 2, 3, 4, 2), length.out = 25)
  v[c(5, 6, 12, 13)] <- 0
  v[25] <- NA
  terra::values(pgm) <- v
  speciesEcoregion <- data.table::data.table(
    year = rep(c(0, 5), each = 6),
    speciesCode = rep(rep(spp, 2), 2),
    ecoregionGroup = rep(rep(ecos, each = 3), 2),
    maxANPP = rep(c(300L, 200L, 250L, 280L, 190L, 240L), 2),
    maxB = rep(c(9000L, 6000L, 7000L, 8500L, 5500L, 6500L), 2)
  )
  list(cohortData = cd, pixelGroupMap = pgm, speciesEcoregion = speciesEcoregion,
       ecos = ecos, spp = spp, burned = c(5L, 6L, 12L, 13L))
}

fastOpsRegen <- function() {
  one <- function(initialB, partial = FALSE) {
    i <- fastOpsRegenInputs()
    new <- data.table::data.table(
      pixelIndex = i$burned,
      speciesCode = i$spp[c(1, 3, 1, 2)],
      ecoregionGroup = i$ecos[c(1, 1, 2, 2)]
    )
    if (partial) {
      ## a pixel that already has a pixelGroup also gets a new cohort: wholesale re-grouping
      new <- rbind(new, data.table::data.table(pixelIndex = 1L, speciesCode = i$spp[3],
                                                ecoregionGroup = i$ecos[1]))
    }
    fastOpsSetSeed(1)
    out <- suppressMessages(updateCohortData(
      new, i$cohortData, i$pixelGroupMap, currentTime = 7, speciesEcoregion = i$speciesEcoregion,
      successionTimestep = 10, initialB = initialB, verbose = 0, doAssertion = TRUE
    ))
    list(cohortData = out$cohortData, pixelGroupMap = terra::values(out$pixelGroupMap, mat = FALSE),
         seed = .Random.seed[-1])
  }
  ## .initiateNewCohorts() itself, with new cohorts in existing pixelGroups (1 is the pixelGroup
  ## whose first row is young); its sumB column is dropped by updateCohortData(), so is not compared
  initiate <- function(initialB) {
    i <- fastOpsRegenInputs()
    new <- data.table::data.table(
      pixelGroup = c(1L, 2L, 3L, 5L),
      speciesCode = i$spp[c(2, 3, 2, 1)],
      ecoregionGroup = i$ecos[c(1, 2, 1, 1)],
      age = 1L
    )
    out <- .initiateNewCohorts(new, i$cohortData, i$pixelGroupMap, currentTime = 7,
                               speciesEcoregion = i$speciesEcoregion, successionTimestep = 10,
                               initialB = initialB)
    if ("sumB" %in% names(out)) out[, sumB := NULL]
    list(cohortData = out, noOldSumB = !"oldSumB" %in% names(i$cohortData))
  }
  list(
    updateInitialB10 = one(10), updateInitialBNA = one(NA), updateInitialBNULL = one(NULL),
    updatePartialInitialB10 = one(10, TRUE), updatePartialInitialBNA = one(NA, TRUE),
    initiateInitialB10 = initiate(10), initiateInitialBNA = initiate(NA),
    initiateInitialBNULL = initiate(NULL)
  )
}

fastOpsSeverity <- function() {
  i <- fastOpsRegenInputs()
  cd <- i$cohortData
  burned <- data.table::data.table(pixelIndex = c(5L, 6L, 12L, 13L, 14L), pixelGroup = c(2L, 2L, 4L, 4L, 1L))
  postfire <- cd[burned, on = "pixelGroup", allow.cartesian = TRUE][, B := as.integer(B / 2)]
  list(
    withB = calcSeverityB(cd, postfire),
    standReplacing = suppressMessages(calcSeverityB(cd, data.table::copy(postfire)[, B := NULL]))
  )
}

fastOpsRmMissing <- function() {
  i <- fastOpsRegenInputs()
  run <- function(cd, pgm) {
    out <- rmMissingCohorts(cd, pgm, doAssertion = FALSE)
    list(cohortData = out$cohortData, pixelGroupMap = terra::values(out$pixelGroupMap, mat = FALSE))
  }
  pgmNothing <- i$pixelGroupMap
  pgmNothing[pgmNothing[] == 0] <- 1
  pgmCDGone <- pgmNothing
  pgmCDGone[pgmCDGone[] == 3] <- 1                      ## pixelGroup 3 is on cohortData only
  cdPGMGone <- data.table::copy(i$cohortData)[pixelGroup != 4L]     ## pixelGroup 4 is on the map only
  list(
    nothingToRemove = run(data.table::copy(i$cohortData), pgmNothing),
    cohortsGoneFromMap = run(data.table::copy(i$cohortData), pgmCDGone),
    pixelGroupGoneFromCohorts = run(cdPGMGone, pgmNothing),
    both = run(data.table::copy(cdPGMGone), pgmCDGone)
  )
}

## A data.table read back from an .rds has a dead self-reference pointer, which identical() sees
fastOpsNoSelfref <- function(x) {
  if (data.table::is.data.table(x)) {
    attr(x, ".internal.selfref") <- NULL
  } else if (is.list(x)) {
    x[] <- lapply(x, fastOpsNoSelfref)
  }
  x
}

fastCohortOpsRun <- function() {
  fastOpsNoSelfref(list(
    leading = fastOpsLeading(), regen = fastOpsRegen(), severity = fastOpsSeverity(),
    rmMissing = fastOpsRmMissing()
  ))
}
