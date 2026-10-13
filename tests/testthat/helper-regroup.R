## Inputs and a runner for the result-identical speed-ups of updateCohortData()'s wholesale
## re-grouping, calcSiteShade() and rmMissingCohorts(). regroupRun() uses only functions that
## already existed before the speed-ups; fixtures/regroup-expected.rds holds what it returned on
## LandR 1.2.0.9054 (04bfdcd8, before them), made with tests/testthat/fixtures/make-regroup-expected.R.

regroupSpp <- function() factor(c("Abie_bal", "Pice_gla", "Pinu_ban", "Popu_tre"))
regroupEcos <- function() factor(c("eco1", "eco2", "eco3"))

## A 12 x 12 map: 20 pixelGroups, 14 emptied (burned) pixels, 6 NA pixels. cohortData is not sorted
## by pixelGroup. pixelGroups 5 and 6 are the same community but differ in `mortality`, so they merge.
regroupInputs <- function(seed = 1) {
  fastOpsSetSeed(seed)
  ecos <- regroupEcos(); spp <- regroupSpp()
  nPG <- 20L
  one <- function(pg) {
    n <- sample(1:4, 1)
    data.table::data.table(pixelGroup = pg, ecoregionGroup = ecos[((pg - 1L) %% 3L) + 1L],
                           speciesCode = spp[sample(4, n)],
                           age = as.integer(sample(c(2, 4, 6, 15, 40, 90, 140), n)),
                           B = as.integer(sample(c(10, 150, 800, 2500, 7000, 100000, 150000), n)))
  }
  cd <- data.table::rbindlist(lapply(seq_len(nPG), one))
  cd <- cd[cd$pixelGroup != 6L, , drop = FALSE]
  six <- data.table::copy(cd[cd$pixelGroup == 5L, , drop = FALSE])
  six[, pixelGroup := 6L]
  cd <- data.table::rbindlist(list(cd, six))
  cd[, `:=`(mortality = pixelGroup * 2L, aNPPAct = pixelGroup * 3L)]
  cd[pixelGroup == 6L, mortality := 99L]
  cd[, sumB := sum(B), by = pixelGroup]
  cd <- cd[sample(nrow(cd)), , drop = FALSE]
  pgm <- terra::rast(terra::ext(0, 12, 0, 12), resolution = 1)
  v <- c(7L, sample(c(1:6, 8:nPG), 143, replace = TRUE))
  v[sample(2:144, 14)] <- 0L
  v[sample(2:144, 6)] <- NA
  terra::values(pgm) <- v
  se <- data.table::CJ(year = c(0, 5), speciesCode = spp, ecoregionGroup = ecos)
  se[, `:=`(maxANPP = as.integer(200 + 10 * as.integer(speciesCode)),
            maxB = as.integer(5000 + 1000 * as.integer(speciesCode)))]
  list(cohortData = cd, pixelGroupMap = pgm, speciesEcoregion = se, ecos = ecos, spp = spp)
}

## new cohorts in emptied pixels (several pixels share a community, so they merge), in occupied pixels
## (which splits the pixelGroup), and "survivors" that updateCohortData() leaves out
regroupNewCohorts <- function(i, withType = TRUE) {
  v <- terra::values(i$pixelGroupMap, mat = FALSE)
  empty <- which(v == 0); occ <- which(v > 0)
  px <- c(empty, empty[1:6], occ[c(1, 2, 3, 10, 11)])
  nEmpty <- length(empty)
  new <- data.table::data.table(
    pixelIndex = px,
    speciesCode = i$spp[c(rep(2:3, length.out = nEmpty), 1, 1, 4, 4, 1, 1, 2, 3, 3, 4, 2)],
    ## a new cohort in an occupied pixel is in the ecoregion of the pixel's pixelGroup
    ecoregionGroup = i$ecos[c(rep(1:2, length.out = nEmpty), 3, 1, 3, 1, 1, 1, 2, 2, 3, 1, 1)[seq_along(px)]]
  )
  ## a pixel with two new cohorts has one ecoregion
  new[, ecoregionGroup := ecoregionGroup[1], by = pixelIndex]
  occupied <- seq_along(px) > nEmpty + 6L
  new[occupied, ecoregionGroup := i$ecos[((v[pixelIndex] - 1L) %% 3L) + 1L]]
  if (withType) {
    new[, type := rep(c("serotiny", "resprouting"), length.out = .N)]
    new <- rbind(new, data.table::data.table(pixelIndex = occ[20:22], speciesCode = i$spp[2],
                                              ecoregionGroup = i$ecos[1], type = "survivor"))
  }
  new
}

regroupUpdate <- function(seed, withType = TRUE, withTable = FALSE, B = FALSE) {
  i <- regroupInputs(seed)
  new <- regroupNewCohorts(i, withType)
  if (B) new[, B := as.integer(seq_len(.N) * 100000L)]
  tfp <- if (withTable) {
    data.table::data.table(pixelIndex = which(terra::values(i$pixelGroupMap, mat = FALSE) == 0))
  } else {
    NULL
  }
  fastOpsSetSeed(seed)
  out <- suppressMessages(updateCohortData(
    new, i$cohortData, i$pixelGroupMap, currentTime = 7, speciesEcoregion = i$speciesEcoregion,
    treedFirePixelTableSinceLastDisp = tfp, successionTimestep = 10, initialB = 10, verbose = 0,
    doAssertion = TRUE
  ))
  list(cohortData = out$cohortData, pixelGroupMap = terra::values(out$pixelGroupMap, mat = FALSE),
       seed = .Random.seed[-1])
}

regroupSiteShade <- function() {
  i <- regroupInputs(3)
  minRelativeB <- data.table::data.table(ecoregionGroup = i$ecos, X1 = 0.15, X2 = c(0.25, 0.3, 0.35),
                                         X3 = 0.5, X4 = 0.7, X5 = 0.9)
  one <- function(cd) {
    cd <- data.table::copy(cd)
    shade <- suppressMessages(calcSiteShade(5, cd, i$speciesEcoregion, minRelativeB))
    list(shade = shade, cohortData = cd)
  }
  cd <- i$cohortData
  young1 <- cd$pixelGroup == 1L & cd$age <= 5L
  list(mixedAges = one(cd),
       ## the first row of pixelGroup 1 is a young cohort: the shade is then read from its row
       youngFirst = one(data.table::rbindlist(list(cd[young1, , drop = FALSE], cd[!young1, , drop = FALSE]))),
       allYoung = one(data.table::copy(cd)[, age := pmin(age, 5L)]),
       withPrevMortality = one(data.table::copy(cd)[, `:=`(prevMortality = 4, sumB = 1L)]))
}

regroupRmMissing <- function() {
  i <- regroupInputs(4)
  run <- function(cd, pgm) {
    out <- rmMissingCohorts(cd, pgm, doAssertion = FALSE)
    list(cohortData = out$cohortData, pixelGroupMap = terra::values(out$pixelGroupMap, mat = FALSE))
  }
  pgm0 <- i$pixelGroupMap
  pgm0[pgm0[] == 0] <- 1
  pgmGone <- pgm0
  pgmGone[pgmGone[] == 3] <- 1                         ## pixelGroup 3 is on cohortData only
  cdGone <- data.table::copy(i$cohortData)
  cdGone <- cdGone[cdGone$pixelGroup != 4L, , drop = FALSE]    ## pixelGroup 4 is on the map only
  list(nothing = run(data.table::copy(i$cohortData), pgm0),
       cohortsGone = run(data.table::copy(i$cohortData), pgmGone),
       pixelsGone = run(cdGone, pgm0), both = run(data.table::copy(cdGone), pgmGone))
}

regroupRun <- function() {
  fastOpsNoSelfref(list(
    update = list(plain = regroupUpdate(1), other = regroupUpdate(2), noType = regroupUpdate(3, withType = FALSE),
                  withTable = regroupUpdate(4, withTable = TRUE), withB = regroupUpdate(5, B = TRUE)),
    siteShade = regroupSiteShade(), rmMissing = regroupRmMissing()
  ))
}
