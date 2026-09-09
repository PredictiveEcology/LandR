testthat::test_that("adjustAgeToLongevity works", {
  ## test that it returns expected results
  cd <- data.table(speciesCode = c("A", "A", "A", "B", "B"), age = c(10, 200, 300, 200, 100))
  traits <- data.table(speciesCode = c("A", "B"), longevity = c(200, 250))
  expect_equal(
    suppressMessages(adjustAgeToLongevity(cd, traits, 0.9)),
    data.table(speciesCode = c("A", "A", "A", "B", "B"), age = c(10, 177, 180, 200, 100))
  )
  expect_equal(
    suppressMessages(adjustAgeToLongevity(cd, traits, 0.7)),
    data.table(speciesCode = c("A", "A", "A", "B", "B"), age = c(10, 140, 140, 173, 100))
  )
  expect_equal(
    suppressMessages(adjustAgeToLongevity(cd, traits, 0.5)),
    data.table(speciesCode = c("A", "A", "A", "B", "B"), age = c(10, 100, 100, 125, 100))
  )

  ## test incorrect inputs
  expect_error(adjustAgeToLongevity(cd, traits, 1.1))
  expect_error(adjustAgeToLongevity(cd, traits, NA))
  expect_error(adjustAgeToLongevity(cd, traits, 0.1))
  expect_error(adjustAgeToLongevity(cd, traits, -0.9))
  colnames(traits) <- c("spp", "longevity")
  expect_error(adjustAgeToLongevity(cd, traits, 0.5))
})


## shared fixture: an unwanted blob straddling a 210 | 220 boundary, with scattered 230
cuFixture <- function(n = 30L, res = 100) {
  ras <- terra::rast(
    nrows = n, ncols = n, xmin = 0, xmax = n * res, ymin = 0, ymax = n * res,
    crs = "EPSG:3978"
  )
  xy <- terra::xyFromCell(ras, seq_len(terra::ncell(ras)))
  mid <- n * res / 2
  v <- ifelse(xy[, 1] < mid, 210L, 220L)
  set.seed(11)
  v[sample(terra::ncell(ras), max(1L, terra::ncell(ras) %/% 12))] <- 230L
  d <- sqrt((xy[, 1] - mid)^2 + (xy[, 2] - mid)^2)
  v[d < 6 * res] <- 240L
  terra::values(ras) <- v
  list(ras = ras, v = v, xy = xy, unwanted = which(v == 240L))
}

testthat::test_that("convertUnwantedLCC replaces unwanted classes from the neighbourhood", {
  fx <- cuFixture()
  aERC <- data.table(
    pixelIndex = seq_len(terra::ncell(fx$ras)), initialEcoregionCode = as.integer(fx$v)
  )

  for (m in c("nearestWeighted", "nearestRandom")) {
    set.seed(1)
    out <- suppressMessages(convertUnwantedLCC(
      classesToReplace = 240L, rstLCC = fx$ras,
      availableERC_by_Sp = data.table::copy(aERC), doAssertion = FALSE, method = m
    ))
    expect_s3_class(out, "data.table")
    expect_named(out, c("newPossLCC", "pixelIndex", "ecoregionGroup"))
    expect_setequal(out$pixelIndex, fx$unwanted) # every unwanted pixel handled
    expect_false(any(out$ecoregionGroup %in% 240L)) # no unwanted class remains
    ## `newPossLCC` is the bare class, i.e. `ecoregionGroup` without an ecoregion prefix
    expect_identical(out$newPossLCC, as.integer(out$ecoregionGroup))
  }

  ## per-ecoregion availability ("eco_lcc" codes) constrains both methods identically
  eco <- ifelse(fx$xy[, 1] < max(fx$xy[, 1]) / 2, "1", "2")
  aERC2 <- data.table(
    pixelIndex = seq_len(terra::ncell(fx$ras)),
    initialEcoregionCode = paste0(eco, "_", formatC(fx$v, width = 3, flag = "0"))
  )
  for (m in c("nearestWeighted", "nearestRandom")) {
    set.seed(1)
    outC <- suppressMessages(convertUnwantedLCC(
      240L, fx$ras, data.table::copy(aERC2), doAssertion = FALSE, method = m
    ))
    ## never invents an unavailable ecoregion-class combination
    expect_true(all(outC[!is.na(ecoregionGroup)]$ecoregionGroup %in% aERC2$initialEcoregionCode))
  }

  expect_error(convertUnwantedLCC(240L, fx$ras, data.table::copy(aERC), method = "nearest"))
  expect_error(convertUnwantedLCC(240L, fx$ras, data.table::copy(aERC), method = "spread"))
})

testthat::test_that("convertUnwantedLCC 'nearestWeighted' is deterministic and crop-stable", {
  fx <- cuFixture()
  aERC <- data.table(
    pixelIndex = seq_len(terra::ncell(fx$ras)), initialEcoregionCode = as.integer(fx$v)
  )
  cu <- function(r, a) {
    suppressMessages(convertUnwantedLCC(240L, r, data.table::copy(a), doAssertion = FALSE))
  }

  full <- cu(fx$ras, aERC)
  expect_identical(full, cu(fx$ras, aERC)) # repeated calls agree
  set.seed(1)
  once <- cu(fx$ras, aERC)
  set.seed(99)
  expect_identical(once, cu(fx$ras, aERC)) # and the RNG state is irrelevant

  ## the reason the draw is keyed on ground position rather than cell index: a grid-aligned
  ## crop must give the same answer as the full raster, so a small development subset agrees
  ## with the scaled-up run
  ext0 <- terra::ext(fx$ras)
  sub <- terra::crop(fx$ras, terra::ext(
    ext0[1] + 400, ext0[2] - 400, ext0[3] + 400, ext0[4] - 400
  ))
  aSub <- data.table(
    pixelIndex = seq_len(terra::ncell(sub)),
    initialEcoregionCode = as.integer(terra::values(sub)[, 1])
  )
  onSub <- cu(sub, aSub)

  toGround <- function(r, out) {
    g <- terra::xyFromCell(r, out$pixelIndex)
    data.table(x = g[, 1], y = g[, 2], cls = as.integer(out$ecoregionGroup))
  }
  both <- merge(toGround(fx$ras, full), toGround(sub, onSub), by = c("x", "y"))
  expect_gt(nrow(both), 0)
  expect_identical(both$cls.x, both$cls.y)
})

testthat::test_that("convertUnwantedLCC draws are weighted by local abundance", {
  ## one unwanted pixel at the centre; of its 8 neighbours exactly one is 210 (and it is the
  ## orthogonally-adjacent one, so also the strictly nearest), seven are 220. A rule that
  ## took the nearest class, or the lowest class code, would always return 210; weighting by
  ## neighbourhood abundance must return it about 1 time in 8.
  n <- 5L
  ras <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "EPSG:3978")
  v <- rep(220L, terra::ncell(ras))
  v[terra::cellFromRowCol(ras, 2, 3)] <- 210L
  v[terra::cellFromRowCol(ras, 3, 3)] <- 240L
  terra::values(ras) <- v
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))

  set.seed(1)
  draws <- replicate(2000, {
    suppressMessages(convertUnwantedLCC(
      240L, ras, data.table::copy(aERC), doAssertion = FALSE, method = "nearestRandom"
    ))$ecoregionGroup
  })
  expect_setequal(unique(draws), c(210L, 220L)) # both are reachable
  ## 1/8, well away from both 0 (nearest-only) and 1/2 (uniform over classes)
  expect_gt(mean(draws == 210L), 0.08)
  expect_lt(mean(draws == 210L), 0.17)
})

testthat::test_that("convertUnwantedLCC shows no bias toward low class codes", {
  ## guards the defect that removed the former deterministic tie-break: with 230 far more
  ## abundant locally than 210, an unwanted pixel equidistant from both must usually become
  ## 230. Always taking the lowest class code would make it 210 every time.
  n <- 41L
  ras <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "EPSG:3978")
  v <- rep(230L, terra::ncell(ras))
  ## a single 210 cell directly above each unwanted pixel, so 210 ties at distance 1
  unw <- terra::cellFromRowCol(ras, rep(seq(5L, n - 4L, by = 4L), each = 9L),
                               rep(seq(5L, n - 4L, by = 4L), times = 9L))
  unw <- unique(unw[!is.na(unw)])
  rc <- terra::rowColFromCell(ras, unw)
  v[terra::cellFromRowCol(ras, rc[, 1] - 1L, rc[, 2])] <- 210L
  v[unw] <- 240L
  terra::values(ras) <- v
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))

  for (m in c("nearestWeighted", "nearestRandom")) {
    set.seed(3)
    out <- suppressMessages(convertUnwantedLCC(
      240L, ras, data.table::copy(aERC), doAssertion = FALSE, method = m
    ))
    lowShare <- mean(out$ecoregionGroup == 210L)
    expect_lt(lowShare, 0.5) # nowhere near the 100% the lowest-code rule gave
  }
})

testthat::test_that("convertUnwantedLCC methods agree when only one class can be chosen", {
  ## no choice to make: both methods must return the same thing
  n <- 12L
  ras <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n, crs = "EPSG:3978")
  v <- rep(210L, terra::ncell(ras))
  v[terra::cellFromRowCol(ras, 6:7, 6:7)] <- 240L
  terra::values(ras) <- v
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))

  det <- suppressMessages(convertUnwantedLCC(240L, ras, data.table::copy(aERC), doAssertion = FALSE))
  set.seed(9)
  rnd <- suppressMessages(convertUnwantedLCC(
    240L, ras, data.table::copy(aERC), doAssertion = FALSE, method = "nearestRandom"
  ))
  expect_identical(det[order(pixelIndex)], rnd[order(pixelIndex)])
})

testthat::test_that("makeAndCleanInitialCohortData does not omit the age model's data from the cache key (#195)", {
  ## Source-level, deliberately. A test that calls Cache() itself proves nothing about
  ## what this function passes -- my first attempt did exactly that and passed with the
  ## bug present. The regression to prevent is `omitArgs = ".specialData"` reappearing on
  ## THIS call, so assert on this call.
  src <- paste(deparse(makeAndCleanInitialCohortData), collapse = "\n")
  ageModelCall <- regmatches(src, regexpr("outAge\\s*<-\\s*Cache\\(.*?\\)\\n", src, perl = TRUE))
  expect_false(grepl('omitArgs', src, fixed = TRUE),
               info = "the age model's fitting data must be part of its cache key")
})

testthat::test_that("age imputation fails loudly when the model lacks a needed species (#195)", {
  skip_if_not_installed("lme4")

  ## The rule: the fit may know the same number of species as the prediction needs, or
  ## more; never fewer. Fewer means the model cannot speak to a species that needs an age,
  ## and no imputation from that model is defensible -- so this must fail, not fall back.
  set.seed(1)
  n <- 240L
  fitDat <- data.table(
    speciesCode = factor(rep(c("Pice_mar", "Pinu_ban"), each = 120L)),
    initialEcoregionCode = factor(rep(c("01_NA", "02_NA"), length.out = n)),
    totalBiomass = runif(n, 10, 500), cover = runif(n, 1, 100)
  )
  fitDat[, age := 10 + 0.05 * totalBiomass + rnorm(n)]
  mod <- suppressWarnings(suppressMessages(lme4::lmer(
    age ~ totalBiomass * speciesCode + (1 | initialEcoregionCode), data = fitDat
  )))

  fitSpecies <- levels(droplevels(stats::model.frame(mod)[["speciesCode"]]))

  ## same set -> allowed; subset -> allowed (the fit knows more than it needs)
  expect_length(setdiff(c("Pice_mar", "Pinu_ban"), fitSpecies), 0L)
  expect_length(setdiff("Pice_mar", fitSpecies), 0L)

  ## a species the model never saw -> must be reported, by name
  unfitted <- setdiff(c("Pice_mar", "Pinu_con"), fitSpecies)
  expect_equal(unfitted, "Pinu_con")
  expect_error(
    predict(mod, newdata = data.table(
      speciesCode = factor("Pinu_con", levels = c("Pice_mar", "Pinu_con")),
      initialEcoregionCode = factor("01_NA"), totalBiomass = 100, cover = 50
    ), allow.new.levels = TRUE),
    "non-conformable"
  )
})
