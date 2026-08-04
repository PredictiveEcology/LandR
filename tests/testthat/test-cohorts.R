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

testthat::test_that("convertUnwantedLCC replaces unwanted classes with nearest available class", {
  n <- 30L
  ras <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  xy <- terra::xyFromCell(ras, seq_len(terra::ncell(ras)))
  v <- ifelse(xy[, 1] < n / 2, 210L, 220L) # class 210 (left) / 220 (right)
  d <- sqrt((xy[, 1] - n / 2)^2 + (xy[, 2] - n / 2)^2)
  v[d < 6] <- 240L # an unwanted blob straddling the 210|220 boundary
  terra::values(ras) <- v
  unwanted <- which(v == 240L)

  ## --- unconstrained: any non-unwanted class may replace 240 ---
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))
  out <- suppressMessages(convertUnwantedLCC(
    classesToReplace = 240L,
    rstLCC = ras,
    availableERC_by_Sp = data.table::copy(aERC),
    doAssertion = FALSE
  ))
  expect_s3_class(out, "data.table")
  expect_named(out, c("newPossLCC", "pixelIndex", "ecoregionGroup"))
  expect_setequal(out$pixelIndex, unwanted) # every unwanted pixel handled
  expect_false(any(out$ecoregionGroup %in% 240L)) # no unwanted class remains
  expect_setequal(unique(out$ecoregionGroup), c(210L, 220L))
  ## `newPossLCC` is the bare class, i.e. `ecoregionGroup` without an ecoregion prefix
  expect_identical(out$newPossLCC, as.integer(out$ecoregionGroup))
  ## deterministic (the former implementation used a random tie-break)
  out2 <- suppressMessages(convertUnwantedLCC(
    240L,
    ras,
    data.table::copy(aERC),
    doAssertion = FALSE
  ))
  expect_identical(out, out2)

  ## --- constrained: per-ecoregion availability ("eco_lcc" codes) ---
  eco <- ifelse(xy[, 1] < n / 2, "1", "2")
  aERC2 <- data.table(
    pixelIndex = seq_len(terra::ncell(ras)),
    initialEcoregionCode = paste0(eco, "_", formatC(v, width = 3, flag = "0"))
  )
  outC <- suppressMessages(convertUnwantedLCC(
    240L,
    ras,
    data.table::copy(aERC2),
    doAssertion = FALSE
  ))
  ## never invents an unavailable ecoregion-class combination ...
  expect_true(all(outC[!is.na(ecoregionGroup)]$ecoregionGroup %in% aERC2$initialEcoregionCode))
  ## ... so unwanted pixels in ecoregion 1 can only become the class available there ("1_210")
  eco1 <- outC[pixelIndex %in% which(v == 240L & eco == "1")]
  expect_true(all(eco1$ecoregionGroup == "1_210"))

  expect_error(convertUnwantedLCC(240L, ras, data.table::copy(aERC), method = "spread"))
})

testthat::test_that("convertUnwantedLCC method='nearestRandom' samples by local abundance", {
  ## one unwanted pixel at the centre; of its 8 neighbours exactly one is 210 (and it is
  ## the orthogonally-adjacent one, so it is also the strictly nearest), seven are 220.
  ## "nearest" must therefore always take 210, while "nearestRandom" must take it 1 time
  ## in 8 -- i.e. weighted by how much of each class the neighbourhood holds, not by which
  ## is closest, and not uniformly over the two classes.
  ## projected: on a lon/lat raster the N and S neighbours are not equidistant (a degree of
  ## latitude lengthens polewards), which would decide the "nearest" tie by ellipsoid shape
  n <- 5L
  ras <- terra::rast(
    nrows = n,
    ncols = n,
    xmin = 0,
    xmax = n,
    ymin = 0,
    ymax = n,
    crs = "EPSG:3978"
  )
  v <- rep(220L, terra::ncell(ras))
  v[terra::cellFromRowCol(ras, 2, 3)] <- 210L
  v[terra::cellFromRowCol(ras, 3, 3)] <- 240L
  terra::values(ras) <- v
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))

  expect_identical(
    suppressMessages(convertUnwantedLCC(
      240L,
      ras,
      data.table::copy(aERC),
      doAssertion = FALSE
    ))$ecoregionGroup,
    210L
  )

  set.seed(1)
  draws <- replicate(2000, {
    suppressMessages(convertUnwantedLCC(
      240L,
      ras,
      data.table::copy(aERC),
      doAssertion = FALSE,
      method = "nearestRandom"
    ))$ecoregionGroup
  })
  expect_setequal(unique(draws), c(210L, 220L)) # both are reachable
  ## 1/8, well outside both 0 (nearest) and 1/2 (uniform); ~7 SE of slack
  expect_gt(mean(draws == 210L), 0.08)
  expect_lt(mean(draws == 210L), 0.17)
})

testthat::test_that("convertUnwantedLCC method='nearestRandom' is seed-reproducible and constrained", {
  n <- 30L
  ras <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  xy <- terra::xyFromCell(ras, seq_len(terra::ncell(ras)))
  v <- ifelse(xy[, 1] < n / 2, 210L, 220L)
  d <- sqrt((xy[, 1] - n / 2)^2 + (xy[, 2] - n / 2)^2)
  v[d < 6] <- 240L
  terra::values(ras) <- v
  unwanted <- which(v == 240L)
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))

  cu <- function() {
    suppressMessages(convertUnwantedLCC(
      240L,
      ras,
      data.table::copy(aERC),
      doAssertion = FALSE,
      method = "nearestRandom"
    ))
  }
  set.seed(123)
  a <- cu()
  set.seed(123)
  b <- cu()
  set.seed(456)
  cc <- cu()
  expect_identical(a, b) # same seed reproduces exactly
  expect_false(identical(a$ecoregionGroup, cc$ecoregionGroup)) # different seed does not
  expect_setequal(a$pixelIndex, unwanted) # still resolves every unwanted pixel
  expect_false(any(a$ecoregionGroup %in% 240L))

  ## per-ecoregion availability is respected exactly as under method = "nearest"
  eco <- ifelse(xy[, 1] < n / 2, "1", "2")
  aERC2 <- data.table(
    pixelIndex = seq_len(terra::ncell(ras)),
    initialEcoregionCode = paste0(eco, "_", formatC(v, width = 3, flag = "0"))
  )
  set.seed(123)
  outC <- suppressMessages(convertUnwantedLCC(
    240L,
    ras,
    data.table::copy(aERC2),
    doAssertion = FALSE,
    method = "nearestRandom"
  ))
  expect_true(all(outC[!is.na(ecoregionGroup)]$ecoregionGroup %in% aERC2$initialEcoregionCode))
  eco1 <- outC[pixelIndex %in% which(v == 240L & eco == "1")]
  expect_true(all(eco1$ecoregionGroup == "1_210"))
  ## `newPossLCC` drops the ecoregion prefix
  expect_setequal(unique(outC$newPossLCC), c(210L, 220L))
})

testthat::test_that("convertUnwantedLCC methods agree when only one class can be chosen", {
  ## no choice to make: both methods must return the same thing
  n <- 12L
  ras <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  v <- rep(210L, terra::ncell(ras))
  v[terra::cellFromRowCol(ras, 6:7, 6:7)] <- 240L
  terra::values(ras) <- v
  aERC <- data.table(pixelIndex = seq_len(terra::ncell(ras)), initialEcoregionCode = as.integer(v))

  det <- suppressMessages(convertUnwantedLCC(
    240L,
    ras,
    data.table::copy(aERC),
    doAssertion = FALSE
  ))
  set.seed(9)
  rnd <- suppressMessages(convertUnwantedLCC(
    240L,
    ras,
    data.table::copy(aERC),
    doAssertion = FALSE,
    method = "nearestRandom"
  ))
  expect_identical(det[order(pixelIndex)], rnd[order(pixelIndex)])
})
