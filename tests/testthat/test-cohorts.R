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
    classesToReplace = 240L, rstLCC = ras, availableERC_by_Sp = data.table::copy(aERC), doAssertion = FALSE
  ))
  expect_s3_class(out, "data.table")
  expect_named(out, c("pixelIndex", "ecoregionGroup"))
  expect_setequal(out$pixelIndex, unwanted) # every unwanted pixel handled
  expect_false(any(out$ecoregionGroup %in% 240L)) # no unwanted class remains
  expect_setequal(unique(out$ecoregionGroup), c(210L, 220L))
  ## deterministic (the former implementation used a random tie-break)
  out2 <- suppressMessages(convertUnwantedLCC(240L, ras, data.table::copy(aERC), doAssertion = FALSE))
  expect_identical(out, out2)

  ## --- constrained: per-ecoregion availability ("eco_lcc" codes) ---
  eco <- ifelse(xy[, 1] < n / 2, "1", "2")
  aERC2 <- data.table(
    pixelIndex = seq_len(terra::ncell(ras)),
    initialEcoregionCode = paste0(eco, "_", formatC(v, width = 3, flag = "0"))
  )
  outC <- suppressMessages(convertUnwantedLCC(240L, ras, data.table::copy(aERC2), doAssertion = FALSE))
  ## never invents an unavailable ecoregion-class combination ...
  expect_true(all(outC[!is.na(ecoregionGroup)]$ecoregionGroup %in% aERC2$initialEcoregionCode))
  ## ... so unwanted pixels in ecoregion 1 can only become the class available there ("1_210")
  eco1 <- outC[pixelIndex %in% which(v == 240L & eco == "1")]
  expect_true(all(eco1$ecoregionGroup == "1_210"))
})
