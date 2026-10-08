## A draw is compared against its own species' ward probability, not against
## the largest ward probability of the species that hit on the previous step.

## 1 x 11 raster, 100 m cells. The receiver (pg 2) is in column 6; pg 9 is
## filler with no source. Species A (short) has its source 4 cells away and B
## (long) has its source 5 cells away, so A is the only species with a source
## on the 400 m ring and B is the only one on the 500 m ring. With
## successionTimestep = 10 A's ward probability at 400 m is ~0.04 and B's at
## 500 m is ~0.34.
.wardScreenFixture <- function() {
  pgm <- terra::rast(xmin = 0, xmax = 1100, ymin = 0, ymax = 100,
                     resolution = c(100, 100), vals = c(1L, rep(9L, 4), 2L, rep(9L, 3), 3L, 9L))
  spTab <- data.table::data.table(
    species = c("A", "B"), speciesCode = 1:2,
    seeddistance_eff = c(50, 500), seeddistance_max = c(600, 2000)
  )
  list(
    pgm = pgm, spTab = spTab,
    dtSrc = data.table::data.table(pixelGroup = c(3L, 1L), speciesCode = c(1L, 2L)),
    dtRcv = data.table::data.table(pixelGroup = c(2L, 2L), speciesCode = 1:2)
  )
}

test_that("long-distance species is screened by its own ward probability, not the previous step's", {
  skip_if_not_installed("terra")
  f <- .wardScreenFixture()
  ts <- 10
  pA <- 1 - (1 - Ward(400, 100, 50, 600, k = 0.95, b = 0.01))^ts
  pB <- 1 - (1 - Ward(500, 100, 500, 2000, k = 0.95, b = 0.01))^ts
  expect_lt(pA, 0.05)
  expect_gt(pB, 0.3)

  nSeeds <- 1500L
  for (useCpp in c(FALSE, TRUE)) {
    gotB <- vapply(seq_len(nSeeds), function(s) {
      set.seed(s)
      out <- LANDISDisp(f$dtSrc, f$dtRcv, f$pgm, f$spTab, successionTimestep = ts,
                        verbose = 0, useCpp = useCpp)
      2L %in% out$speciesCode
    }, logical(1))
    ## binomial 4-sigma band around the species' own probability
    tol <- 4 * sqrt(pB * (1 - pB) / nSeeds)
    expect_lt(abs(mean(gotB) - pB), tol,
              label = paste("|B success rate - pB|, useCpp =", useCpp))
  }
})
