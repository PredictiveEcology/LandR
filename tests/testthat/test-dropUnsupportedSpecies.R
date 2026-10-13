## minSpeciesEcoregionShare: a species rare in an ecoregion is removed from all of its pixels there.

test_that("dropUnsupportedSpecies removes species below the share, per ecoregion (not per lcc)", {
  ## ecoregion A: 20 pixels over two land-cover classes; Thuj_pli in 1 of them (5%), Pseu_men in all.
  ## ecoregion B: 10 pixels; Thuj_pli in 5 (50%).
  cd <- data.table::data.table(
    pixelIndex = c(1:20, 1, 21:30, 21:25),
    speciesCode = c(rep("Pseu_men", 20), "Thuj_pli", rep("Pseu_men", 10), rep("Thuj_pli", 5)),
    initialEcoregionCode = c(rep(c("A_210", "A_220"), each = 10), "A_210",
                             rep("B_210", 10), rep("B_210", 5))
  )
  expect_identical(LandR:::dropUnsupportedSpecies(cd, 0), cd)
  out <- LandR:::dropUnsupportedSpecies(cd, 0.10)
  expect_identical(nrow(out[speciesCode == "Thuj_pli" & startsWith(as.character(initialEcoregionCode), "A")]), 0L)
  expect_identical(nrow(out[speciesCode == "Thuj_pli" & startsWith(as.character(initialEcoregionCode), "B")]), 5L)
  expect_identical(nrow(out[speciesCode == "Pseu_men"]), 30L)
  ## at 5% exactly, Thuj_pli stays in A
  expect_identical(nrow(LandR:::dropUnsupportedSpecies(cd, 0.05)), nrow(cd))
})

test_that("ecoregion labels with underscores keep everything before the land-cover suffix", {
  cd <- data.table::data.table(pixelIndex = c(1L, 2L, 3L, 3L), speciesCode = c("A", "A", "A", "B"),
                               initialEcoregionCode = c("BC_CWH_210", "BC_CWH_220", "BC_MH_210", "BC_MH_210"))
  out <- LandR:::dropUnsupportedSpecies(cd, 0.6)
  ## BC_CWH: A in 2 of 2 pixels (kept). BC_MH: A in 1 of 1 pixel, B in 1 of 1 (both kept)
  expect_identical(nrow(out), 4L)
})
