testthat::test_that("sppColors works", {
  ## confirm standardized
  sppSet1 <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Pice_mar")]
  SppColors1 <- sppColors(sppEquiv = sppSet1, sppEquivCol = "LandR", newVals = "Mixed")
  expect_true(SppColors1[["Pice_mar"]] == "#00479E") ## standardized Pice_mar is blue

  ## confirm standardized even with popu_spp
  sppSet2 <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Pice_mar")]
  sppSet2[LandR %in% c("Popu_tre", "Popu_bal"), LandR := "Popu_spp"]
  SppColors2 <- sppColors(sppSet2, "LandR", "Mixed")
  expect_true(SppColors2[["Popu_spp"]] == "#FFACFD")

  ## confirm it doesn't use standardized
  sppSet3 <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Arbu_men")]
  SppColors3 <- sppColors(sppSet3, "LandR", "Mixed", palette = "Accent")
  expect_true(SppColors3["Popu_tre"] != "#FF00B6")

  ## confirm the fact that many rows with fewer unique values still uses standardized
  sppSet4 <- rbind(sppSet2, sppEquivalencies_CA[LandR %in% c("Pseu_men"),])
  SppColors4 <- sppColors(sppSet4, "LandR", newVals = "Mixed")
  expect_true(SppColors4["Pseu_men"] == "#720055")
})

testthat::test_that("sppColors falls back to the palette when colours are not distinct", {
  ## The "enough distinct colours" test read
  ## `length(unique(sppEquiv[[sppEquivCol]] <= length(unique(sppEquiv$colorHex))))`,
  ## which compares names to a number and then takes the length of the result: 1 or 2,
  ## both truthy. Two species sharing one colorHex therefore both got that colour.
  sppSet <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Pice_mar")]
  sppSet[, colorHex := "#00479E"] ## standardized Pice_mar blue, for both species

  sppCols <- sppColors(sppSet, "LandR", newVals = "Mixed")

  expect_false(any(duplicated(sppCols[c("Popu_tre", "Pice_mar")])))
  expect_true(sppCols[["Popu_tre"]] != "#00479E")
})

testthat::test_that("sppColors adds a newVals colour only when one is asked for", {
  ## `length(newVals == 1)` was a misplaced bracket; with a length-2 `newVals` it would have
  ## named a single gray with two names, but the check above keeps that case off this path.
  sppSet <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Pice_mar")]

  ## three Popu_bal rows (the species and two varieties) collapse to one name
  expect_identical(names(sppColors(sppSet, "LandR")), unique(sppSet$LandR))
  expect_identical(
    names(sppColors(sppSet, "LandR", newVals = "Mixed")), c(unique(sppSet$LandR), "Mixed")
  )
})
