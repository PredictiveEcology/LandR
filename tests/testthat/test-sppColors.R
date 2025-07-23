testthat::test_that("sppColors works", {

  #confirm standardized
  sppSet1 <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Pice_mar")]
  SppColors1 <- sppColors(sppEquiv = sppSet1, sppEquivCol = "LandR", newVals = "Mixed")
  expect_true(SppColors1[["Pice_mar"]] == "#00479E") #standardized Pice_mar is blue
  #confirm standardized even with popu_spp
  sppSet2 <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Pice_mar")]
  sppSet2[LandR %in% c("Popu_tre", "Popu_bal"), LandR := "Popu_spp"]
  SppColors2 <- sppColors(sppSet2, "LandR", "Mixed")
  expect_true(SppColors2[["Popu_spp"]] == "#FFACFD")

  #confirm it doesn't use standardized
  sppSet3 <- sppEquivalencies_CA[LandR %in% c("Popu_tre", "Popu_bal", "Arbu_men")]
  SppColors3 <- sppColors(sppSet3, "LandR", "Mixed", palette = "Accent")
  expect_true(SppColors3["Popu_tre"] != "#FF00B6")
  #confirm the fact that many rows with fewer unique values still uses standardized
  sppSet4 <- rbind(sppSet2, sppEquivalencies_CA[LandR %in% c("Pseu_men"),])
  SppColors4 <- sppColors(sppSet4, "LandR", newVals = "Mixed")
  expect_true(SppColors4["Pseu_men"] == "#720055")
})
