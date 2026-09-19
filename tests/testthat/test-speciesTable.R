testthat::test_that("speciesTable has correct column types", {
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_offline()

  test_dir <- withr::local_tempdir("test_speciesTable")

  withr::defer({
    reproducible::clearCache(ask = FALSE)
  })

  testthat::expect_no_warning({
    sTraw <- getSpeciesTable(dPath = test_dir)
  })
  testthat::expect_no_error(assertSpeciesTableRaw(sTraw))

  sT <- data.table::copy(sTraw)
  testthat::expect_no_warning({
    sT <- prepSpeciesTable(sT)
  })
  testthat::expect_no_error(assertSpeciesTable(sT))
})
