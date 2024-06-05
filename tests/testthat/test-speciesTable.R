test_that("speciesTable has correct column types", {
  skip_if_offline()

  test_dir <- file.path(tempdir(), "test_speciesTable")
  dir.create(test_dir, recursive = TRUE)
  on.exit({
    reproducible::clearCache(ask = FALSE)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  expect_no_warning({
    sTraw <- getSpeciesTable(dPath = test_dir)
  })
  expect_no_error(assertSpeciesTableRaw(sTraw))

  sT <- data.table::copy(sTraw)
  expect_no_warning({
    sT <- prepSpeciesTable(sT)
  })
  expect_no_error(assertSpeciesTable(sT))
})
