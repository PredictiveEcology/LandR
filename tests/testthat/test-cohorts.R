testthat::test_that("adjustAgeToLongevity works", {
  # test that it returns expected results
  cd <- data.table(
    speciesCode = c("A", "A", "A", "B", "B"),
    age = c(10, 200, 300, 200, 100)
  )
  traits <- data.table(
    speciesCode = c("A", "B"),
    longevity = c(200, 250)
  )
  expect_equal(
    suppressMessages(adjustAgeToLongevity(cd, traits, 0.9)),
    data.table(
      speciesCode = c("A", "A", "A", "B", "B"),
      age = c(10, 177, 180, 200, 100)
    )
  )
  expect_equal(
    suppressMessages(adjustAgeToLongevity(cd, traits, 0.7)),
    data.table(
      speciesCode = c("A", "A", "A", "B", "B"),
      age = c(10, 140, 140, 173, 100)
    )
  )
  expect_equal(
    suppressMessages(adjustAgeToLongevity(cd, traits, 0.5)),
    data.table(
      speciesCode = c("A", "A", "A", "B", "B"),
      age = c(10, 100, 100, 125, 100)
    )
  )

  # test incorrect inputs
  expect_error(adjustAgeToLongevity(cd, traits, 1.1))
  expect_error(adjustAgeToLongevity(cd, traits, NA))
  expect_error(adjustAgeToLongevity(cd, traits, 0.1))
  expect_error(adjustAgeToLongevity(cd, traits, -0.9))
  colnames(traits) <- c("spp", "longevity")
  expect_error(adjustAgeToLongevity(cd, traits, 0.5))
})
