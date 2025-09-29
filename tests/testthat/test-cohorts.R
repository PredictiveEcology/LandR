testthat::test_that("adjustAgeToLongevity works", {
  # test that it returns expected results
  cd <- data.table(
    speciesCode = c("A", "A", "B"),
    age = c(10, 20, 20)
  )
  traits <- data.table(
    speciesCode = c("A", "B"),
    longevity = c(19, 30)
  )
  expect_equal(
    adjustAgeToLongevity(cd, traits, 1),
    data.table(
      speciesCode = c("A", "A", "B"),
      age = c(10, 19, 20)
    )
  )
  expect_equal(
    adjustAgeToLongevity(cd, traits, 0.7),
    data.table(
      speciesCode = c("A", "A", "B"),
      age = c(10, round(19 * 0.7), 20)
    )
  )
  expect_equal(
    adjustAgeToLongevity(cd, traits, 0.5),
    data.table(
      speciesCode = c("A", "A", "B"),
      age = c(10, round(19 * 0.5), 15)
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
