test_that("no spelling errors found", {
  skip_on_os(c("macOS", "solaris", "windows"))

  if (isTRUE(Sys.info()[["sysname"]] == "Linux") && isTRUE(Sys.getenv("TRAVIS" == "true")))
    if (requireNamespace("spelling", quietly = TRUE))
      spelling::spell_check_test(vignettes = TRUE, error = FALSE, skip_on_cran = TRUE)
})
