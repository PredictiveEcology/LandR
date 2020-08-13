test_that("no spelling errors found", {
  skip_on_cran()
  skip_on_os(c("macOS", "solaris", "windows"))

  if (isTRUE(Sys.info()[["sysname"]] == "Linux"))
    if (requireNamespace("spelling", quietly = TRUE))
      spelling::spell_check_package("../../../LandR")
})
