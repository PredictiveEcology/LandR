withr::local_options(list(
  warnPartialMatchArgs = TRUE,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE
))

library(testthat)
library(LandR)

test_check("LandR")
