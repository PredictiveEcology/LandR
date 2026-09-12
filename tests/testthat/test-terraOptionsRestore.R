## `terraOptions()` and `rasterOptions()` are process-wide. A function that sets one
## and returns without restoring it silently reconfigures the rest of the session --
## which is what happened to a long-lived worker that had set memmax = 4 and later
## found it at 1, set by makePickellStack().
##
## A behavioural test would need a Pickell raster fixture, so this asserts the
## invariant at the source level: every function in this package that *sets* a terra
## or raster option also registers a restore.

test_that("functions that set terra/raster options restore them on exit", {
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "source not available (installed package)")

  files <- list.files(rdir, pattern = "[.]R$", full.names = TRUE)
  offenders <- character()

  for (f in files) {
    src <- readLines(f, warn = FALSE)
    code <- grep("^\\s*#", src, invert = TRUE, value = TRUE)
    ## a *setting* call passes a named argument; reading passes none, or print/default
    sets <- grep("(terraOptions|rasterOptions)\\s*\\([^)]*=", code, value = TRUE)
    sets <- grep("print\\s*=|default\\s*=", sets, invert = TRUE, value = TRUE)
    if (!length(sets)) next
    ## and the same file must register a restore
    restores <- grep("on\\.exit\\(.*(terraOptions|rasterOptions)", code, value = TRUE)
    if (!length(restores)) offenders <- c(offenders, basename(f))
  }

  expect_equal(offenders, character(0),
               info = paste("these set a terra/raster option with no on.exit restore:",
                            paste(offenders, collapse = ", ")))
})
