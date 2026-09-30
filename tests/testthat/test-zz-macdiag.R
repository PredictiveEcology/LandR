## TEMPORARY diagnostic (PR #257, DO NOT MERGE): print the inputs and per-step
## trace of spiralLoopCpp for one tiny case, so macOS and Linux logs can be diffed.
## Fails on purpose so the text appears in the R CMD check log.
test_that("macdiag: spiralLoopCpp inputs and trace", {
  fix <- makeLANDISDispFixture(size = "tiny", fixtureSeed = 11L, successionTimestep = 1L)
  ns <- asNamespace("LandR")
  orig <- get("spiralLoopCpp", envir = ns)
  captured <- NULL
  unlockBinding("spiralLoopCpp", ns)
  assign("spiralLoopCpp", function(...) { captured <<- list(...); orig(...) }, envir = ns)
  on.exit({ assign("spiralLoopCpp", orig, envir = ns); lockBinding("spiralLoopCpp", ns) })

  fp <- function(x) {
    if (is.matrix(x)) x <- as.vector(x)
    if (is.logical(x) && length(x) == 1) return(as.character(x))
    v <- as.numeric(x); ok <- !is.na(v)
    sprintf("len=%d nNA=%d sum=%.17g wsum=%.17g head=%s", length(v), sum(!ok),
            sum(v[ok]), sum(v[ok] * seq_along(v)[ok]),
            paste(format(head(v, 6), digits = 17), collapse = ","))
  }
  lines <- c(sprintf("platform=%s RNGkind=%s", R.version$platform, paste(RNGkind(), collapse = "/")))
  set.seed(42); lines <- c(lines, paste("unif:", paste(format(runif(4), digits = 17), collapse = ",")))
  set.seed(42); lines <- c(lines, paste("runifC:", paste(format(SpaDES.tools::runifC(4), digits = 17), collapse = ",")))
  withr::local_options(LandR.LANDISDisp.debug = TRUE)
  trace <- utils::capture.output(out <- runLANDISDispOnFixture(fix, runSeed = 42L, useCpp = TRUE))
  lines <- c(lines, vapply(names(captured), function(n) paste0(n, ": ", fp(captured[[n]])), ""))
  lines <- c(lines, sprintf("nOut=%d", NROW(out)), "trace:", head(grep("^\\[cpp\\]", trace, value = TRUE), 120))
  fail(paste(lines, collapse = "\n"))
})
