## Wrapper: load Google auth, then run devtools::test() sequentially so the
## auth state propagates (testthat's parallel mode would spawn fresh workers
## that don't inherit the auth tokens).
eval(parse(file = "~/googledriveAuthentication.R")) |> options()
cat("[auth] drive_user:\n"); print(googledrive::drive_user())

Sys.setenv(LANDR_SLOW_TESTS = "1")
## Force sequential testthat. Parallel workers are fresh callr processes that
## don't inherit the parent's options() — including the gargle OAuth state.
## Sys.setenv is propagated, options() is not. We override the package's
## Config/testthat/parallel: true via the env var.
Sys.setenv(TESTTHAT_PARALLEL = "false")
testthat::test_local(stop_on_failure = FALSE,
                     reporter = c("summary", "fail"))
