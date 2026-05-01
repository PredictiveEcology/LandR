## Pre-test Google Drive authentication.
##
## On CI we set GOOGLEDRIVE_AUTH to the *contents* of a service-account JSON
## key (org/repo secret). On local dev the env var is unset and `drive_auth()`
## falls back to its normal interactive flow, so this file is a no-op.
##
## testthat picks up `setup-*.R` files automatically and runs them once
## before any test in this directory.
.gdAuth <- Sys.getenv("GOOGLEDRIVE_AUTH", "")
if (nzchar(.gdAuth) && requireNamespace("googledrive", quietly = TRUE)) {
  .gdJsonPath <- tempfile(fileext = ".json")
  writeLines(.gdAuth, .gdJsonPath)
  try(googledrive::drive_auth(path = .gdJsonPath, cache = FALSE), silent = TRUE)
}
