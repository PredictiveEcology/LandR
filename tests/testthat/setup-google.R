## Pre-test Google Drive authentication.
##
## On CI we set GOOGLEDRIVE_AUTH to the *contents* of a service-account JSON
## key (org/repo secret). On local dev the env var is unset and `drive_auth()`
## falls back to its normal interactive flow, so this file is a no-op.
##
## testthat picks up `setup-*.R` files automatically and runs them once
## before any test in this directory.
.gdAuth <- Sys.getenv("GOOGLEDRIVE_AUTH", "")
cat(sprintf("[setup-google] GOOGLEDRIVE_AUTH length: %d\n", nchar(.gdAuth)))
if (nzchar(.gdAuth) && requireNamespace("googledrive", quietly = TRUE)) {
  .gdJsonPath <- tempfile(fileext = ".json")
  writeLines(.gdAuth, .gdJsonPath)
  .sz <- tryCatch(file.info(.gdJsonPath)$size, error = function(e) NA_integer_)
  cat(sprintf("[setup-google] wrote %s bytes to %s\n",
              if (is.na(.sz)) "NA" else as.character(.sz), .gdJsonPath))
  ## Quick JSON sanity-check (don't print contents — these are credentials).
  .ok <- tryCatch(jsonlite::fromJSON(.gdJsonPath), error = function(e) e)
  if (inherits(.ok, "error")) {
    cat(sprintf("[setup-google] JSON parse FAILED: %s\n", conditionMessage(.ok)))
  } else {
    cat(sprintf("[setup-google] JSON parsed OK (type=%s, has_private_key=%s)\n",
                if (is.null(.ok$type)) "?" else .ok$type,
                !is.null(.ok$private_key)))
  }
  .res <- tryCatch({
    googledrive::drive_auth(path = .gdJsonPath, cache = FALSE)
    "success"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[setup-google] drive_auth: %s\n", .res))
} else {
  cat("[setup-google] env var empty or googledrive not installed; skipping\n")
}
