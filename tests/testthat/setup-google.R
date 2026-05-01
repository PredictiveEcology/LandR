## Pre-test Google Drive authentication.
##
## On CI, the workflow's "Stage Google service-account credentials" step
## writes the org secret to a JSON file and sets `GOOGLEDRIVE_AUTH_FILE` to
## the path. We avoid putting the JSON content directly in an env var because
## multi-line env vars don't survive the GitHub Actions runner setup on
## Windows (PowerShell can't set them, so the var arrives empty).
##
## On local dev neither var is set, and `drive_auth()` falls back to its
## normal interactive flow — this file is a no-op.
##
## testthat picks up `setup-*.R` files automatically and runs them once
## (per worker, in parallel mode) before any test in this directory.
.gdAuthFile <- Sys.getenv("GOOGLEDRIVE_AUTH_FILE", "")
cat(sprintf("[setup-google] GOOGLEDRIVE_AUTH_FILE: '%s'\n", .gdAuthFile))
if (nzchar(.gdAuthFile) && file.exists(.gdAuthFile) &&
    requireNamespace("googledrive", quietly = TRUE)) {
  cat(sprintf("[setup-google] file size: %d bytes\n",
              file.info(.gdAuthFile)$size))
  .res <- tryCatch({
    googledrive::drive_auth(path = .gdAuthFile, cache = FALSE)
    "success"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[setup-google] drive_auth: %s\n", .res))
} else {
  cat("[setup-google] no SA credentials staged; skipping\n")
}
