drive_error <- function(status = "Client error: (404) Not Found") {
  simpleError(paste0(
    "Could not access the Google Drive resource:\n",
    "  https://drive.google.com/file/d/1nXPS3bpFUESYieNfXO25OKlZJEgqtRnD\n",
    "  (open it in a browser to confirm it exists and is shared 'Anyone with the link')\n",
    "  ", status
  ))
}

test_that(".isGoogleDriveAccessError recognizes access failures but not other errors", {
  expect_true(LandR:::.isGoogleDriveAccessError(drive_error()))
  expect_true(LandR:::.isGoogleDriveAccessError(simpleError("Client error: (404) Not Found")))
  expect_true(LandR:::.isGoogleDriveAccessError(simpleError("HTTP 403 Forbidden")))
  expect_true(LandR:::.isGoogleDriveAccessError(simpleError("Permission denied")))

  ## must not swallow unrelated failures and mislabel them as an access problem
  expect_false(LandR:::.isGoogleDriveAccessError(simpleError("cannot open connection")))
  expect_false(LandR:::.isGoogleDriveAccessError(simpleError("subscript out of bounds")))
  expect_false(LandR:::.isGoogleDriveAccessError(simpleError("non-numeric argument")))
})

test_that(".withSCANFIAccess explains an access failure and keeps the original error", {
  err <- tryCatch(
    LandR:::.withSCANFIAccess(
      stop(drive_error()),
      what = "the SCANFI stand age map",
      dataYear = 2020,
      dataVersion = "V2",
      urlArg = "ageURL",
      alternatives = c("KNN", "NTEMS")
    ),
    error = function(e) conditionMessage(e)
  )

  ## the guidance is wrapped for the console, so match on content not layout
  flat <- gsub("[[:space:]]+", " ", err)

  expect_match(flat, "Could not download the SCANFI stand age map (V2 2020)", fixed = TRUE)
  expect_match(flat, "permissions problem rather than a broken link", fixed = TRUE)

  ## the three ways out
  expect_match(flat, "https://opendata.nfis.org/", fixed = TRUE)
  expect_match(flat, "googledrive::drive_auth()", fixed = TRUE)
  expect_match(flat, "`ageURL`", fixed = TRUE)
  expect_match(flat, 'dataSource = "KNN" or dataSource = "NTEMS"', fixed = TRUE)

  ## the underlying error is still reported, not swallowed
  expect_match(flat, "Client error: (404) Not Found", fixed = TRUE)
})

test_that(".withSCANFIAccess omits the bullets that do not apply", {
  err <- tryCatch(
    LandR:::.withSCANFIAccess(stop(drive_error()), what = "the SCANFI species layers"),
    error = function(e) conditionMessage(e)
  )

  expect_match(err, "https://opendata.nfis.org/", fixed = TRUE)
  expect_false(grepl("dataSource =", err, fixed = TRUE))
  expect_false(grepl("pass it via", err, fixed = TRUE))
})

test_that(".withSCANFIAccess leaves unrelated errors alone", {
  expect_error(
    LandR:::.withSCANFIAccess(stop("something else broke"), what = "the SCANFI stand age map"),
    "something else broke"
  )
  ## and does not add SCANFI guidance to them
  err <- tryCatch(
    LandR:::.withSCANFIAccess(stop("something else broke"), what = "the SCANFI stand age map"),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("opendata.nfis.org", err, fixed = TRUE))
})

test_that(".withSCANFIAccess is transparent when it succeeds or is disabled", {
  expect_identical(LandR:::.withSCANFIAccess(42L, what = "x"), 42L)

  ## `enabled = FALSE` must not intercept: used by functions serving several sources
  expect_error(
    LandR:::.withSCANFIAccess(stop(drive_error()), what = "x", enabled = FALSE),
    "Could not access the Google Drive resource"
  )
  expect_identical(LandR:::.withSCANFIAccess(7L, what = "x", enabled = FALSE), 7L)
})
