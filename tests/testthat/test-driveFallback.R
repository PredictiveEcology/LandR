## The Google Drive fallback in prepSpeciesLayers_SCANFI()/_KNN(). On
## 2026-09-07 a primary source that was itself a Drive URL was pinged with
## RCurl::url.exists(); under packet loss the ping failed, the fallback ran,
## found no folder of the expected name, and drive_ls(character(0)) killed the
## job two hours in with only "Parent specified via `path` is invalid: Does not
## exist" to show for it.

test_that(".isGoogleDriveUrl recognises Drive links and nothing else", {
  expect_true(LandR:::.isGoogleDriveUrl("https://drive.google.com/file/d/15T4HIFeqzwp0TuOuxmYoexuXdLFnCZBi"))
  expect_true(LandR:::.isGoogleDriveUrl("https://drive.google.com/drive/folders/1zuHRIDWIzKyWcvcgG-p3bXA0Rek3xmaQ"))
  expect_false(LandR:::.isGoogleDriveUrl("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/x/"))
  expect_false(LandR:::.isGoogleDriveUrl(NULL))
  expect_false(LandR:::.isGoogleDriveUrl(NA_character_))
  expect_false(LandR:::.isGoogleDriveUrl(c("https://drive.google.com/a", "https://drive.google.com/b")))
})

test_that(".driveFolderLink fails loudly, naming the folder and where it looked", {
  dt <- data.table::data.table(name = c("SCANFIForestAttributes_2010", "readme.txt"),
                               id = c("1aaa", "1bbb"))
  expect_error(LandR:::.driveFolderLink(dt, "SCANFIForestAttributes_2020", "https://drive.google.com/drive/folders/X"),
               "no folder named 'SCANFIForestAttributes_2020'.*folders/X.*SCANFIForestAttributes_2010")
})

test_that(".driveFolderLink returns a folder URL for a present folder, without a network call", {
  dt <- data.table::data.table(name = c("readme.txt", "SCANFIForestAttributes_2020"),
                               id = c("1bbb", "1ccc"))
  expect_equal(LandR:::.driveFolderLink(dt, "SCANFIForestAttributes_2020", "https://x"),
               "https://drive.google.com/drive/folders/1ccc")
})

test_that("neither prepSpeciesLayers_* references an unbound `name`/`id` any more", {
  skip_if_not_installed("codetools")
  for (f in list(prepSpeciesLayers_SCANFI, prepSpeciesLayers_KNN)) {
    vars <- codetools::findGlobals(f, merge = FALSE)$variables
    expect_false(any(c("name", "id") %in% vars))
  }
})
