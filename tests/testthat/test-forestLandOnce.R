## fireSense calls prepInputs_SCANFI_LCC_FAO() once per dataYear on one study area. The
## forest-land inputs -- the FAO layer and each forestLandYears land cover -- are the same for
## every one of those calls, so they must be prepared once, not once per dataYear. In production
## each FAO preparation cost ~12 min, five times per study area. It must hold with reproducible's
## Cache on and off: SpaDES.core's `spades.useCache = "eventsOnly"` sets
## options(reproducible.useCache = FALSE) for the run, and fireSense_dataPrepFit wraps each
## dataYear's call in its own Cache(userTags = c("makeFireSenseLCC", dy)), which is then FALSE
## and, being the outer call, would override an inner useCache = TRUE.
##
## No network: `prepInputs` is replaced by one that counts its calls and postProcesses a small
## local file with the real reproducible::postProcessTo().

## Land-cover and FAO values on the 4 x 2 target window (row-major).
##   cell:     1    2    3    4    5    6    7    8
.flBase <- c( 50, 100,  20, 210,  33,  50,  50,  50)
.flFAO  <- c(  2,   1,   0,   1,   0,   1,   1,   0)
.flYear <- function(y) {
  v <- .flBase
  if (y == 1990) v[2] <- NA  ## no data in 1990 only: must not hide cell 2's FAO in other years
  if (y == 1985) v[3] <- 210 ## a lake edge that was treed in 1985
  if (y == 2015) v[c(6, 8)] <- c(220, 210)
  v
}

## Sources are larger than the target, so the crop is real; the target cells are the middle
## 4 x 2 block of a 6 x 4 grid.
.flWriteSource <- function(vals, dir, name) {
  r <- terra::rast(nrows = 4, ncols = 6, xmin = 0, xmax = 6000, ymin = 0, ymax = 4000,
                   crs = "EPSG:3978", vals = 50)
  m <- matrix(terra::values(r, mat = FALSE), nrow = 4, byrow = TRUE)
  m[2:3, 2:5] <- matrix(vals, nrow = 2, byrow = TRUE)
  terra::values(r) <- as.vector(t(m))
  f <- file.path(dir, name)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

## mode: "optionTRUE"  -- Cache on, each dataYear's call inside the module's per-year Cache()
##       "optionFALSE" -- Cache off (eventsOnly), called directly
##       "nestedFALSE" -- Cache off, inside the module's per-year Cache() with userTags, which is
##                        how fireSense_dataPrepFit.R runs it under eventsOnly
.flRunDataYears <- function(dataYears, mode) {
  withr::local_package("terra")
  srcDir <- withr::local_tempdir("flSrc_")
  dPath <- withr::local_tempdir("flInputs_")
  withr::local_options(list(
    reproducible.useCache = mode == "optionTRUE",
    reproducible.cachePath = withr::local_tempdir("flCache_"),
    reproducible.destinationPath = dPath,
    reproducible.verbose = -1
  ))

  scanfiYears <- LandR:::.scanfi_v2_years
  src <- lapply(scanfiYears, function(y) {
    .flWriteSource(.flYear(y), srcDir, LandR:::.scanfiLCCFAOSource(y, "V2")$targetFile)
  })
  names(src) <- vapply(scanfiYears, function(y) LandR:::.scanfiLCCFAOSource(y, "V2")$targetFile,
                       character(1))
  src[["CA_FAO_forest_2022.zip"]] <- .flWriteSource(.flFAO, srcDir, "CA_FAO_forest_2022.tif")

  counts <- new.env()
  testthat::local_mocked_bindings(
    prepInputs = function(url = NULL, targetFile = NULL, destinationPath = NULL,
                          method = NULL, overwrite = NULL, writeTo = NULL, ...) {
      key <- if (is.null(targetFile)) basename(url) else targetFile
      counts[[key]] <- sum(counts[[key]], 1L)
      reproducible::postProcessTo(terra::rast(src[[key]]), ...)
    },
    .package = "LandR"
  )

  ## makeFireSenseLCC()'s `to` and `maskTo`: the target window, with one masked-out cell (7)
  to <- terra::rast(nrows = 2, ncols = 4, xmin = 1000, xmax = 5000, ymin = 1000, ymax = 3000,
                    crs = "EPSG:3978", vals = c(1, 1, 1, 1, 1, 1, NA, 1))

  outs <- lapply(dataYears, function(y) {
    ## as fireSenseUtils::makeFireSenseLCC() calls it
    prep <- function(y) {
      prepInputs_SCANFI_LCC_FAO(year = y, disturbedCode = 240, overwrite = TRUE,
                                destinationPath = dPath, cropTo = to, writeTo = NULL,
                                maskTo = to)
    }
    out <- if (mode == "optionFALSE") {
      prep(y)
    } else {
      reproducible::Cache(prep(y), userTags = c("makeFireSenseLCC", y))
    }
    terra::values(out, mat = FALSE)
  })
  names(outs) <- dataYears
  list(values = outs, counts = unlist(as.list(counts)))
}

test_that("forest-land inputs are prepared once across dataYears, whatever the Cache setting", {
  skip_if_not_installed("withr")
  dataYears <- c(1990, 2015, 2020)
  forestLandYears <- LandR:::.defaultForestLandYears(LandR:::.scanfi_v2_years)
  expect_equal(forestLandYears, c(1985, 1995, 2005, 2015, 2025))

  ## Hand-derived from the rule, identical to what the code produced before the reuse:
  ## forest land = FAO 1|2, or treed in a forestLandYear other than the dataYear itself.
  expected <- list(
    "1990" = c(240,  NA, 240, 210, 33, 240, NA, 240),
    "2015" = c(240, 240, 240, 210, 33, 220, NA, 210),
    "2020" = c(240, 240, 240, 210, 33, 240, NA, 240)
  )

  for (mode in c("optionTRUE", "optionFALSE", "nestedFALSE")) {
    res <- .flRunDataYears(dataYears, mode = mode)
    info <- mode

    ## the FAO layer: once, not once per dataYear
    expect_equal(res$counts[["CA_FAO_forest_2022.zip"]], 1L, info = info)
    ## each forest-land year: once as a forest-land input, plus once more if it is also a
    ## dataYear (its own land cover, which is not a forest-land preparation)
    for (y in forestLandYears) {
      tf <- LandR:::.scanfiLCCFAOSource(y, "V2")$targetFile
      expect_equal(res$counts[[tf]], 1L + (y %in% dataYears), info = paste(info, "; year", y))
    }

    expect_equal(res$values, expected, info = info)
  }
})
