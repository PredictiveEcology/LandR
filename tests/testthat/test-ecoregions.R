## makeEcoregionMap() must read ecoregionMap's raw cell values as mapcodes, not the
## SpatRaster's active category. as.data.table(<SpatRaster>) returns the ACTIVE category's
## labels, and terra changes which category is active on a file round trip (writeRaster()
## then rast(), which Cache() and terraOptions(todisk = TRUE) both do): the first text
## column becomes active instead of "mapcode". Build the categorical raster the same way
## ecoregionProducer() does (>= 10 codes, several text category columns, so label order
## matters), force it through a file round trip, then check every pixel's ecoregionGroup
## still matches the ecoregion_lcc of its original raw mapcode.

test_that("makeEcoregionMap() recovers the correct ecoregionGroup after a raster file round trip", {
  nCodes <- 20L
  levs <- data.frame(
    ID = seq_len(nCodes),
    mapcode = seq_len(nCodes),
    ecoregion = sprintf("er%02d", rep(1:5, each = 4)),
    landcover = sprintf("lc%d", rep(1:4, times = 5)),
    ecoregion_lcc = paste0(sprintf("er%02d", rep(1:5, each = 4)), "_", sprintf("lc%d", rep(1:4, times = 5))),
    stringsAsFactors = TRUE
  )

  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, crs = "EPSG:3978")
  origRaw <- rep(seq_len(nCodes), length.out = terra::ncell(r))
  terra::values(r) <- origRaw
  levels(r) <- levs

  ## sanity check on the premise: a file round trip changes the active category away from
  ## "mapcode" (terra 1.9.46 makes the first text column, "ecoregion", active instead)
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f)
  rRead <- terra::rast(f)
  expect_false(identical(terra::activeCat(rRead), terra::activeCat(r)))

  ecoregionTable <- data.table::as.data.table(levs)[, .(active = "yes", mapcode, ecoregion, landcover, ecoregion_lcc)]
  data.table::setnames(ecoregionTable, "ecoregion_lcc", "ecoregionGroup")
  ecoregionTable[, ecoregionName := ecoregion]

  ecoregionFiles <- list(ecoregionMap = rRead, ecoregion = ecoregionTable)
  pixelCohortData <- data.table::data.table(ecoregionGroup = unique(ecoregionTable$ecoregionGroup))

  result <- makeEcoregionMap(ecoregionFiles, pixelCohortData)

  expectedGroup <- as.character(levs$ecoregion_lcc[match(origRaw, levs$mapcode)])

  actualRaw <- terra::values(result, mat = FALSE)
  cats1 <- terra::cats(result)[[1]]
  actualGroup <- as.character(cats1$ecoregionGroup[match(actualRaw, cats1$ID)])

  expect_identical(actualGroup, expectedGroup)
})

## makeEcoregionMap() writes each pixel's ecoregionGroup as its factor level index (alphabetical
## order), so the category table's IDs must be those same indices. It numbered the table's rows
## instead, in the order the groups first appear in `ecoregionFiles$ecoregion`. prepEcoregions()
## takes that order from the ecoregion polygons (e.g. "02_*" rows before "01_*"), so a reader that
## matches cell values to `ID`, as pemisc::factorValues2() does, got another group's labels. The
## test above lists its groups alphabetically, so the two orders agree there.
test_that("makeEcoregionMap() category IDs match the cell values when groups are not alphabetical", {
  ## ecoregionProducer() numbers mapcodes by ecoregion_lcc, alphabetically
  groups <- c("01_210", "01_220", "02_210", "02_220", "03_210", "03_220")
  ecoNames <- c("01" = "138", "02" = "139", "03" = "140")
  ## rows in the order prepEcoregions() leaves them: ecoregion 02, then 03, then 01
  ord <- c(3:6, 1:2)
  eco <- substr(groups[ord], 1, 2)
  ecoregionTable <- data.table::data.table(
    active = "yes",
    mapcode = ord,
    ecoregion = factor(eco),
    landcover = factor(substr(groups[ord], 4, 6)),
    ecoregionGroup = factor(groups[ord]),
    ecoregionName = factor(unname(ecoNames[eco]))
  )

  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6, ymin = 0, ymax = 6, crs = "EPSG:3978")
  origRaw <- rep(seq_along(groups), length.out = terra::ncell(r))
  terra::values(r) <- origRaw
  ecoregionFiles <- list(ecoregionMap = r, ecoregion = ecoregionTable)
  ## a group with no cohorts is dropped, so the cell values are not the mapcodes either
  pixelCohortData <- data.table::data.table(ecoregionGroup = factor(setdiff(groups, "02_220")))

  result <- makeEcoregionMap(ecoregionFiles, pixelCohortData)

  vals <- terra::values(result, mat = FALSE)
  expected <- groups[origRaw]
  expected[expected == "02_220"] <- NA
  expect_identical(as.character(pemisc::factorValues2(result, vals, att = "ecoregionGroup")), expected)
  expect_identical(as.character(pemisc::factorValues2(result, vals, att = "ecoregionName")),
                   unname(ecoNames[substr(expected, 1, 2)]))

  cats1 <- terra::cats(result)[[1]]
  expect_setequal(cats1$ID, unique(na.omit(vals)))
  expect_identical(anyDuplicated(cats1$ID), 0L)
})

test_that("makeEcoregionMap() stops when an ecoregionGroup has two ecoregionNames", {
  ecoregionTable <- data.table::data.table(
    active = "yes",
    mapcode = 1:3,
    ecoregion = factor(c("01", "01", "02")),
    landcover = factor(c("210", "210", "210")),
    ecoregionGroup = factor(c("01_210", "01_210", "02_210")),
    ecoregionName = factor(c("138", "139", "140"))
  )
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3, crs = "EPSG:3978")
  terra::values(r) <- rep(1:3, length.out = terra::ncell(r))
  ecoregionFiles <- list(ecoregionMap = r, ecoregion = ecoregionTable)
  pixelCohortData <- data.table::data.table(ecoregionGroup = factor(c("01_210", "02_210")))

  expect_error(makeEcoregionMap(ecoregionFiles, pixelCohortData), "more than one.*01_210")
})

## ecoregionProducer() used raster::factorValues(), which on a SpatRaster returns only the ACTIVE
## category. A file round trip makes the first text column active, so a category table with a text
## column ahead of `ecoregionName` returned the wrong labels and the join on "ecoregionName" failed.
test_that("ecoregionProducer() reads ecoregionName by cell value, whatever the active category", {
  skip_if_not_installed("fasterize")
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6, ymin = 0, ymax = 6, crs = "EPSG:3978")
  ids <- rep(1:3, length.out = terra::ncell(r))
  lcc <- terra::rast(r, vals = rep(c(210L, 220L), length.out = terra::ncell(r)))
  rtm <- terra::rast(r, vals = 1L)
  ecoregionTable <- data.table::data.table(
    ID = factor(c("1", "2", "3")),
    ecoregionName = factor(c("138", "139", "140"))
  )

  ## as prepEcoregions() builds it: ecoregionName is the only label column
  base <- terra::rast(r, vals = ids)
  levels(base) <- data.frame(ID = 1:3, ecoregionName = c("138", "139", "140"))

  ## another text column ahead of ecoregionName, then a file round trip
  multi <- terra::rast(r, vals = ids)
  levels(multi) <- data.frame(ID = 1:3, code = c("a", "b", "c"), ecoregionName = c("138", "139", "140"))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(multi, f)
  multi <- terra::rast(f)
  ## premise: ecoregionName is not the active category
  expect_false(identical(names(terra::cats(multi)[[1]])[terra::activeCat(multi) + 1L], "ecoregionName"))

  expected <- ecoregionProducer(list(base, lcc), rasterToMatch = rtm, ecoregionTable = ecoregionTable)
  got <- ecoregionProducer(list(multi, lcc), rasterToMatch = rtm, ecoregionTable = ecoregionTable)

  expect_identical(terra::values(got$ecoregionMap, mat = FALSE), terra::values(expected$ecoregionMap, mat = FALSE))
  expect_identical(got$ecoregion, expected$ecoregion)
  expect_identical(as.character(got$ecoregion$ecoregion_lcc), c("1_210", "1_220", "2_210", "2_220", "3_210", "3_220"))
})

## prepEcoregions() must normalize a supplied categorical ecoregion raster the same way for terra
## and raster: the raster's category table and ecoregionTable both (ID, ecoregionName), padded IDs.
## The SpatRaster branch used to rename only the table (and not pad IDs), and the RasterLayer branch
## renamed nothing, so a label column not literally named ecoregionName failed the join in
## ecoregionProducer(). 12 ecoregions, so ID padding ("01" vs "1") matters.
test_that("prepEcoregions() gives the same result for terra and raster categorical inputs", {
  withr::local_options(list(reproducible.useCache = FALSE, reproducible.cachePath = withr::local_tempdir()))
  nEco <- 12L
  r <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 12, crs = "EPSG:3978")
  ids <- rep(seq_len(nEco), length.out = terra::ncell(r))
  names_ <- as.character(130 + seq_len(nEco))
  lcc <- terra::rast(r, vals = rep(c(210L, 220L, 230L), length.out = terra::ncell(r)))
  rtm <- terra::rast(r, vals = 1L)

  mk <- function(df) {
    x <- terra::rast(r, vals = ids)
    levels(x) <- df
    x
  }
  ## as prepEcoregions() builds it from polygons: ecoregionName is the only label column
  base <- mk(data.frame(ID = seq_len(nEco), ecoregionName = names_))
  ## label column named something else (SpatRaster)
  spat <- mk(data.frame(ID = seq_len(nEco), ECOREGION = names_))
  ## several label columns; the active one holds the ecoregion labels
  multi <- mk(data.frame(ID = seq_len(nEco), code = letters[seq_len(nEco)], ECOREGION = names_))
  terra::activeCat(multi) <- "ECOREGION"
  ## the same raster as a RasterLayer (raster has no active category: first label column)
  rl <- raster::raster(spat)

  run <- function(eco, lccIn = lcc, rtmIn = rtm) {
    out <- prepEcoregions(ecoregionRst = eco, ecoregionLayer = NULL, rasterToMatchLarge = rtmIn,
                          rstLCCAdj = lccIn, pixelsToRm = NULL, cacheTags = "test")
    m <- out$ecoregionMap
    if (inherits(m, "Raster")) m <- terra::rast(m) ## rast() on a SpatRaster would drop its values
    list(map = terra::values(m, mat = FALSE), table = out$ecoregion)
  }
  expected <- run(base)
  expect_false(anyNA(expected$map)) ## the comparison below is of real values
  expect_identical(run(spat), expected)
  expect_identical(run(multi), expected)
  expect_identical(run(rl, raster::raster(lcc), raster::raster(rtm)), expected)
  ## the padded IDs joined: every ecoregion present, none dropped by the final na.omit()
  expect_setequal(as.character(expected$table$ecoregionName), names_)
})
