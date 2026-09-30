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
