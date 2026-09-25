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
