## a 1-row raster of land-cover codes
lccRow <- function(v) {
  r <- terra::rast(
    nrows = 1,
    ncols = length(v),
    xmin = 0,
    xmax = length(v),
    ymin = 0,
    ymax = 1,
    crs = "EPSG:3978"
  )
  terra::values(r) <- v
  r
}

## the simulation land cover (LCC) is kept, except where the current-condition one (CC) is water
## (18) or barren (16), which make the pixel non-vegetated, or urban (17), which is re-imputed
remapRow <- function(lccClasses, ccClasses) {
  remap <- data.table::as.data.table(expand.grid(LCC = c(NA, lccClasses), CC = c(NA, ccClasses)))
  remap[, newLCC := LCC]
  remap[CC %in% c(16L, 18L), newLCC := NA_integer_]
  remap[CC %in% 17L, newLCC := 99L]
  remap[]
}

test_that("overlayLCCs() applies the remap table's non-vegetated pixels", {
  withr::local_options(reproducible.rasterRead = "terra::rast")
  lcc <- lccRow(c(210L, 220L, 210L, 30L))
  cc <- lccRow(c(1L, 18L, 16L, 1L))
  out <- overlayLCCs(
    list(CC = cc, LCC = lcc),
    forestedList = list(CC = 1L, LCC = c(210L, 220L)),
    outputLayer = "LCC",
    remapTable = remapRow(c(30L, 210L, 220L), c(1L, 16L, 18L)),
    classesToReplace = 99L,
    availableERC_by_Sp = NULL
  )
  expect_equal(terra::values(out, mat = FALSE), c(210L, NA, NA, 30L))
})

test_that("overlayLCCs() keeps the non-vegetated pixels beside a re-imputed one", {
  withr::local_options(reproducible.rasterRead = "terra::rast")
  ## the urban pixel's nearest neighbours are 210s, which the remap makes non-vegetated, so the
  ## only class it can take is the 220 three pixels away
  lcc <- lccRow(c(210L, 230L, 210L, 210L, 220L))
  cc <- lccRow(c(18L, 17L, 18L, 18L, 1L))
  out <- overlayLCCs(
    list(CC = cc, LCC = lcc),
    forestedList = list(CC = 1L, LCC = c(210L, 220L, 230L)),
    outputLayer = "LCC",
    remapTable = remapRow(c(210L, 220L, 230L), c(1L, 17L, 18L)),
    classesToReplace = 99L,
    availableERC_by_Sp = NULL
  )
  expect_equal(terra::values(out, mat = FALSE), c(NA, 220L, NA, NA, 220L))
})
