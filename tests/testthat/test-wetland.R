## SCANFI land cover has no wetland classes; CWIM3A supplies the site, and wetlandToLCC() turns
## it into the NTEMS codes 80/81. No network here: the CWIM source is a small synthetic file
## with the same codes and the same NoData value (15).

test_that("wetlandToLCC recodes wet pixels by whether they are treed", {
  lcc <- c(210, 220, 230, 240, 50, 100, 40, 20, 80, 81, 210, NA)
  wet <- c(1,   1,   1,   1,   1,  1,   1,  1,  1,  1,  0,   1)
  expect_identical(
    wetlandToLCC(lcc, wet),
    c(81, 81, 81, 81, 80, 80, 80, 20, 80, 81, 210, NA)
  )
})

test_that("wetlandToLCC leaves dry and unknown-site pixels alone", {
  lcc <- c(210, 50, 240)
  expect_identical(wetlandToLCC(lcc, c(0, 0, 0)), lcc)
  expect_identical(wetlandToLCC(lcc, c(NA, NA, NA)), lcc)
})

test_that("240 counts as treed: disturbed forest land is still forest ground", {
  expect_identical(wetlandToLCC(240, 1), 81)
  expect_identical(wetlandToLCC(240, 1, treedClasses = c(210, 220, 230)), 80)
})

test_that("wetlandToLCC works on rasters and keeps the geometry", {
  withr::local_package("terra")
  lcc <- rast(nrows = 2, ncols = 2, vals = c(210, 50, 20, 230))
  wet <- rast(lcc, vals = c(1, 1, 1, 0))
  out <- wetlandToLCC(lcc, wet)
  expect_true(compareGeom(out, lcc))
  expect_identical(as.vector(values(out)), c(81, 80, 20, 230))
})

test_that("wetlandToLCC refuses mismatched inputs", {
  expect_error(wetlandToLCC(c(210, 220), 1:3), "same number of cells")
})

## A synthetic CWIM3A: EPSG:3979 at 10 m, 4 x 2 target cells of 240 m (24 x 24 source cells
## each). Target columns, left to right: all bog; all NoData; all shallow water; and a split
## cell whose wet fraction is set exactly.
## Called from inside tests that have already attached terra; attaching it again here would
## make the two local_package() calls detach it twice.
makeCwim <- function(wetFracLastCol) {
  n <- 24
  src <- rast(nrows = 2 * n, ncols = 4 * n, xmin = 0, xmax = 4 * 240, ymin = 0, ymax = 2 * 240,
              crs = "EPSG:3979")
  m <- matrix(NA_integer_, nrow = 2 * n, ncol = 4 * n)
  m[, 1:n] <- 1L                                 ## bog
  m[, (2 * n + 1):(3 * n)] <- 5L                 ## shallow water
  k <- round(wetFracLastCol * n)
  last <- (3 * n + 1):(4 * n)
  m[, last] <- 15L                               ## NoData ...
  if (k > 0) m[, last[seq_len(k)]] <- 4L         ## ... except k columns of swamp
  m[is.na(m)] <- 15L
  values(src) <- as.vector(t(m))
  f <- tempfile(fileext = ".tif")
  writeRaster(src, f, datatype = "INT1U", NAflag = 15, overwrite = TRUE)
  f
}

test_that("prepInputs_CWIM: bog is wet, NoData and shallow water are not", {
  withr::local_package("terra")
  f <- makeCwim(wetFracLastCol = 0.75)
  to <- rast(nrows = 2, ncols = 4, xmin = 0, xmax = 960, ymin = 0, ymax = 480, crs = "EPSG:3979")
  values(to) <- 1
  out <- prepInputs_CWIM(to, url = f)
  expect_true(compareGeom(out, to))
  expect_identical(names(out), "wetland")
  expect_identical(as.vector(values(out)), rep(c(1, 0, 0, 1), 2))
})

test_that("prepInputs_CWIM: the wet fraction decides a mixed cell", {
  withr::local_package("terra")
  to <- rast(nrows = 2, ncols = 4, xmin = 0, xmax = 960, ymin = 0, ymax = 480, crs = "EPSG:3979")
  values(to) <- 1
  lastCol <- function(frac, thr) {
    v <- as.vector(values(prepInputs_CWIM(to, url = makeCwim(frac), wetThreshold = thr)))
    v[c(4, 8)]
  }
  expect_identical(lastCol(0.25, 0.5), c(0, 0))
  expect_identical(lastCol(0.50, 0.5), c(1, 1))   ## the threshold is inclusive
  expect_identical(lastCol(0.25, 0.2), c(1, 1))
})

test_that("prepInputs_CWIM: shallow water counts only if asked", {
  withr::local_package("terra")
  to <- rast(nrows = 2, ncols = 4, xmin = 0, xmax = 960, ymin = 0, ymax = 480, crs = "EPSG:3979")
  values(to) <- 1
  out <- prepInputs_CWIM(to, url = makeCwim(0), wetClasses = 1:5)
  expect_identical(as.vector(values(out))[c(3, 7)], c(1, 1))
})

test_that("prepInputs_CWIM keeps `to`'s NA cells NA and projects across CRSs", {
  withr::local_package("terra")
  f <- makeCwim(1)
  ## same Lambert family as SCANFI (lat_0 = 0) rather than EPSG:3979 (lat_0 = 49)
  crsS <- "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  tmpl <- rast(nrows = 2, ncols = 4, xmin = 0, xmax = 960, ymin = 0, ymax = 480,
               crs = "EPSG:3979", vals = 1)
  toS <- project(tmpl, crsS, res = 240)
  toS[1] <- NA
  out <- prepInputs_CWIM(toS, url = f)
  expect_true(compareGeom(out, toS))
  expect_true(is.na(values(out)[1]))
  expect_true(all(values(out)[-1] %in% c(0, 1, NA)))
})

test_that("prepInputs_CWIM wants a raster grid", {
  expect_error(prepInputs_CWIM(to = "not a raster", url = "x.tif"), "must be a SpatRaster")
})
