## ---------------------------------------------------------------------------
## Precompute the artifacts backing vignettes/scanfi-summary.Rmd
##
## Run this manually on a workstation; it is NOT part of the build and must
## never run in CI. It reads a national forest raster end to end (tens of GB for
## SCANFI at 30 m) and downloads the ecostratification polygons, neither of
## which a CI runner can do.
##
## What it writes (all small, all committed):
##
##   inst/extdata/<prefix>_counts_ecozone.csv.gz  exact value-frequency table,
##                                                computed from the FULL raster
##   inst/extdata/<prefix>_display.tif            coarsened raster, for maps only
##   inst/extdata/ecozones.gpkg                   simplified ecozone polygons
##   inst/extdata/vignette-artifacts.dcf          hashes + provenance
##
## The vignette recomputes every summary statistic, histogram and map from these
## at build time, so nothing it reports can drift from the code in R/. The one
## thing that *can* drift is the counts table itself; `.write_artifact_manifest()`
## records a hash of the code that produced it, and
## tests/testthat/test-vignette-artifacts.R fails when they diverge.
##
## Re-run this when:
##   - tests/testthat/test-vignette-artifacts.R reports a stale counts hash, or
##   - the upstream data product is revised (new SCANFI version or year).
## ---------------------------------------------------------------------------

## ---- configuration --------------------------------------------------------

DATA_SOURCE <- Sys.getenv("LANDR_VIGNETTE_SOURCE", "KNN")
DATA_YEAR <- as.integer(Sys.getenv("LANDR_VIGNETTE_YEAR", "2011"))
DATA_VERSION <- Sys.getenv("LANDR_VIGNETTE_VERSION", "V2")

## Target resolution of the committed display raster. The statistics come from
## the full-resolution raster; this exists only so the vignette can draw maps
## without shipping a multi-GB file.
DISPLAY_RES_M <- 5000

## Simplification tolerance for the committed polygons (metres). Well below what
## is visible at the display resolution.
POLYGON_TOLERANCE_M <- 2000

ID_COL <- "ZONE_NAME"
ECOZONE_URL <- "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"

out_dir <- file.path("inst", "extdata")
prefix <- tolower(paste0(DATA_SOURCE, "_age_", DATA_YEAR))

## ---- setup ----------------------------------------------------------------

stopifnot(
  "run from the package root" = dir.exists("R") && file.exists("DESCRIPTION")
)

library(sf)
pkgload::load_all(".", quiet = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Process on disk with bounded RAM (see ?calc_raster_counts). `tempdir` MUST
## point at fast *local* storage -- never NFS, and never a tmpfs such as /tmp on
## systems where that is RAM-backed.
scratch <- Sys.getenv("LANDR_VIGNETTE_SCRATCH", tempfile("landr_precompute_"))
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(memmax = 16, todisk = TRUE, tempdir = scratch)

message("scratch: ", scratch)

## ---- inputs ---------------------------------------------------------------

message("fetching ", DATA_SOURCE, " stand age (", DATA_YEAR, ") ...")
age <- prepInputsStandAgeMap(
  dataSource = DATA_SOURCE,
  dataYear = DATA_YEAR,
  dataVersion = DATA_VERSION,
  destinationPath = scratch
)

message("fetching ecozones ...")
ecozones <- reproducible::prepInputs(
  url = ECOZONE_URL,
  destinationPath = scratch
) |>
  sf::st_make_valid()

## ---- 1. exact counts, from the FULL-resolution raster ---------------------
##
## This is the expensive step, and the only one that touches every pixel.

message("computing value-frequency table (this is the slow part) ...")
counts <- calc_raster_counts(raster = age, polygons = ecozones, polygon_id = ID_COL)

counts_file <- file.path(out_dir, paste0(prefix, "_counts_ecozone.csv.gz"))
utils::write.csv(counts, gzfile(counts_file), row.names = FALSE)
message("wrote ", counts_file, " (", nrow(counts), " rows, ",
        round(file.size(counts_file) / 1024), " KB)")

gc()

## ---- 2. coarsened display raster ------------------------------------------

fact <- max(1L, as.integer(round(DISPLAY_RES_M / terra::res(age)[1])))
message("aggregating for display (factor ", fact, ") ...")

display <- terra::aggregate(age, fact = fact, fun = "mean", na.rm = TRUE)

display_file <- file.path(out_dir, paste0(prefix, "_display.tif"))
terra::writeRaster(
  display,
  display_file,
  overwrite = TRUE,
  datatype = "INT2U",
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9", "TILED=YES")
)
message("wrote ", display_file, " (",
        paste(dim(display)[1:2], collapse = "x"), " cells, ",
        round(file.size(display_file) / 1024), " KB)")

gc()

## ---- 3. simplified polygons -----------------------------------------------

polygons_file <- file.path(out_dir, "ecozones.gpkg")

ecozones |>
  subset(select = ID_COL) |>
  sf::st_transform(sf::st_crs(age)) |>
  sf::st_simplify(dTolerance = POLYGON_TOLERANCE_M) |>
  sf::st_make_valid() |>
  sf::st_write(polygons_file, delete_dsn = TRUE, quiet = TRUE)

message("wrote ", polygons_file, " (", round(file.size(polygons_file) / 1024), " KB)")

## ---- 4. manifest ----------------------------------------------------------

.write_artifact_manifest(
  file.path(out_dir, "vignette-artifacts.dcf"),
  provenance = list(
    DataSource = DATA_SOURCE,
    DataYear = DATA_YEAR,
    DataVersion = if (identical(DATA_SOURCE, "SCANFI")) DATA_VERSION else NA,
    PolygonSource = ECOZONE_URL,
    CountsFile = basename(counts_file),
    DisplayFile = basename(display_file),
    DisplayResolution = paste0(round(terra::res(display)[1]), " m"),
    PolygonFile = basename(polygons_file)
  )
)

message("wrote ", file.path(out_dir, "vignette-artifacts.dcf"))
message("\nDone. Commit inst/extdata/, then rebuild the vignette to check it.")
