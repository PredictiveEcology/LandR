## ---------------------------------------------------------------------------
## SCANFI v3 land cover
##
## Unlike V1/V2, which are distributed pre-converted to Canada LCC codes through
## Google Drive, V3 is published directly over https as one Cloud-Optimized GeoTIFF
## per year, with its own 20-class legend. This file reads a study-area window of
## that COG via GDAL's /vsicurl driver (so the ~3.4 GB national file is never
## downloaded whole) and applies the crosswalk to Canada LCC codes on the fly.
## ---------------------------------------------------------------------------

## SCANFI v3 publishes one landcover layer per year, 1985-2025
.scanfi_v3_years <- 1985:2025

## The v3 COG server (an S3-backed https endpoint) returns 403 to requests without a
## browser-like User-Agent; this is set on the GDAL config before every read.
.scanfiV3UserAgent <- paste0(
  "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 ",
  "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
)

#' Crosswalk from SCANFI v3 land-cover codes to Canada LCC codes
#'
#' SCANFI v3's landcover legend (see the Open Canada record
#' <https://open.canada.ca/data/en/dataset/50f132f9-f312-4951-bb9f-9ea99580f29f>) has 20
#' classes and no wetland class. This table maps each to the Canada LCC code `LandR` and
#' `fireSenseUtils` use elsewhere (the same codes NTEMS and SCANFI V2 use). Code 4, burn
#' scars, has no Canada LCC equivalent -- it is not shrubland, and collapsing it into 50
#' would blur a recent burn with land that has always been shrubby -- so it keeps its own
#' code, 60, which [.applyForestLand()] still relabels to `disturbedCode` (240) wherever the
#' pixel is forest land. Cropland, urban and road (17-19) have no analogue in the
#' 20-230 Canada LCC scheme used here and are mapped to 0 (no data), matching how V2 and
#' NTEMS handle land uses outside that scheme. Code 255 (the source's NoData value) maps
#' to `NA`.
#'
#' @format A `data.frame` with one row per SCANFI v3 code (1-20, plus 255 for NoData) and
#'   columns `scanfiV3` (the source code), `lcc` (the Canada LCC code, or `NA`), and
#'   `description`.
#'
#' @export
scanfiV3ToCanadaLCC <- data.frame(
  scanfiV3 = c(1:20, 255L),
  lcc = c(
    20, 30, 33, 60, 40, 100, 50, 50, 220, 230,
    210, 210, 210, 210, 210, 210, 0, 0, 0, 31,
    NA
  ),
  description = c(
    "Water", "Rock", "Soil", "Burn scars (SCANFI v3 only; not a Canada LCC code)",
    "Lichen", "Herbaceous", "Low shrubs", "Tall shrubs", "Treed broadleaf",
    "Treed mixed", "Treed coniferous", "Treed coniferous with lichen",
    "Treed coniferous with rock/soil", "Treed coniferous with herbs",
    "Treed coniferous with low shrub", "Treed coniferous with tall shrub",
    "Cropland", "Urban", "Road", "Snow/Ice",
    "No data"
  ),
  stringsAsFactors = FALSE
)

#' Where one year of SCANFI v3 land cover comes from
#'
#' @param year A year SCANFI v3 publishes, `r min(.scanfi_v3_years)`-`r max(.scanfi_v3_years)`.
#'
#' @return The https URL of that year's national Cloud-Optimized GeoTIFF.
#'
#' @keywords internal
.scanfiV3Url <- function(year) {
  if (!(year %in% .scanfi_v3_years)) {
    stop("SCANFI V3 Landcover does not exist for this year")
  }
  paste0(
    "https://download-telecharger.services.geo.ca/pub/nrcan_rncan/Forests_Foret/",
    "SCANFI/v3/cog_SCANFI_landcover_", year, "_v3_20260528.tif"
  )
}

#' Recode a SCANFI v3 land-cover raster to Canada LCC codes
#'
#' @param lcc A `SpatRaster` of SCANFI v3 codes (1-20, `NA` for NoData).
#'
#' @return A `SpatRaster` of Canada LCC codes, per [scanfiV3ToCanadaLCC].
#'
#' @keywords internal
.applySCANFIv3Crosswalk <- function(lcc) {
  terra::subst(
    lcc,
    from = scanfiV3ToCanadaLCC$scanfiV3,
    to = scanfiV3ToCanadaLCC$lcc,
    datatype = "INT1U",
    NAflag = 255
  )
}

#' Read one year of SCANFI v3 land cover, windowed to a study area
#'
#' Reads the national COG through GDAL's `/vsicurl` driver, so only the blocks
#' overlapping the study area are fetched -- the ~3.4 GB national file is never downloaded whole.
#' The study area is given the same way as to [reproducible::postProcessTo()] (and so to
#' `prepInputs()` for SCANFI V1/V2): `to`, `cropTo`, `maskTo`, `projectTo`. At least one is
#' required. The result is always a new raster (in memory or a temporary file), never one that
#' still points at the remote file.
#'
#' @param year A year SCANFI v3 publishes.
#' @param to,cropTo,maskTo,projectTo The study area, as in [reproducible::postProcessTo()].
#' @param method Resampling method, as in [reproducible::postProcessTo()].
#' @param what Passed to `.withSCANFIv3Access()`, for its error message.
#' @param url The file to read; defaults to that year's NRCan COG. An `http(s)` URL is read
#'   through `/vsicurl`; anything else is read as a local file (used by the tests).
#'
#' @return A `SpatRaster` of raw SCANFI v3 codes for the study area.
#'
#' @keywords internal
.readSCANFIv3 <- function(year, to = NULL, cropTo = NULL, maskTo = NULL, projectTo = NULL,
                          method = "near", what = "the SCANFI v3 land cover map",
                          url = .scanfiV3Url(year)) {
  studyArea <- list(to = to, cropTo = cropTo, maskTo = maskTo, projectTo = projectTo)
  studyArea <- studyArea[!vapply(studyArea, is.null, logical(1))]
  if (!length(studyArea)) {
    stop("SCANFI V3 is read as a study-area window: pass `to`, `cropTo`, `maskTo` or `projectTo` ",
         "(or `rasterToMatch` to prepInputs_SCANFI_LCC_FAO()).")
  }
  src <- if (grepl("^https?://", url)) paste0("/vsicurl/", url) else url

  oldUA <- Sys.getenv("GDAL_HTTP_USERAGENT", unset = NA)
  Sys.setenv(GDAL_HTTP_USERAGENT = .scanfiV3UserAgent)
  on.exit({
    if (is.na(oldUA)) {
      Sys.unsetenv("GDAL_HTTP_USERAGENT")
    } else {
      Sys.setenv(GDAL_HTTP_USERAGENT = oldUA)
    }
  }, add = TRUE)

  .withSCANFIv3Access({
    r <- terra::rast(src)
    r <- do.call(reproducible::postProcessTo, c(list(from = r, method = method), studyArea))
    ## A window that needed no change can come back still pointing at the source file; Cache()
    ## would then store (and mangle) the remote path instead of the data.
    if (any(terra::sources(r) %in% src)) {
      r <- terra::writeRaster(r, tempfile(fileext = ".tif"), datatype = "INT1U", NAflag = 255)
    }
    r
  }, what = what, dataYear = year)
}
