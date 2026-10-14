#' Resolve SCANFI Google Drive ids to the PredictiveEcology mirror
#'
#' LandR addresses SCANFI v2 files by their Google Drive ids, and some of those ids no longer
#' resolve for an anonymous user -- the 2020 land cover and 2020 stand age return 404, which
#' stops the default SCANFI path of `Biomass_borealDataPrep`. The same files are mirrored on
#' the PredictiveEcology object store on arbutus (Digital Research Alliance of Canada). This
#' builds a [reproducible::makeUrlRemap()] hook from the manifest shipped with LandR, so a
#' Drive id -- or a SCANFI species *folder* -- is fetched from the mirror instead, with no
#' Google login and no change at any call site.
#'
#' LandR sets `options(reproducible.urlRemap = scanfiUrlRemap())` when it is loaded, but only
#' if no remap is set already and `getOption("LandR.scanfiMirror", TRUE)` is `TRUE`. A remap
#' of your own is never replaced; if you use one, include these rows in it (the manifest is
#' `system.file("extdata", "arbutus_manifest_SCANFI_v2.csv", package = "LandR")`).
#'
#' @return A remap function for `reproducible.urlRemap`, or `NULL` if this `reproducible` has
#'   no `makeUrlRemap()`.
#'
#' @export
scanfiUrlRemap <- function() {
  if (!"makeUrlRemap" %in% getNamespaceExports("reproducible")) {
    return(NULL)
  }
  manifest <- tryCatch(.scanfiManifest(), error = function(e) NULL)
  if (is.null(manifest)) {
    return(NULL)
  }
  reproducible::makeUrlRemap(manifest)
}

## Every url of the manifest starts with this.
.arbutusPrefix <- "https://object-arbutus.alliancecan.ca/"
.isArbutusUrl <- function(url) {
  length(url) > 0 && all(startsWith(as.character(url), .arbutusPrefix))
}

## The one table of SCANFI v2 files on the arbutus mirror: filename, key, url, Drive id, type.
.scanfiManifest <- function() {
  f <- system.file("extdata", "arbutus_manifest_SCANFI_v2.csv", package = "LandR")
  if (!nzchar(f)) {
    stop("the SCANFI v2 manifest is not installed with LandR")
  }
  utils::read.csv(f, stringsAsFactors = FALSE)
}

## File-name stem of each national SCANFI v2 layer: SCANFI_<stem>_<year>_v2_<date>.tif
.scanfiV2Layers <- c(
  age = "age_median",
  biomass = "att_biomass",
  height = "att_height",
  closure = "att_closure",
  lccRaw = "att_nfiLandcover", ## the 1-8 codes that convert_SCANFI_LCC_codes() recodes
  lcc = "att_nfiLandcover_CanadaLCCclassCodes"
)

## Rows of the manifest that are files of one SCANFI v2 `year`, with names matching `pattern`.
.scanfiV2Rows <- function(year, pattern) {
  year <- as.character(as.integer(year))
  if (!(year %in% .scanfi_v2_years)) {
    stop("SCANFI V2 data is currently available for ",
         paste(.scanfi_v2_years, collapse = ", "), " only")
  }
  m <- .scanfiManifest()
  m <- m[m$type == "file" & startsWith(m$key, paste0("SCANFI_v2/", year, "/")) &
           grepl(pattern, m$filename), , drop = FALSE]
  if (!nrow(m)) {
    stop("no SCANFI V2 file matching '", pattern, "' for ", year, " in the arbutus manifest")
  }
  m
}

#' Where one national SCANFI v2 layer of one year comes from
#'
#' @param layer one of `names(.scanfiV2Layers)`.
#' @param year a year in `.scanfi_v2_years`.
#'
#' @return A list with the arbutus `url` and the `targetFile` name.
#'
#' @keywords internal
.scanfiV2Source <- function(layer, year) {
  layer <- match.arg(layer, names(.scanfiV2Layers))
  m <- .scanfiV2Rows(year, paste0("^SCANFI_", .scanfiV2Layers[[layer]], "_[0-9]{4}_v2_[0-9]+[.]tif$"))
  list(url = m$url, targetFile = m$filename)
}

## The species-layer files (name, url) of one SCANFI v2 year on the mirror.
.scanfiV2SpeciesFiles <- function(year) {
  m <- .scanfiV2Rows(year, "^SCANFI_sps")
  data.table::data.table(name = m$filename, url = m$url)
}

## Called from .onLoad(): install the SCANFI mirror unless the user opted out or has a remap.
.setScanfiMirror <- function() {
  if (!isTRUE(getOption("LandR.scanfiMirror", TRUE))) {
    return(invisible(FALSE))
  }
  if (!is.null(getOption("reproducible.urlRemap"))) {
    return(invisible(FALSE))
  }
  remap <- tryCatch(scanfiUrlRemap(), error = function(e) NULL)
  if (is.null(remap)) {
    return(invisible(FALSE))
  }
  options(reproducible.urlRemap = remap)
  invisible(TRUE)
}
