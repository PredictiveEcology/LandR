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
#' LandR applies this remap itself, only for the duration of each SCANFI download and only if
#' `getOption("reproducible.urlRemap")` is `NULL` and `getOption("LandR.scanfiMirror", TRUE)` is
#' `TRUE`. The global option is never changed. A remap of your own always wins; if you use one,
#' include these rows in it (the manifest is
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
  manifest <- .scanfiManifest()
  if (is.null(manifest)) {
    return(NULL)
  }
  reproducible::makeUrlRemap(manifest)
}

## The manifest of SCANFI v2 files on the mirror (NULL if it is not installed).
.scanfiManifest <- function() {
  f <- system.file("extdata", "arbutus_manifest_SCANFI_v2.csv", package = "LandR")
  if (!nzchar(f)) {
    return(NULL)
  }
  utils::read.csv(f, stringsAsFactors = FALSE)
}

## Evaluate `expr` with the effective SCANFI remap as `reproducible.urlRemap`, for the duration of
## the call only: LandR never sets the option globally, and a remap of the user's always wins.
.withScanfiRemap <- function(expr) {
  if (is.null(getOption("reproducible.urlRemap")) && isTRUE(getOption("LandR.scanfiMirror", TRUE))) {
    remap <- scanfiUrlRemap()
    if (!is.null(remap)) {
      old <- options(reproducible.urlRemap = remap)
      on.exit(options(old), add = TRUE)
    }
  }
  expr
}

## Will LandR's own remap send these SCANFI V2 `url`s (Drive ids or folders) to the mirror?
.remapsToMirror <- function(url, dataVersion) {
  if (!identical(dataVersion, "V2") || !is.null(getOption("reproducible.urlRemap")) ||
      !isTRUE(getOption("LandR.scanfiMirror", TRUE)) || !length(url)) {
    return(FALSE)
  }
  known <- c(.scanfiManifest()$id, names(attr(scanfiUrlRemap(), "byDir")))
  all(sub(".*/", "", as.character(url)) %in% known |
        startsWith(as.character(url), "https://object-arbutus."))
}
