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
  f <- system.file("extdata", "arbutus_manifest_SCANFI_v2.csv", package = "LandR")
  if (!nzchar(f) || !"makeUrlRemap" %in% getNamespaceExports("reproducible")) {
    return(NULL)
  }
  reproducible::makeUrlRemap(utils::read.csv(f, stringsAsFactors = FALSE))
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
