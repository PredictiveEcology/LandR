#' Obtain a SCANFI structural attribute layer
#'
#' SCANFI publishes canopy height and canopy closure alongside the biomass layer that
#' [prepRawBiomassMap()] fetches. They describe how much structure a pixel carries, independent
#' of what species carry it, which is what makes them useful as controls: two pixels of the same
#' height and closure hold comparable stands, so a difference in biomass between them is a
#' difference in composition rather than in site quality.
#'
#' @param attribute `"height"` (metres) or `"closure"` (percent crown closure).
#' @param year one of `r paste(LandR:::.scanfi_v2_years, collapse = ", ")`.
#' @param dataVersion SCANFI product version. Only `"V2"` publishes these layers.
#' @param ... passed to [reproducible::prepInputs()] and [reproducible::Cache()]; pass one of
#'   `to`, `cropTo` or `rasterToMatch` -- the national rasters are large.
#'
#' @return a `SpatRaster`
#'
#' @export
prepInputs_SCANFI_structure <- function(attribute = c("height", "closure"), year = 2020,
                                        dataVersion = "V2", ...) {
  attribute <- match.arg(attribute)
  if (!identical(dataVersion, "V2")) {
    stop("SCANFI height and closure layers are published for V2 only")
  }
  if (!(as.integer(year) %in% .scanfi_v2_years)) {
    stop("SCANFI V2 ", attribute, " is available for ",
         paste(.scanfi_v2_years, collapse = ", "), " only")
  }
  Args <- list(...)
  if (is.null(Args$url)) {
    Args$url <- .scanfiV2Source(attribute, year)$url
  }
  ## as in prepRawBiomassMap(): a failed prepInputs often leaves a corrupt partial download, and
  ## a successful one is held by the Cache below rather than re-fetched.
  if (is.null(Args$overwrite)) Args$overwrite <- TRUE
  if (is.null(Args$method)) Args$method <- "bilinear"

  do.call(prepInputs, args = Args) |>
    Cache(quick = "writeTo",
          .functionName = paste0("prepInputs_SCANFI_", attribute),
          omitArgs = c("destinationPath", "targetFile", "stable"))
}
