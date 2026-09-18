## SCANFI V2 Google Drive ids for the two structural attribute layers. Same provenance as the
## biomass ids in `prepRawBiomassMap()`; all nine V2 years exist for both.
.scanfiStructureIds <- list(
  height = c(
    "1985" = "1KFkOklhVOzFaaS-t1MtWlg-MxMWNS41g",
    "1990" = "1sTU6ODgKs6ZVAbYCmb05qm0BTPhZ0260",
    "1995" = "1G5L96tnxFOliJIz5BC0953xZedqiS7sN",
    "2000" = "1KS8fLVywjBorvL8kICo7pUsgG8G5q2Oq",
    "2005" = "1UtUaT_E6x7z4sZAApjKQpZOSKh_ruGH6",
    "2010" = "1E27EvIDLXHWz14RCd0t6NlV91ob9NXpu",
    "2015" = "11htblc85pwhz6z65uaVy7IxwRO_EwQsv",
    "2020" = "1krOB9PIV1OnJS2N8MzfrR8RIpljkZ1uN",
    "2025" = "1rzdbJppv_exhGaHqEqNXPoUkyIhpKlPY"
  ),
  closure = c(
    "1985" = "1bACZ3SxvXp4nFPyiTUzJgrRi7ngXhTPr",
    "1990" = "1xvIR8iL5pTHgfnzbzhuw4BmVC7T93F34",
    "1995" = "1x9yRwiPFJgN01I8y8OkduqEFHGaN8iLB",
    "2000" = "18DzUWk_QsSTfkuJgIfisdT5r1hbG6ilE",
    "2005" = "1ML-iCtN2um70RxYVN1PYlktPGmcDdep2",
    "2010" = "1qdP_5fA_49t6cNjewnjDq71C7xsEo2AN",
    "2015" = "18owQ9Ix6h-q-rQjiplp6Ls7v5DMHO-zv",
    "2020" = "1yrGovqTd-qVeawaxVmWv2E0q43GfDLkh",
    "2025" = "1jDsWKdB5I6e3ELrssKvzc7IqnRJo6-jY"
  )
)

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
    Args$url <- paste0("https://drive.google.com/file/d/",
                       .scanfiStructureIds[[attribute]][[as.character(as.integer(year))]])
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
