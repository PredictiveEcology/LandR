## Public cloud-optimised GeoTIFF of the Canadian Wetland Inventory Map v3A (Open Government
## Licence - Canada). Listed at open.canada.ca, dataset 87127901-bd6d-46de-9142-e1362d980174.
.cwimURL <- paste0(
  "https://datacube-prod-data-public.s3.amazonaws.com/store/land/wetlands/",
  "wetland-inventory/canadian-wetland-inventory-v3a-classification.tif"
)

#' Wetland (site) layer from the Canadian Wetland Inventory Map
#'
#' SCANFI's land cover has no wetland classes, so a land-cover map built from it cannot say
#' which forest stands on wet ground -- the distinction NTEMS draws with classes 80
#' (wetland) and 81 (treed wetland). This builds that site layer from the Canadian Wetland
#' Inventory Map v3A (CWIM3A; Mahdianpari et al. 2021): 10 m, national, 2016-2020 imagery.
#'
#' Only the study window is read: the source is a cloud-optimised GeoTIFF, so `terra` fetches
#' the tiles it needs rather than the national raster.
#'
#' CWIM3A codes are 1 Bog, 2 Fen, 3 Marsh, 4 Swamp, 5 Shallow water; 15 is NoData, and CWIM3A
#' maps nothing but those five classes, so NoData means "not wetland". Shallow water is not
#' counted as wet by default: it is open water, and it coincides with land-cover water (class
#' 20) almost exactly -- on a 60 km boreal test window, 6,593 of 6,764 CWIM shallow-water cells
#' were VLCE2 water.
#'
#' The 10 m cells are summarised onto `to` as the *fraction* of each target cell that is wet;
#' a cell is wet when that fraction is at least `wetThreshold`. For a 0/1 layer that is the
#' majority rule, without `"mode"`'s arbitrary tie-breaking.
#'
#' @param to A `SpatRaster` giving the grid to return, and the study window. Its `NA` cells
#'   are `NA` in the result.
#' @param url The CWIM3A source. A URL is read through GDAL's `/vsicurl/`; a local path is
#'   read directly (useful for a pre-cut regional copy).
#' @param wetClasses CWIM3A codes counted as wet. Default: Bog, Fen, Marsh, Swamp.
#' @param wetThreshold Minimum wet fraction of a target cell for it to be wet.
#' @param writeTo Optional file to write the result to.
#'
#' @return A single-layer `SpatRaster` on `to`'s grid, named `"wetland"`: `1` wet, `0` not wet.
#'
#' @references Mahdianpari, M., et al. (2021). The Third Generation of Pan-Canadian Wetland
#'   Map at 10 m Resolution Using Multisource Earth Observation Data on Cloud Computing
#'   Platform. IEEE JSTARS 14, 8789-8803.
#'
#' @seealso [wetlandToLCC()]
#' @export
prepInputs_CWIM <- function(to, url = .cwimURL, wetClasses = 1:4, wetThreshold = 0.5,
                            writeTo = NULL) {
  if (!inherits(to, "SpatRaster")) {
    stop("`to` must be a SpatRaster: it defines the grid the wetland layer is returned on")
  }
  stopifnot(is.numeric(wetThreshold), length(wetThreshold) == 1,
            wetThreshold > 0, wetThreshold <= 1)

  src <- url
  if (grepl("^https?://", url)) {
    src <- paste0("/vsicurl/", url)
    ## Without these GDAL lists the whole S3 "directory" before opening one file.
    old <- Sys.getenv(c("GDAL_DISABLE_READDIR_ON_OPEN", "CPL_VSIL_CURL_ALLOWED_EXTENSIONS"),
                      unset = NA)
    Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
               CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif")
    on.exit({
      for (nm in names(old)) {
        if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(old[nm]))
      }
    }, add = TRUE)
  }
  cw <- terra::rast(src)

  ## The study window in CWIM's CRS, widened by two target cells so reprojection has source
  ## data right to the edge. CWIM (EPSG:3979) and the LandR grids are both Lambert conformal
  ## conic but not always with the same latitude of origin, so the window is projected rather
  ## than reused as coordinates.
  win <- terra::ext(terra::project(terra::as.polygons(terra::ext(to), crs = terra::crs(to)),
                                   terra::crs(cw)))
  win <- terra::extend(win, 2 * max(terra::res(to)))
  x <- terra::crop(cw, win)

  ## A plain vector function through app(): `%in%` on a SpatRaster is not dispatched from
  ## package code in every terra version (1.9.46 hands `ifel()` a bare logical).
  wet <- terra::app(x, fun = function(v) as.integer(!is.na(v) & v %in% wetClasses))
  frac <- terra::project(wet, to, method = "average")
  out <- terra::ifel(frac >= wetThreshold, 1L, 0L)
  out <- terra::mask(out, to)
  names(out) <- "wetland"

  if (!is.null(writeTo)) {
    out <- terra::writeRaster(out, writeTo, overwrite = TRUE, datatype = "INT1U")
  }
  out
}

#' Add wetland classes to a land-cover map
#'
#' Marks the pixels a site layer calls wet with the NTEMS wetland codes: `wetTreedCode` (81)
#' where the land cover is treed, `wetCode` (80) otherwise. Pixels that are not wet are
#' returned unchanged, and so are water and any pixel already carrying a wetland code.
#'
#' "Treed" includes 240 by default: the disturbed-forest code LandR assigns to forest land
#' that is not currently treed, which is still forest ground.
#'
#' @param lcc Land cover: a `SpatRaster` or a numeric vector.
#' @param wet Site layer of the same shape: non-zero means wet, `0` or `NA` not. Typically
#'   [prepInputs_CWIM()].
#' @param treedClasses Land-cover codes that become `wetTreedCode` when wet.
#' @param keepClasses Codes never overwritten (default: water).
#' @param wetCode,wetTreedCode Codes for wet non-treed and wet treed pixels.
#'
#' @return `lcc`, with wet pixels recoded, in the same form it was given.
#'
#' @seealso [prepInputs_CWIM()]
#' @export
wetlandToLCC <- function(lcc, wet, treedClasses = c(210, 220, 230, 240),
                         keepClasses = 20, wetCode = 80L, wetTreedCode = 81L) {
  isRas <- inherits(lcc, "SpatRaster")
  if (isRas && inherits(wet, "SpatRaster")) {
    terra::compareGeom(lcc, wet)
  }
  v <- if (isRas) terra::values(lcc, mat = FALSE) else lcc
  w <- if (inherits(wet, "SpatRaster")) terra::values(wet, mat = FALSE) else wet
  if (length(v) != length(w)) {
    stop("`lcc` and `wet` must have the same number of cells (", length(v), " vs ", length(w), ")")
  }

  isWet <- !is.na(w) & w != 0 & !is.na(v)
  untouched <- v %in% c(keepClasses, wetCode, wetTreedCode)
  toTreed <- isWet & !untouched & v %in% treedClasses
  toWet <- isWet & !untouched & !v %in% treedClasses

  out <- v
  out[toTreed] <- wetTreedCode
  out[toWet] <- wetCode

  if (isRas) {
    r <- lcc
    terra::values(r) <- out
    r
  } else {
    out
  }
}
