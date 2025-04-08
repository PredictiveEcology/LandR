#' Conditionally Create a Tiled Virtual Raster (VRT)
#'
#' @description
#' Evaluates the size (number of cells) of an input `SpatRaster`. If it exceeds
#' a specified threshold, the function converts it into a tiled Virtual Raster
#' (.vrt). Otherwise, it can optionally write the raster to a standard file format.
#' This is primarily useful for managing very large raster datasets that might
#' exceed memory limits or benefit from tiled access patterns.
#'
#' @param ras A `terra::SpatRaster` object. While the function works with
#'   in-memory rasters, VRT creation is most relevant for large, file-backed
#'   rasters.
#' @param bitExponent An integer defining the cell count threshold. If
#'   `terra::ncell(ras)` exceeds `2^bitExponent`, the raster will be processed
#'   as a VRT. Defaults to `30` (approximately 1 billion cells).
#' @param destinationPath Character string. The directory path where output files
#'   (tiles, VRT, or directly written raster) should be saved. If `NULL` (the default)
#'   and files need to be written (due to tiling or `writeTo` being set),
#'   the function attempts to write relative to the current working directory or
#'   based on `writeTo` if it includes a path. It is strongly recommended to
#'   provide an explicit path. The directory will be created if it doesn't exist.
#' @param writeTo Optional character string. The base filename (without extension)
#'   for the output.
#'   - If tiling occurs (raster exceeds threshold), this name is used for the
#'     `.vrt` file and as a prefix for the individual tile files (e.g.,
#'     `filename_1.tif`, `filename_2.tif`, ...). The extension `.vrt` is appended
#'     automatically to the VRT filename.
#'   - If tiling does *not* occur, this filename (with its original or appropriate
#'     extension) is used when writing the raster directly via `terra::writeRaster`.
#'   - If `NULL` and tiling occurs, a base name is derived from `basename(sources(ras))`.
#'     if the raster is on disk - otherwise `tempfile()` is used
#'   - If `NULL` and tiling does *not* occur, the raster is *not* written to disk
#'     (unless it was already file-backed).
#' @param overwrite Logical. Controls whether existing output files (tiles, .vrt file,
#'   or directly written raster file) should be overwritten. Passed directly to
#'   `terra::makeTiles`, `terra::vrt`, and `terra::writeRaster`. If `NULL` (the default),
#'   the behavior depends on the defaults of the underlying `terra` functions
#'   (typically `FALSE`, preventing overwrite).
#'
#' @details
#' Large raster datasets can be cumbersome to work with. Virtual Rasters (VRT)
#' provide a way to handle them efficiently by defining a virtual dataset composed
#' of multiple smaller tiles.
#'
#' This function automates the process:
#' 1. It checks if the total number of cells in `ras` (`terra::ncell(ras)`)
#'    is greater than `2^bitExponent`.
#' 2. If the threshold is exceeded:
#'    - It calculates the necessary number of tiles (`ndiv`) to ensure each tile is
#'      roughly below the cell threshold.
#'    - It calls `terra::makeTiles` to generate the individual raster tile files
#'      (saved in `destinationPath` with names derived from `writeTo`).
#'    - It calls `terra::vrt` to create a single `.vrt` file in `destinationPath`
#'      that references these tiles.
#'    - The function returns a `SpatRaster` object pointing to this `.vrt` file.
#' 3. If the threshold is *not* exceeded:
#'    - If `writeTo` is provided, the function calls `terra::writeRaster` to save
#'      the input `ras` directly to `file.path(destinationPath, writeTo)`. The
#'      returned `SpatRaster` points to this new file.
#'    - If `writeTo` is `NULL`, the original `ras` object is returned unmodified
#'      (no files are written by this function).
#'

#' @return A `terra::SpatRaster` object.
#'   - If tiling occurred, it represents the VRT dataset (`.vrt` file).
#'   - If tiling did not occur but `writeTo` was provided, it represents the
#'     raster written directly to disk.
#'   - If tiling did not occur and `writeTo` was `NULL`, it is the original
#'     input `ras` object.
#' @importFrom tools file_path_sans_ext
#' @importFrom terra makeTiles vrt ncell sources writeRaster nrow ncol rast inMemory
#' @importFrom reproducible checkPath .suffix
#'
#' @export
buildVRT <- function(ras, bitExponent = 30,
                     destinationPath = NULL,  writeTo = NULL,
                     overwrite = NULL) {

  #confirm if raster needs to be written
  if (ncell(ras) > 2^bitExponent) {

    if (is.null(writeTo)) {
      if (!inMemory(ras)) {
        rasBasename <- basename(sources(ras))
      } else {
        rasBasename <- tempfile(pattern = ".tif")
      }
    } else {
      rasBasename <- writeTo
    }
    filenameNoExt <- file_path_sans_ext(rasBasename)

    message("creating tiles for raster due to raster size...")
    #this will make at least 2 tiles when the raster is anywhere close 2^30
    ndiv <- ceiling(ncell(ras)/c(2^bitExponent)) #with equal columns, 4 tiles minimum
    nrows <- ceiling(nrow(ras)/ndiv)
    ncols <- ceiling(ncol(ras)/ndiv)

    if (is.null(destinationPath)) {
      destinationPath <- "." #get the folder
    }

    tileNames <- file.path(destinationPath, .suffix(rasBasename, "_")) #tiles will be _1, _2, etc
    ras <- terra::makeTiles(ras, filename = tileNames,
                            y = c(nrows, nrows),
                            na.rm = TRUE, overwrite = overwrite)
    vrtFilename <- file.path(destinationPath, paste0(filenameNoExt, ".vrt"))
    ras <- terra::vrt(ras, filename = vrtFilename, overwrite = overwrite)

  } else if (!is.null(writeTo)) {
    ras <- writeRaster(ras, filename = file.path(destinationPath, writeTo))
  }
  return(ras)

}
