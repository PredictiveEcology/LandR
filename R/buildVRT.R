#' Check if a SpatRaster should be converted to tiled virtual raster
#'
#' @param ras a file-backed raster
#' @param bitExponent 2 exponent this number sets a threshold whereby rasters
#' containing more cells than the threshold are converted to virtual rasters via tiles
#' @param destinationPath the output directory for raster, if written
#' @param writeTo filename - extension will be modified to .vrt if tiled
#' @param overwrite overte files

#' @returns the raster as a VRT
#' @importFrom tools file_path_sans_ext
#' @importFrom terra makeTiles vrt ncell sources
#' @importFrom reproducible checkPath .suffix
#'
#' @export
buildVRT <- function(ras, bitExponent = 30,
                     destinationPath = NULL,  writeTo = NULL,
                     overwrite = NULL) {

  #confirm if raster needs to be written
  if (ncell(ras) > 2^bitExponent) {
    if (!is.null(writeTo)) {
      rasBasename <- writeTo
    } else if (!inMemory(ras)) {
      rasBasename <- basename(sources(ras))
    }
    filenameNoExt <- file_path_sans_ext(rasBasename)

    message("creating tiles for raster due to raster size...")
    #this will make at least 2 tiles when the raster is anywhere close 2^30
    ndiv <- ceiling(ncell(ras)/c(2^bitExponent)) #with equal columns, 4 tiles minimum
    nrows <- ceiling(nrow(ras)/ndiv)
    ncols <- ceiling(ncol(ras)/ndiv)


    tileNames <- file.path(destinationPath, .suffix(writeTo, "_")) #tiles will be _1, _2, etc
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
