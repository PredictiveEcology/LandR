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
#'     This may fail or produce unexpected results for in-memory rasters.
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
#' The `reproducible::checkPath` function is used to ensure `destinationPath` exists,
#' creating it if necessary.
#'
#' @return A `terra::SpatRaster` object.
#'   - If tiling occurred, it represents the VRT dataset (`.vrt` file).
#'   - If tiling did not occur but `writeTo` was provided, it represents the
#'     raster written directly to disk.
#'   - If tiling did not occur and `writeTo` was `NULL`, it is the original
#'     input `ras` object.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom terra makeTiles vrt ncell sources writeRaster nrow ncol rast inMemory
#' @importFrom reproducible checkPath .suffix
#'
#' @export
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Create a temporary directory
#' tempDir <- tempdir()
#'
#' # --- Example 1: Small raster (no tiling) ---
#' rasSmall <- rast(nrows = 100, ncols = 100, vals = 1:10000)
#'
#' # Case 1a: No writing
#' rasOut1a <- buildVRT(rasSmall, bitExponent = 10, destinationPath = tempDir)
#' print(rasOut1a) # Original in-memory raster returned
#' terra::sources(rasOut1a) # No source file
#'
#' # Case 1b: Writing directly (no tiling)
#' rasOut1b <- buildVRT(rasSmall, bitExponent = 10, destinationPath = tempDir,
#'                     writeTo = "small_raster.tif", overwrite = TRUE)
#' print(rasOut1b)
#' terra::sources(rasOut1b) # Points to small_raster.tif
#' file.exists(file.path(tempDir, "small_raster.tif")) # TRUE
#'
#' # --- Example 2: Large raster (triggers tiling) ---
#' # Create a larger dummy raster (e.g., > 2^15 cells)
#' rasLarge <- rast(nrows = 200, ncols = 200) # 40,000 cells
#' values(rasLarge) <- 1:ncell(rasLarge)
#' # Save it to disk first, as makeTiles works best with file-backed rasters
#' largeRasPath <- file.path(tempDir, "large_temp.tif")
#' writeRaster(rasLarge, largeRasPath, overwrite = TRUE)
#' rasLargeFile <- rast(largeRasPath)
#'
#' # Build VRT (using a small bitExponent for demonstration)
#' rasOut2 <- buildVRT(rasLargeFile, bitExponent = 15, # Threshold = 32768 cells
#'                    destinationPath = tempDir,
#'                    writeTo = "large_raster_vrt", overwrite = TRUE)
#'
#' print(rasOut2) # SpatRaster representing the VRT
#' terra::sources(rasOut2) # Points to large_raster_vrt.vrt
#'
#' # Check created files
#' file.exists(file.path(tempDir, "large_raster_vrt.vrt"))   # TRUE
#' # Tile files (names might vary slightly based on terra version/exact dims)
#' list.files(tempDir, pattern = "large_raster_vrt_.*\\.tif")
#'
#' # Clean up
#' unlink(tempDir, recursive = TRUE)
#' }
#'
buildVRT <- function(ras, bitExponent = 30,
                     destinationPath = NULL, writeTo = NULL,
                     overwrite = NULL) {

  # Default overwrite to FALSE if NULL for safety, consistent with terra typically
  # Although terra::vrt/makeTiles might handle NULL, explicitly setting avoids ambiguity
  if (is.null(overwrite)) {
    overwrite <- FALSE
  }

  # Determine if tiling is needed based on cell count threshold
  needsTiling <- terra::ncell(ras) > (2^bitExponent)

  if (needsTiling) {
    message("Raster exceeds cell threshold (", prettyNum(2^bitExponent, big.mark = ","),
            " cells). Creating tiled VRT...")

    # Ensure destination path exists if writing tiles/VRT
    if (is.null(destinationPath)) {
      # If no destinationPath, attempt to use current dir or path from writeTo
      # However, it's cleaner if the user provides it. Warn if ambiguous.
      if (!is.null(writeTo) && dirname(writeTo) != ".") {
        destinationPath <- dirname(writeTo)
        writeTo <- basename(writeTo) # Ensure writeTo is just the basename
        warning("`destinationPath` not provided. Using directory from `writeTo`: ", destinationPath)
      } else {
        destinationPath <- tempdir()
        warning("`destinationPath` not provided. Using temp directory: ", destinationPath)
      }
    }
    destinationPath <- reproducible::checkPath(destinationPath, create = TRUE)

    # Determine base filename for outputs
    if (!is.null(writeTo)) {
      rasBasename <- tools::file_path_sans_ext(basename(writeTo)) # Use provided base name
    } else if (!terra::inMemory(ras)) {
      # Derive from source file if available and writeTo not given
      rasBasename <- tools::file_path_sans_ext(basename(terra::sources(ras)[1]))
      message("`writeTo` not provided. Deriving base name from source: ", rasBasename)
    } else {
      # Cannot reliably determine name if in memory and writeTo is NULL
      #TODO: message
      rasBasename <- basename(tempfile())
    }

    # Define VRT filename
    vrtFilename <- file.path(destinationPath, paste0(rasBasename, ".vrt"))

    # Define base filename pattern for tiles
    # Use .suffix to handle potential existing files if overwrite = FALSE
    tileBase <- file.path(destinationPath, rasBasename)
    # Note: terra::makeTiles adds underscores and numbers automatically.
    # We provide the base path/name pattern.
    # Example: tileNames could be path/to/output/basename_ (terra adds 1.tif, 2.tif...)
    # Let's just pass the base filename to makeTiles, it handles suffixing.
    tileFilenamePattern <- file.path(destinationPath, rasBasename)

    # Calculate approximate number of divisions needed
    # Aim for at least 2 divisions if close to threshold for better tiling potential
    ndiv <- ceiling(terra::ncell(ras) / (2^bitExponent))
    # Ensure at least 2 divisions for meaningful tiling, unless raster is extremely wide/tall
    ndiv <- max(2, ndiv)

    # Calculate tile dimensions (equal number of rows/cols per tile, approximately)
    # Note: terra::makeTiles parameter 'y' expects c(nrows, ncols) tile dimensions
    # Alternative strategy: divide raster into ndiv x ndiv grid (approx)
    nrowsTile <- ceiling(terra::nrow(ras) / ndiv)
    ncolsTile <- ceiling(terra::ncol(ras) / ndiv)

    message("Creating tiles with approximate dimensions: ", nrowsTile, " rows x ", ncolsTile, " cols.")

    # Create the tiles on disk
    # filename arg in makeTiles is the *pattern* for output tiles
    rasTiles <- terra::makeTiles(ras, filename = paste0(tileFilenamePattern, ".tif"), # Ensure extension for tiles
                                 y = c(nrowsTile, ncolsTile),
                                 na.rm = TRUE, # Often useful
                                 overwrite = overwrite)

    # Create the VRT file referencing the tiles
    # rasTiles is now a vector of the tile filenames
    rasOut <- terra::vrt(rasTiles, filename = vrtFilename, overwrite = overwrite)
    message("VRT created: ", vrtFilename)

  } else {
    # Raster is below threshold, no tiling needed
    if (!is.null(writeTo)) {
      # Write the raster directly if writeTo is specified
      if (is.null(destinationPath)) {
        if (dirname(writeTo) != ".") {
          destinationPath <- dirname(writeTo)
          writeTo <- basename(writeTo)
          warning("`destinationPath` not provided. Using directory from `writeTo`: ", destinationPath)
        } else {
          destinationPath <- getwd()
          warning("`destinationPath` not provided. Using current working directory: ", destinationPath)
        }
      }
      destinationPath <- reproducible::checkPath(destinationPath, create = TRUE)
      outFilename <- file.path(destinationPath, writeTo)
      message("Raster below threshold. Writing directly to: ", outFilename)
      rasOut <- terra::writeRaster(ras, filename = outFilename, overwrite = overwrite)
    } else {
      # No tiling, no writing requested - return original raster
      message("Raster below threshold. No output file requested (`writeTo` is NULL). Returning original SpatRaster.")
      rasOut <- ras
    }
  }

  return(rasOut)
}
