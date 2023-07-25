utils::globalVariables(c(
  "allPres", "allPresFac", "pixel", "variable"
))

maskTo <- utils::getFromNamespace("maskTo", "reproducible")
projectTo <- utils::getFromNamespace("projectTo", "reproducible")

#' Make a species factor raster
#'
#' This will download all KNN layers in (Boreal) Forest of Canada, and make
#' a factor raster at resolution provided by `res` (larger is faster).
#'
#' @param year Default (and only implemented) is 2011. This will download the 2011 KNN data layers
#'
#' @param dPath A character string indicating where to download all the KNN layers
#'
#' @param res The resolution (one dimension, in m) for the resulting raster
#'
#' @param minPctCover An integer indicating what percent cover a species must have
#'   in a pixel to be considered present in that pixel.
#'
#' @return A `SpatRaster` object with 2 layers: `"speciesPresent"` is a factor, with
#' a legend (i.e., it is numbers on a map, that correspond to a legend) and
#' `"numberSpecies"` which represents the number of species in each pixel.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("googledrive", quietly = TRUE)) {
#'   # Make the dataset
#'   speciesPresent <- speciesPresentFromKNN(dPath = "~/data/KNN")
#'
#'   # To upload this:
#'   speciesPresentRas <- raster::stack(speciesPresent)[[1]]
#'   fn <- "SpeciesPresentInCanadianForests.tif"
#'   writeRaster(speciesPresentRas, file = fn)
#'   zipFn <- gsub(".tif", ".zip", fn)
#'   zip(files = dir(pattern = fn), zipFn)
#'   out <- googledrive::drive_put(zipFn)
#'   driveID <- "1Oj78jJBeha5L6XDBBdWDAfimgNjYc9UD"
#'
#'   # Get species list
#'   sa <- LandR::randomStudyArea(size = 1e11)
#'   species <- LandR::speciesInStudyArea(sa)
#' }
#' }
#'
#' @export
#' @importFrom data.table as.data.table setorderv melt.data.table
#' @importFrom reproducible asPath Cache
#' @importFrom sf st_crs
#' @rdname speciesPresent
speciesPresentFromKNN <- function(year = 2011, dPath = asPath("."), res = 2000, minPctCover = 10) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop("Suggested package 'terra' not installed.\n",
         "Install it using `install.packages('terra')`.")
  }

  studyAreaED <- Cache(
    prepInputs,
    url =  "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    destinationPath = dPath, #fun = quote(SA_ERIntersect(x = targetFilePath, studyArea)),
    overwrite = FALSE
  )

  bf <- Cache(prepInputs, url = borealForestURL, fun = "forestOutline")

  opts <- options("reproducible.useTerra" = TRUE)
  on.exit(opts)
  studyAreaER <- Cache(
    prepInputs,
    url =  "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecodistrict_shp.zip",
    destinationPath = dPath,
    fun = "terra::vect",
    overwrite = TRUE
  )
  sa <- maskTo(studyAreaER, bf)
  sa <- projectTo(sa, sf::st_crs(bf))

  allForestedStk <- Cache(loadAndAggregateKNN, dPath, res, sa)
  allForestedStk <- round(allForestedStk, 0)
  allForestedStk[allForestedStk <= minPctCover] <- 0

  numSp <- sum(allForestedStk > 0)

  mat <- terra::values(allForestedStk)
  dt <- as.data.table(mat)
  dt[, pixel := 1:.N]
  dt2 <- melt(dt, measure.vars = setdiff(colnames(dt), "pixel"), na.rm = TRUE, id.vars = "pixel")
  dt2 <- dt2[value != 0]
  setorderv(dt2, c("pixel", "variable"))
  dt3 <- dt2[, list(allPres = paste(variable, collapse = "__")), by = "pixel"]
  dt3[, allPresFac := factor(allPres)]

  # Create a new empty rast
  speciesPres <- terra::rast(allForestedStk[[1]])
  # fill it with the integer values
  speciesPres[dt3$pixel] <- as.integer(dt3$allPresFac)
  names <- unique(dt3$allPresFac)
  numerics <- as.integer(names)
  # assign the levels
  levels(speciesPres) <- data.frame(ID = numerics, spGroup = names)

  return(c(speciesPres, numSp))
}

#' Get species list in a given study area for boreal forest of Canada
#'
#' `speciesInStudyArea` defaults to use a url of a dataset uploaded to Google Drive that is
#' from Canadian Boreal Forests, but a different factor raster can be passed e.g.,
#' from `speciesPresentFromKNN`.
#'
#' @template studyArea
#'
#' @param url A url to get a `speciesPresence` raster e.g., from `peciesPresentFromKNN`
#'
#' @param speciesPresentRas A factor raster where the character string is a string of
#'   species names, separated by 2 underscores, sorted alphabetically. Can be produced
#'   with `speciesPresentFromKNN`
#'
#' @return A named list of length 2: `speciesRas` is a factor `RasterLayer`
#' and `speciesList` is a character string containing the unique, sorted
#' species on the `speciesRas`, for convenience.
#'
#' @export
#' @importFrom pemisc factorValues2
#' @importFrom raster raster
#' @importFrom reproducible preProcess
#' @rdname speciesPresent
speciesInStudyArea <- function(studyArea, url = NULL, speciesPresentRas = NULL) {
  if (is.null(speciesPresentRas)) {
    if (is.null(url)) {
      url <- "https://drive.google.com/file/d/1Oj78jJBeha5L6XDBBdWDAfimgNjYc9UD/"
    }
    speciesPres <- preProcess(url = url)
    speciesPresRas <- raster::raster(speciesPres$targetFilePath)
  }

  if (getOption("reproducible.useTerra", TRUE) && requireNamespace("terra")) {
    bb <- postProcessTerra(speciesPresRas, studyArea = studyArea)
  } else {
    bb <- postProcess(x = speciesPresRas, studyArea = studyArea)
  }
  rasLevs <- raster::levels(speciesPresRas)[[1]]
  rasLevs <- rasLevs[rasLevs$ID %in% na.omit(getValues(bb)), ]
  levels(bb) <- rasLevs
  bb <- raster::deratify(bb)
  speciesCommunities <- na.omit(factorValues2(bb, bb[], att = "category"))
  species <- as.character(speciesCommunities)
  species <- unique(unlist(strsplit(species, "__")))
  return(list(speciesRas = bb, speciesList = species))
}

#' @keywords internal
forestOutline <- function(x) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop("Suggested package 'terra' not installed.\n",
         "Install it using `install.packages('terra')`.")
  }

  x1 <- terra::vect(x)
  bf2 <- terra::simplifyGeom(x1, tolerance = 5000)
  bf3 <- terra::makeValid(bf2)
  bf4 <- terra::aggregate(bf3)
  bf5 <- terra::makeValid(bf4)
  bf6 <- terra::buffer(bf5, 6000)
  bf7 <- terra::aggregate(bf6)
  bf8 <- terra::buffer(bf7, -6000)
}

## TODO: randomized URL changes
borealForestURL <- "https://d278fo2rk9arr5.cloudfront.net/downloads/boreal.zip"

#' @importFrom sf as_Spatial st_as_sf st_intersects st_read st_transform
#' @keywords internal
SA_ERIntersect <- function(x, studyArea) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop("Suggested package 'terra' not installed.\n",
         "Install it using `install.packages('terra')`.")
  }

  x <- sf::st_read(x)
  sa_sf <- sf::st_as_sf(studyArea)
  ecoregions <- sf::st_transform(x, sf::st_crs(sa_sf))
  studyAreaER <- sf::st_intersects(ecoregions, sa_sf, sparse = FALSE)
  terra::vect(sf::as_Spatial(ecoregions[studyAreaER,]))
}

#' @importFrom reproducible postProcessTerra
#' @keywords internal
loadAndAggregateKNN <- function(dPath, res, sa) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop("Suggested package 'terra' not installed.\n",
         "Install it using `install.packages('terra')`.")
  }

  ll <- terra::rast(loadkNNSpeciesLayers(dPath))
  llCoarse <- terra::aggregate(ll, res/250)
  postProcessTerra(from = llCoarse, to = sa, method = "near")
}
