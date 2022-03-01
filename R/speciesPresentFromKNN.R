#' Make a Species factor raster
#'
#' This will download all KNN layers in (Boreal) Forest of Canada, and make
#' a factor raster at resolution provided by \code{res} (larger is faster).
#' @param year Default (and only implemented) is 2011. This will download the 2011 KNN data layers
#' @param dPath A character string indicating where to download all the KNN layers
#' @param res The resolution (one dimension, in m) for the resulting raster
#' @param minPctCover An integer indicating what percent cover a species must have
#'   in a pixel to be considered present in that pixel.
#' @rdname speciesPresent
#' @return
#' A terra rast object with 2 layers, "speciesPresent" is a factor rast, with
#' a legend (i.e., it is numbers on a map, that correspond to a legend) and
#' "numberSpecies" which represents the number of species in each pixel.
#' @examples
#' # Make the dataset
#' speciesPresent <- speciesPresentFromKNN("~/data/KNN")
#'
#' \dontrun{
#'   # To upload this:
#'   speciesPresentRas <- raster::stack(speciesPresent)[[1]]
#'   fn <- "SpeciesPresentInCanadianForests.tif"
#'   writeRaster(speciesPresentRas, file = fn)
#'   out <- googledrive::drive_put(fn)
#'   driveID <- "1GrV5LjXS_N4iMwGYfBjUZ4elT2DkuKyJ"
#' }
#'
#' # Get species list
#' sa <- LandR::randomStudyArea(size = 1e11)
#' species <- LandR::speciesInStudyArea(sa)
#'
speciesPresentFromKNN <- function(year = 2011, dPath = asPath("."), res = 2000, minPctCover = 10) {
  studyAreaED <- Cache(prepInputs, url =  "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
                       destinationPath = dPath, #fun = quote(SA_ERIntersect(x = targetFilePath, studyArea)),
                       overwrite = FALSE)

  bf <- Cache(prepInputs, url = borealForestURL, fun = "forestOutline")

  opts <- options("reproducible.useTerra" = TRUE)
  on.exit(opts)
  studyAreaER <- Cache(prepInputs, url =  "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecodistrict_shp.zip",
                       destinationPath = dPath, fun = "terra::vect",
                       overwrite = TRUE)
  sa <- reproducible:::maskTo(studyAreaER, bf)
  sa <- reproducible:::projectTo(sa, sf::st_crs(bf))

  allForestedStk <- Cache(loadAndAggregateKNN, dPath, res, sa)
  allForestedStk <- round(allForestedStk, 0)
  allForestedStk[allForestedStk <= minPctCover] <- 0

  numSp <- sum(allForestedStk > 0)

  mat <- terra::values(allForestedStk)
  dt <- as.data.table(mat)
  dt[, pixel := 1:.N]
  dt2 <- melt(dt, measure.vars = setdiff(colnames(dt), "pixel"), na.rm = TRUE,
              id.vars = "pixel")
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
#' \code{speciesInStudyArea} defaults to use a url of a dataset uploaded
#' to googledrive that is from
#' Canadian Boreal Forests, but a different factor raster can be passed e.g.,
#' from \code{speciesPresentFromKNN}.
#'
#' @rdname speciesPresent
#' @return
#' A named list of length 2: \code{speciesRas} is a factor \code{RasterLayer}
#' and \code{speciesList} is a character string containing the unique, sorted
#' species on the \code{speciesRas}, for convenience.
#' @export
#' @param studyArea a vector map (e.g., SpatialPolygonsDataFrame)
#' @param url A url to get a speciesPresence raster e.g., from \code{peciesPresentFromKNN}
#' @param speciesPresentRas A factor raster where the character string is a string of
#'   species names, separated by 2 underscores, sorted alphabetically. Can be produced
#'   with \code{speciesPresentFromKNN}
#'
speciesInStudyArea <- function(studyArea, url = NULL, speciesPresentRas = NULL) {
  if (is.null(speciesPresentRas)) {
    if (is.null(url))
      url <- "https://drive.google.com/file/d/1GrV5LjXS_N4iMwGYfBjUZ4elT2DkuKyJ"
    speciesPres <- preProcess(url = url)
    speciesPresRas <- raster::raster(speciesPres$targetFilePath)
  }

  bb <- postProcess(x = speciesPresRas, studyArea = studyArea)
  speciesCommunities <- na.omit(factorValues2(bb, bb[], att = "category"))
  species <- as.character(speciesCommunities)
  species <- unique(unlist(strsplit(species, "__")))
  return(list(speciesRas = bb, speciesList = species))
}

forestOutline <- function(x) {
  x1 <- terra::vect(x)
  bf2 <- terra::simplifyGeom(x1, tolerance = 5000)
  bf3 <- terra::makeValid(bf2)
  bf4 <- terra::aggregate(bf3)
  bf5 <- terra::makeValid(bf4)
  bf6 <- terra::buffer(bf5, 6000)
  bf7 <- terra::aggregate(bf6)
  bf8 <- terra::buffer(bf7, -6000)
}

borealForestURL <- "https://d278fo2rk9arr5.cloudfront.net/downloads/boreal.zip"

SA_ERIntersect <- function(x, studyArea) {
  x <- sf::st_read(x)
  sa_sf <- sf::st_as_sf(studyArea)
  ecoregions <- sf::st_transform(x, sf::st_crs(sa_sf))
  studyAreaER <- sf::st_intersects(ecoregions, sa_sf, sparse = FALSE)
  terra::vect(sf::as_Spatial(ecoregions[studyAreaER,]))
}


loadAndAggregateKNN <- function(dPath, res, sa) {
  ll <- terra::rast(loadkNNSpeciesLayers(dPath))
  llCoarse <- terra::aggregate(ll, res/250)
  postProcessTerra(from = llCoarse, to = sa, method = "near")
}
