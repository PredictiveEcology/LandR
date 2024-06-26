% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speciesPresentFromKNN.R
\name{speciesPresentFromKNN}
\alias{speciesPresentFromKNN}
\alias{speciesInStudyArea}
\title{Make a species factor raster}
\usage{
speciesPresentFromKNN(
  year = 2011,
  dPath = asPath("."),
  res = 2000,
  minPctCover = 10
)

speciesInStudyArea(
  studyArea,
  url = NULL,
  speciesPresentRas = NULL,
  dPath = getOption("reproducible.destinationPath")
)
}
\arguments{
\item{year}{Default (and only implemented) is 2011. This will download the 2011 KNN data layers}

\item{dPath}{Passed to \code{destinationPath} in \code{preProcess}.}

\item{res}{The resolution (one dimension, in m) for the resulting raster}

\item{minPctCover}{An integer indicating what percent cover a species must have
in a pixel to be considered present in that pixel.}

\item{studyArea}{A \verb{SpatialPolygons*} object used as the principle study region,
passed to \code{\link[reproducible:prepInputs]{reproducible::prepInputs()}}.}

\item{url}{A url to get a \code{speciesPresence} raster e.g., from \code{peciesPresentFromKNN}}

\item{speciesPresentRas}{A factor raster where the character string is a string of
species names, separated by 2 underscores, sorted alphabetically. Can be produced
with \code{speciesPresentFromKNN}}
}
\value{
A \code{SpatRaster} object with 2 layers: \code{"speciesPresent"} is a factor, with
a legend (i.e., it is numbers on a map, that correspond to a legend) and
\code{"numberSpecies"} which represents the number of species in each pixel.

A named list of length 2: \code{speciesRas} is a factor \code{RasterLayer}
and \code{speciesList} is a character string containing the unique, sorted
species on the \code{speciesRas}, for convenience.
}
\description{
This will download all KNN layers in (Boreal) Forest of Canada, and make
a factor raster at resolution provided by \code{res} (larger is faster).

\code{speciesInStudyArea} defaults to use a url of a dataset uploaded to Google Drive that is
from Canadian Boreal Forests, but a different factor raster can be passed e.g.,
from \code{speciesPresentFromKNN}.
}
\examples{
\dontrun{
if (requireNamespace("googledrive", quietly = TRUE)) {
  # Make the dataset
  speciesPresent <- speciesPresentFromKNN(dPath = "~/data/KNN")

  # To upload this:
  speciesPresentRas <- terra::rast(speciesPresent)[[1]]
  fn <- "SpeciesPresentInCanadianForests.tif"
  writeRaster(speciesPresentRas, file = fn)
  zipFn <- gsub(".tif", ".zip", fn)
  zip(files = dir(pattern = fn), zipFn)
  out <- googledrive::drive_put(zipFn)
  driveID <- "1Oj78jJBeha5L6XDBBdWDAfimgNjYc9UD"

  # Get species list
  sa <- LandR::randomStudyArea(size = 1e11)
  species <- LandR::speciesInStudyArea(sa)
}
}

}
