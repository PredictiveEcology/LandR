% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepSpeciesLayers.R
\name{NAcover2zero}
\alias{NAcover2zero}
\title{Convert \code{NA} values in \code{speciesLayers} to zeros}
\usage{
NAcover2zero(speciesLayers, rasterToMatch)
}
\arguments{
\item{speciesLayers}{A \code{RasterStack} or \code{RasterLayer} that
should contain species cover data in the study area}

\item{rasterToMatch}{A \code{RasterLayer} objects to use as the template for all subsequent
raster operations (i.e., the one used throughout the simulation).}
}
\value{
the \code{speciesLayers} with \code{0} in pixels that had \code{NA}
}
\description{
Pixels that are \code{NA} but are inside \code{rasterToMatch} may need to be converted to \code{0},
as they can could still potentially be forested
}
