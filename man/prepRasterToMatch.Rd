% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepInputObjects.R
\name{prepRasterToMatch}
\alias{prepRasterToMatch}
\title{Create \code{rasterToMatch} and \code{rasterToMatchLarge}}
\usage{
prepRasterToMatch(
  studyArea,
  studyAreaLarge,
  rasterToMatch,
  rasterToMatchLarge,
  destinationPath,
  templateRas,
  studyAreaName,
  cacheTags
)
}
\arguments{
\item{studyArea}{A \verb{SpatialPolygons*} object used as the principle study region,
passed to \code{\link[reproducible:prepInputs]{reproducible::prepInputs()}}.}

\item{studyAreaLarge}{same as \code{studyArea}, but larger and completely
covering it.}

\item{rasterToMatch}{A \code{RasterLayer} objects to use as the template for all subsequent
raster operations (i.e., the one used throughout the simulation).}

\item{rasterToMatchLarge}{the \code{rasterToMatchLarge} object from \code{sim}}

\item{destinationPath}{path to data directory where objects will be downloaded or saved to}

\item{templateRas}{a template raster used to make \code{rasterToMatch}
and/or \code{rasterToMatchLarge}. Must match \code{studyAreaLarge}.}

\item{studyAreaName}{character string giving the study area name.}

\item{cacheTags}{a character vector of strings to pass to \code{Cache(userTags)}}
}
\description{
\code{rasterToMatch} and \code{rasterToMatchLarge} raster layers are created
from \code{studyArea} and \code{studyAreaLarge} polygons (respectively)
using a template raster (often \code{rawBiomassMap})
}
