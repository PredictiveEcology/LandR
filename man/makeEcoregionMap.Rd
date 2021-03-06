% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoregions.R
\name{makeEcoregionMap}
\alias{makeEcoregionMap}
\title{Make the \code{ecoregionMap} raster}
\usage{
makeEcoregionMap(ecoregionFiles, pixelCohortData)
}
\arguments{
\item{ecoregionFiles}{A list with two objects: the \code{ecoregionMap} and a table summarizing
its information per \code{pixelID}.}

\item{pixelCohortData}{The full \code{cohortData} \code{data.table}}
}
\value{
A raster with ecoregion codes.
}
\description{
Creates a raster of ecoregion codes per pixel.
Unlike \code{ecoregionProducer}, this fills the raster with pixel information contained in
\code{pixelCohortData}.
}
