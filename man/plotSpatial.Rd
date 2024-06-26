% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terra-raster-utils.R
\name{plotSpatial}
\alias{plotSpatial}
\title{Create a \code{ggplot} of a raster or \code{sf} object.}
\usage{
plotSpatial(x, plotTitle, limits = NULL, field = NULL)
}
\arguments{
\item{x}{\code{SpatRaster}, \code{RasterLayer}, \code{SpatVector} or \code{sf} object}

\item{plotTitle}{character. A title for the plot passed to \code{ggplot::labs(title = plotTitle)}.}

\item{limits}{TODO}

\item{field}{character. If \code{x} is \code{sf} or \code{SpatVector}, a field to plot.}
}
\description{
Can be used with \code{SpaDES.core::Plots}.
}
