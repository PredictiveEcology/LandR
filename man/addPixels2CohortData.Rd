% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohorts.R
\name{addPixels2CohortData}
\alias{addPixels2CohortData}
\title{Generate \code{cohortData} table per pixel:}
\usage{
addPixels2CohortData(
  cohortData,
  pixelGroupMap,
  cohortDefinitionCols = c("pixelGroup", "age", "speciesCode"),
  doAssertion = getOption("LandR.assertions", TRUE)
)
}
\arguments{
\item{cohortData}{A \code{data.table} with columns: \code{pixelGroup}, \code{ecoregionGroup},
\code{speciesCode}, and optionally \code{age}, \code{B}, \code{mortality}, \code{aNPPAct}, and \code{sumB}.}

\item{pixelGroupMap}{A \code{RasterLayer} with pixel values equal to a pixel group
number that corresponds exactly to \code{pixelGroup} column in \code{cohortData}.}

\item{cohortDefinitionCols}{the columns in \code{cohortData} that define unique cohorts}

\item{doAssertion}{A logical indicating whether some internal tests should be run to
ensure the function is running correctly.
Default: \code{getOption("LandR.assertions", TRUE)}.}
}
\value{
An expanded \code{cohortData} \code{data.table} with a new \code{pixelIndex} column.
}
\description{
Generate \code{cohortData} table per pixel:
}
