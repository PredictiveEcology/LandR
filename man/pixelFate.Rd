% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohorts.R
\name{pixelFate}
\alias{pixelFate}
\title{Create or amend data to a pixelFateDT object}
\usage{
pixelFate(
  pixelFateDT,
  fate = NA_character_,
  pixelsRemoved = 0,
  runningPixelTotal = NA_integer_
)
}
\arguments{
\item{pixelFateDT}{A pixelFateDT data.table (has 3 columns, fate, pixelsRemoted, and
runningPixelTotal)}

\item{fate}{A character string (length 1) describing in words the change}

\item{pixelsRemoved}{A numeric indicating how many pixels were removed due to the
\code{fate}}

\item{runningPixelTotal}{an optional numeric with new, running total. If not supplied,
it will be calculated from the last row of pixelFateDT runningTotal minus the
\code{pixelsRemoved}}
}
\value{
A pixelFateDT object, updated with one extra row.
}
\description{
Create or amend data to a pixelFateDT object
}