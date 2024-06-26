% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohorts.R
\name{makeAndCleanInitialCohortData}
\alias{makeAndCleanInitialCohortData}
\title{Generate initial \code{cohortData} table}
\usage{
makeAndCleanInitialCohortData(
  inputDataTable,
  sppColumns,
  imputeBadAgeModel = quote(lme4::lmer(age ~ B * speciesCode + cover * speciesCode + (1 |
    initialEcoregionCode))),
  minCoverThreshold,
  doAssertion = getOption("LandR.assertions", TRUE),
  doSubset = TRUE
)
}
\arguments{
\item{inputDataTable}{A \code{data.table} with columns described above.}

\item{sppColumns}{A vector of the names of the columns in \code{inputDataTable} that
represent percent cover by species, rescaled to sum up to 100\%\%.}

\item{imputeBadAgeModel}{statistical model used to impute ages in pixels with missing
data or with cover == 0. If set to NULL no imputation will be attempted, and pixels with
missing age are excluded.}

\item{minCoverThreshold}{minimum total cover percentage necessary to consider the pixel
vegetated, or a cohort present in a pixel.}

\item{doAssertion}{A logical indicating whether some internal tests should be run to
ensure the function is running correctly.
Default: \code{getOption("LandR.assertions", TRUE)}.}

\item{doSubset}{Turns on/off subsetting. Defaults to \code{TRUE}.}
}
\value{
a \code{cohortData} \code{data.table} with attribute \code{"imputedPixID"}
(a vector of pixel IDs that suffered imputation).
}
\description{
Takes a single \code{data.table} input, which has the following columns in addition to
others that will be labelled with species name, and contain percent cover of each:
}
\details{
\itemize{
\item \code{pixelIndex} (integer)
\item \code{age} (integer)
\item \code{logAge} (numeric)
\item \code{initialEcoregionCode} (factor)
\item \code{totalBiomass} (integer)
\item \code{lcc} (integer)
\item \code{rasterToMatch} (integer)
\item \code{speciesCode} (factor)
\item \code{cover} (integer)
\item \code{coverOrig} (integer)
\item \code{B} (integer)
}

Several data correction/imputation operations are also performed. Namely, age is imputed
in pixels where age data is missing (but not cover) and where \code{cover == 0} but \code{age > 0},
total biomass is zeroed if \code{age == 0}, and age is zeroed if \code{biomass == 0}.
}
\author{
Eliot McIntire
}
