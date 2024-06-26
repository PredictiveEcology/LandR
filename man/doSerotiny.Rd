% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serotiny-resprouting.R
\name{doSerotiny}
\alias{doSerotiny}
\title{Activate serotiny after a (fire) disturbance}
\usage{
doSerotiny(
  burnedPixelCohortData,
  postFirePixelCohortData,
  postFireRegenSummary = NULL,
  species,
  sufficientLight,
  speciesEcoregion,
  currentTime,
  treedFirePixelTableSinceLastDisp,
  calibrate = FALSE
)
}
\arguments{
\item{burnedPixelCohortData}{An expanded \code{cohortData} \code{data.table} with pixel-level
cohort information on burnt pixels and the following (optional) columns:
\code{severity} - fire severity in that pixel calculated based on fire behaviour properties;
\code{firetolerance} - species-level fire tolerance;
\code{severityToleranceDif} - the difference between \code{severity} and \code{firetolerance}.}

\item{postFirePixelCohortData}{an empty \code{cohortData}-like table with columns:
\code{age}, \code{B}, \code{mortality}, \code{aNPPAct}, and \code{sumB} removed; and \code{pixelIndex} added.}

\item{postFireRegenSummary}{a data.table summarizing for which species serotiny/resprouting were
activated and in how many pixels, for each year. Only necessary if \code{calibrate = TRUE}.}

\item{species}{a \code{data.table} with species traits such as longevity, shade tolerance, etc.
Must have column \code{speciesCode}, with species names/IDs. The following is a
list of default trait columns:
\itemize{
\item "species" same as "speciesCode" -- species ID name
\item "speciesCode"
\item "Area" -- inherited from LANDIS-II default table, the Canadian ecoregion
from which traits where derived. Not used during the simulation
\item "firetolerance" -- module\emph{relative} (to other species) fire tolerance
\item "growthcurve" and "mortalityshape" -- growth curve shape parameters.
\item "longevity" -- maximum species age
\item "postfireregen" -- post-fire regeneration strategy ("serotiny", "resprout" or "none")
\item "resproutprob" -- probability of resprouting
\item "resproutage_min" -- minimum age at which species is capable of resprouting
\item "resproutage_max" -- maximum age at which species is capable of resprouting
\item "seeddistance_eff" -- effective dispersal distance
\item "seeddistance_max" -- maximum dispersal distance
\item "shadetolerance" -- \emph{relative} (to other species) shade tolerance
\item "sexualmature" -- age at sexual maturity
Known optional parameters added/needed by some modules (the user may add others for their own modules):
\item "inflationFactor" -- \code{Biomass_speciesParameters} module: inflation factor for \code{maxB}
\item "growthCurveSource" -- \code{Biomass_speciesParameters} module: how "growthcurve" was estimated
\item "mANPPproportion" -- \code{Biomass_speciesParameters} module: multiplication factor to calculate \code{maxANPP} from \code{maxB}
\item "thermokarsttol" -- \code{Biomass_disturbances} module: proportion of biomass surviving after thermokarst (i.e. permafrost thaw). Applied equally across cohorts.
Parameters inherited from LANDIS-II default table, but not used in LandR at the moment:
\item "leaflongevity"
\item "wooddecayrate"
\item "leafLignin"
\item "hardsoft"
Please see the LANDIS-II Biomass Succession Extension v3.2.1 manual (Scheller and Miranda 2015)
for further detail.
}}

\item{sufficientLight}{a \code{data.frame} (not \code{data.table}!) containing probability
of establishment, given a site's light conditions (\code{X0}-\code{X5}) for each level
of a species shade tolerance (\code{1}-\code{5}).}

\item{speciesEcoregion}{A \code{data.table} with \code{species}-\code{ecoregion}-specific species trait values.
Ecoregion refers to "ecolocation", a categorical variable grouping sites with similar biophysical
characteristics. The table should have at least the following columns: \code{speciesCode} and
\code{ecoregionGroup}, character representation of species and ecoregion groups respectively,
\code{maxB} the maximum biomass for the species in a given 'ecoregion', \code{maxANPP} the maximum
aboveground net primary productivity and \code{SEP} the species establishment probability.
May contain columns \code{inflationFactor}  (used to adjust \code{maxB}) and \code{mANPPproportion}
(used to calculate \code{maxANPP}).}

\item{currentTime}{integer. The current simulation time obtained with \code{time(sim)}}

\item{treedFirePixelTableSinceLastDisp}{a vector of pixels that burnt and were forested
in the previous time step.}

\item{calibrate}{logical. Determines whether to output \code{postFirePixelCohortData}.
Defaults to FALSE}
}
\value{
A list of objects:
\code{postFirePixelCohortData}, a \code{data.table} with the cohorts that undergo serotiny;
\code{serotinyPixel}, a vector of pixels where serotiny was activated;
\code{postFireRegenSummary}, the updated \code{postFireRegenSummary}, if \code{calibrate = TRUE}.
}
\description{
Activate serotiny after a (fire) disturbance
}
