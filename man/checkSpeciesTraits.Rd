% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepInputObjects.R
\name{checkSpeciesTraits}
\alias{checkSpeciesTraits}
\title{Check if all species in have trait values}
\usage{
checkSpeciesTraits(speciesLayers, species, sppColorVect)
}
\arguments{
\item{speciesLayers}{A \code{RasterStack} or \code{RasterLayer} that
should contain species cover data in the study area}

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

\item{sppColorVect}{A named vector of colours to use for plotting.
The names must conform with species name convention used
(see \link{sppEquivalencies_CA} for possible naming
conventions) and should also contain a colour for 'Mixed',
when a mixed forest type is supported  (see \code{vegLeadingProportion}
parameter in \code{\link[=vegTypeGenerator]{vegTypeGenerator()}} for details on mixed
forest types).}
}
\value{
A \code{list} with the \code{speciesLayers} and \code{sppColorVect}
containing only the species that have trait values in \code{species}
}
\description{
Check if all species in have trait values
}
