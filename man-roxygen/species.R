#' @param species a `data.table` with species traits such as longevity, shade tolerance, etc.
#'   Must have column `speciesCode`, with species names/IDs. The following is a
#'   list of default trait columns:
#'   - "species" same as "speciesCode" -- species ID name
#'   - "speciesCode"
#'   - "Area" -- inherited from LANDIS-II default table, the Canadian ecoregion
#'     from which traits where derived. Not used during the simulation
#'   - "firetolerance" -- module*relative* (to other species) fire tolerance
#'   - "growthcurve" and "mortalityshape" -- growth curve shape parameters.
#'   - "longevity" -- maximum species age
#'   - "postfireregen" -- post-fire regeneration strategy ("serotiny", "resprout" or "none")
#'   - "resproutprob" -- probability of resprouting
#'   - "resproutage_min" -- minimum age at which species is capable of resprouting
#'   - "resproutage_max" -- maximum age at which species is capable of resprouting
#'   - "seeddistance_eff" -- effective dispersal distance
#'   - "seeddistance_max" -- maximum dispersal distance
#'   - "shadetolerance" -- *relative* (to other species) shade tolerance
#'   - "sexualmature" -- age at sexual maturity
#'   Known optional parameters added/needed by some modules (the user may add others for their own modules):
#'   - "inflationFactor" -- `Biomass_speciesParameters` module: inflation factor for `maxB`
#'   - "growthCurveSource" -- `Biomass_speciesParameters` module: how "growthcurve" was estimated
#'   - "mANPPproportion" -- `Biomass_speciesParameters` module: multiplication factor to calculate `maxANPP` from `maxB`
#'   - "thermokarsttol" -- `Biomass_disturbances` module: proportion of biomass surviving after thermokarst (i.e. permafrost thaw). Applied equally across cohorts.
#'   Parameters inherited from LANDIS-II default table, but not used in LandR at the moment:
#'   - "leaflongevity"
#'   - "wooddecayrate"
#'   - "leafLignin"
#'   - "hardsoft"
#'  Please see the LANDIS-II Biomass Succession Extension v3.2.1 manual (Scheller and Miranda 2015)
#'  for further detail.
