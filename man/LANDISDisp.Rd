% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seedDispersalLANDIS.R
\name{LANDISDisp}
\alias{LANDISDisp}
\title{Simulate a LANDIS-II dispersal process on a landscape.}
\usage{
LANDISDisp(
  dtSrc,
  dtRcv,
  pixelGroupMap,
  speciesTable,
  dispersalFn = Ward,
  b = 0.01,
  k = 0.95,
  plot.it = FALSE,
  successionTimestep,
  verbose = getOption("LandR.verbose", TRUE),
  ...
)
}
\arguments{
\item{dtSrc}{data.table}

\item{dtRcv}{data.table}

\item{pixelGroupMap}{A \code{RasterLayer} with pixel values equal to a pixel group
number that corresponds exactly to \code{pixelGroup} column in \code{cohortData}.}

\item{speciesTable}{A species traits table, with \strong{at least} the following columns:
\itemize{
\item \code{speciesCode} an character representation of species;
\item \code{Area} a character of geographical area within a species range for which trait values are relevant;
\item \code{longevity} species longevity in years, \code{sexualmature} age in years at sexual maturity;
\item \code{shadetolerance} numeric value between 1-5 of relative shade tolerance (with respect to other species);
\item \code{seeddistance_eff} the "effective" seed dispersal distance;
\item \code{seeddistance_max} a numeric with the maximum seed dispersal distance;
\item \code{mortalityshape} and \code{growthcurve}: growth curve shape parameters.
Other columns (e.g. fire-related traits) can also be included depending on LandR modules in use.
Please see the LANDIS-II Biomass Succession Extension v3.2.1 manual (Scheller and Miranda 2015)
for further detail.
}}

\item{dispersalFn}{An expression that can take a "dis" argument. See details.
Default is "Ward" (temporarily unused, as it is hard coded inside Rcpp function)}

\item{b}{LANDIS Ward seed dispersal calibration coefficient (set to 0.01 in LANDIS)}

\item{k}{LANDIS Ward seed dispersal the probability that seed will disperse within
the effective distance (e.g., 0.95)}

\item{plot.it}{Deprecated. If TRUE, then plot the raster at every interaction,
so one can watch the \code{LANDISDisp} event grow.}

\item{successionTimestep}{The time between successive seed dispersal events.
In LANDIS-II, this is called "Succession Timestep".}

\item{verbose}{Controls message output. Defaults to \code{getOption("LandR.verbose")}}

\item{...}{Additional parameters. Currently none.}
}
\value{
A numeric vector of raster pixel indices, in the same resolution and extent as
\code{seedSrc} raster.
}
\description{
Simulate seed dispersal using user defined function. This is a "receiving pixel" focused
dispersal approach.
It is the "potentially receiving" cell that looks around itself for potential seed sources.
If it finds a single seed source, that passes the probability function described by the
\code{dispersalFn}.
If this passes a comparison to a uniform random draw, then the receiving cell is deemed to have
a "successful" dispersal for that species.
This function can therefore only be used for a relatively specific situation
where there is a yes/no returned for each potential receiving cell, i.e., not abundance.
This function is also not cumulative, i.e,. there is no higher abundance of seeds received if
a receiving cell has lots of seed sources around it vs. a single seed source.
The difference will come with a higher probability of successfully receiving a "seed".
}
\details{
\code{dispersalFn} (temporarily unused as code is converted to Rcpp -- the
default \code{dispersalFn} is hard coded within the \code{spiralSeedDispersal}
function that uses C++) must be an expression that returns a probability
distribution. Because it is a dispersal kernel, it must be a probability
distribution. The expression that can take an argument named "dis" (without
quotes) as this will be calculated internally and represents the distance
from the initial (receiving) pixel and all active pixels within that cluster
of active pixels. \code{SpaDES} includes the \code{\link[=Ward]{Ward()}} kernel as
defined in the LANDIS-II documentation.
}
\examples{
if (require("googledrive")) {
  seed <- sample(1e6, 1)
  set.seed(seed)
  library(data.table)

  # keep this here for interactive testing with a larger raster
  rasterTemplate <- LandR:::rasterRead(terra::ext(0, 2500, 0, 2500), res = 100)

  # make a pixelGroupMap
  pgs <- 4 # make even just because of approach below requires even
  pixelGroupMap <- SpaDES.tools::randomPolygons(rasterTemplate, numTypes = pgs)
  pixelGroupMap[1:100] <- NA # emulate a mask at the start

  # Make a receive pixels table -- need pixelGroup and species
  nSpecies <- 3
  maxNSpeciesPerPixel <- min(5, nSpecies)
  rcvSpByPG <- lapply(seq_len(pgs / 2), function(pg) {
    data.table(speciesCode = sample(nSpecies, size = sample(maxNSpeciesPerPixel, 1)))
  })
  seedReceive <- rbindlist(rcvSpByPG, idcol = "pixelGroup")

  # Make a source pixels table -- need pixelGroup and species
  srcSpByPG <- lapply(seq_len(pgs / 2), function(pg) {
    data.table(speciesCode = sample(nSpecies, size = sample(maxNSpeciesPerPixel, 1)))
  })
  seedSource <- rbindlist(srcSpByPG, idcol = "pixelGroup")
  # make source pixels not same pixelGroups as receive
  seedSource[, pixelGroup := pixelGroup + pgs / 2]

  # Get a species table -- if using in Canada, can use this
  speciesTable <- getSpeciesTable(dPath = tempdir())
  speciesTable <- speciesTable[Area == "BSW"]
  speciesTable[, speciesCode := as.factor(LandisCode)]
  speciesTable[, seeddistance_eff := SeedEffDist]
  speciesTable[, seeddistance_max := SeedMaxDist]

  speciesTable <- speciesTable
  speciesTable <- data.table(speciesTable)[, speciesCode := seq_along(LandisCode)]
  seedReceiveFull <- speciesTable[seedReceive, on = "speciesCode"]
  output <- LANDISDisp(
    dtRcv = seedReceiveFull, plot.it = interactive(),
    dtSrc = seedSource,
    speciesTable = speciesTable,
    pixelGroupMap,
    verbose = TRUE,
    successionTimestep = 10
  )
  # Summarize
  output[, .N, by = speciesCode]

  ## Plot the maps
  if (interactive()) {
    library(quickPlot)
    clearPlot()
    spMap <- list()
    spMap$pixelGroupMap <- pixelGroupMap
    for (sppp in unique(output$speciesCode)) {
      spppChar <- paste0("Sp_", sppp)
      spMap[[spppChar]] <- LandR:::rasterRead(pixelGroupMap)
      ss <- unique(seedSource[speciesCode == sppp], on = c("pixelGroup", "speciesCode"))
      spMap[[spppChar]][pixelGroupMap[] \%in\% ss$pixelGroup] <- 1

      receivable <- LandR:::rasterRead(pixelGroupMap)
      srf <- unique(seedReceiveFull[speciesCode == sppp], on = c("pixelGroup", "speciesCode"))
      receivable[pixelGroupMap[] \%in\% srf$pixelGroup] <- 1

      forest <- which(!is.na(pixelGroupMap[]))
      src <- which(!is.na(spMap[[spppChar]][]))
      recvable <- which(!is.na(receivable[]))
      rcvd <- output[speciesCode == sppp]$pixelIndex

      spMap[[spppChar]][forest] <- 0
      spMap[[spppChar]][recvable] <- 2
      spMap[[spppChar]][src] <- 1
      spMap[[spppChar]][rcvd] <- 3
      spMap[[spppChar]][intersect(src, rcvd)] <- 4

      levels(spMap[[spppChar]]) <- data.frame(ID = 0:4,
                                              type = c("OtherForest", "Source", "Didn't receive",
                                                       "Received", "Src&Rcvd"))
    }
    Plot(spMap, cols = "Set2")

    # A summary
    rr <- apply(rast(spMap)[[-1]][] + 1, 2, tabulate)
    rownames(rr) <- raster::levels(spMap[[2]])[[1]][,"type"][1:NROW(rr)]
    # next line only works if there are some places that are both source and potential to receive
    # rr <- rbind(rr, propSrcRcved = round(rr[5,]/ (rr[5,]+rr[2,]), 2))
 }
}

}
\author{
Eliot McIntire
}
