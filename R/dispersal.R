if (getRversion() >= "3.1.0") {
  utils::globalVariables(
    c("abund", "abundActive", "abundSettled", "direction", "distance", "distGrp",
      "from", "indFull", "lenRec", "lenSrc", "mags", "meanNumNeighs",
      "newDirs", "newMags", "prop", "srcAbundActive", "sumAbund", "sumAbund2")
  )
}

#' An alternative spread function -- conceived for insects
#'
#' This is built with \code{\link[SpaDES.tools]{spread2}} and
#' is still experimental. This one differs from other attempts
#' in that it treats the advection and dispersal as mathematical
#' vectors that are added together. They are "rounded" to pixel
#' centres.
#'
#' @param start Raster indices from which to initiate dispersal
#' @param rasQuality A raster with habitat quality. Currently, must
#'   be scaled from 0 to 1, i.e., a probability of "settling"
#' @param rasAbundance A raster where each pixel represents the number
#'   of "agents" or pseudo-agents contained. This number of agents, will
#'   be spread horizontally, and distributed from each pixel
#'   that contains a non-zero non NA value.
#' @param advectionDir A single number or \code{RasterLayer} in degrees
#'   from North = 0 (though it will use radians if all values are
#'   \code{abs(advectionDir) > 2 * pi)}. This indicates
#'   the direction of advective forcing (i.e., wind).
#' @param advectionMag A single number or \code{RasterLayer} in distance units of the
#'   \code{rasQuality}, e.g., meters, indicating the relative forcing that will
#'   occur. It is imposed on the total event, i.e., if the \code{meanDist} is
#'   \code{10000}, and \code{advectionMag} is \code{5000}, then the expected
#'   distance (i.e., 63\% of agents) will have settled by \code{15000} map units.
#' @param meanDist A single number indicating the mean distance parameter in map units
#'    (not pixels), for a negative exponential distribution
#'    dispersal kernel (e.g., \code{dexp}). This will mean that 63% of agents will have
#'    settled at this \code{meanDist} (still experimental)
#' @param verbose Numeric. With increasing numbers above 0, there will be more
#'     messages produced. Currently, only 0, 1, or 2+ are discinct.
#' @param plot.it Numeric. With increasing numbers above 0, there will be plots
#'     produced during iterations. Currently, only 0, 1, or 2+ are discinct.
#' @param minNumAgents Single numeric indicating the minimum number of agents
#'    to consider all dispersing finished. Default is 50
#' @param saveStack If provided as a character string, it will save each iteration
#'   as part of a \code{rasterStack} to disk upon exit.
#'
#' @return
#' A \code{data.table} with all information used during the spreading
#'
#' @examples
#'
#' #########################################################
#' # Simple case, no variation in rasQuality, numeric advectionDir and advectionMag
#' #########################################################
#' library(raster)
#' library(quickPlot)
#' maxDim <- 10000
#' ras <- raster::raster(extent(c(0, maxDim, 0, maxDim)), res = 100, vals = 0)
#' rasQuality <- raster(ras)
#' rasQuality[] <- 1
#' rasAbundance <- raster(rasQuality)
#' rasAbundance[] <- 0
#' # startPixel <- middlePixel(rasAbundance)
#' startPixel <- sample(seq(ncell(rasAbundance)), 30)
#' rasAbundance[startPixel] <- 1000
#' advectionDir <- 70
#' advectionMag <- 4 * res(rasAbundance)[1]
#' meanDist <- 2600
#'
#' dev() # don't use Rstudio windows, which is very slow
#' clearPlot()
#' out <- spread3(rasAbundance = rasAbundance,
#'         rasQuality = rasQuality,
#'         advectionDir = advectionDir,
#'         advectionMag = advectionMag,
#'         meanDist = meanDist, verbose = 2,
#'         plot.it = 2)
#'
#' # Test the dispersal kernel
#' out[, disGroup:=round(distance/100)*100]
#' freqs <- out[, .N, by = "disGroup"]
#' freqs[, cumSum := cumsum(N) ]
#' #plot(freqs$disGroup, freqs$cumSum) # can plot the distance X number
#' #abline(v = advectionMag + meanDist)
#' # should be 0.63:
#' freqs[disGroup == advectionMag + meanDist, cumSum]/
#'    tail(freqs,1)[, cumSum]
#'
#' #########################################################
#' ### The case of variable quality raster
#' #########################################################
#' library(sf) # needed to use fasterize
#' library(SpaDES.tools) # for gaussMap
#' library(fasterize) # faster than raster::rasterize
#' a <- randomStudyArea(size = 1e9)
#' ras <- raster(extent(a), res = 100)
#' mask <- fasterize(st_as_sf(a), ras)
#' rasQuality <- gaussMap(ras)
#' crs(rasQuality) <- crs(a)
#' rasQuality[is.na(mask)] <- NA
#' # rescale so min is 0.75 and max is 1
#' rasQuality[] <- rasQuality[] / (maxValue(rasQuality) * 4 ) + 3/4
#' rasAbundance <- raster(rasQuality)
#' rasAbundance[] <- 0
#' startPixel <- sample(seq(ncell(rasAbundance)), 300)
#' rasAbundance[startPixel] <- 1000
#' advectionDir <- 90
#' advectionMag <- 4 * res(rasAbundance)[1]
#' clearPlot()
#' out <- spread3(rasAbundance = rasAbundance,
#'         rasQuality = rasQuality,
#'         advectionDir = advectionDir,
#'         advectionMag = advectionMag,
#'         meanDist = 600, verbose = 2,
#'         plot.it = 1)
#'
#' ###############################################################################
#' ### The case of variable quality raster, raster for advectionDir & advectionMag
#' ###############################################################################
#' library(raster)
#' library(quickPlot)
#' library(SpaDES.tools)
#' maxDim <- 10000
#' ras <- raster::raster(extent(c(0, maxDim, 0, maxDim)), res = 100, vals = 0)
#' rasQuality <- raster(ras)
#' rasQuality[] <- 1
#' rasAbundance <- raster(rasQuality)
#' rasAbundance[] <- NA
#' # startPixel <- middlePixel(rasAbundance)
#' startPixel <- sample(seq(ncell(rasAbundance)), 25)
#' rasAbundance[startPixel] <- 1000
#'
#' # raster for advectionDir
#' advectionDir <- gaussMap(ras)
#' crs(advectionDir) <- crs(rasQuality)
#' # rescale so min is 0.75 and max is 1
#' advectionDir[] <- advectionDir[] / (maxValue(advectionDir)) * 180
#'
#' # raster for advectionMag
#' advectionMag <- gaussMap(ras)
#' crs(advectionMag) <- crs(rasQuality)
#' # rescale so min is 0.75 and max is 1
#' advectionMag[] <- advectionMag[] / (maxValue(advectionMag)) * 600
#'
#' dev() # don't use Rstudio windows, which is very slow
#' clearPlot()
#' Plot(advectionDir, title = "Wind direction", cols = "Reds")
#' Plot(advectionMag, title = "Wind speed", cols = "Blues")
#' out <- spread3(rasAbundance = rasAbundance,
#'         rasQuality = rasQuality,
#'         advectionDir = advectionDir,
#'         advectionMag = advectionMag,
#'         meanDist = 2600, verbose = 2,
#'         plot.it = 1)
#' Plot(rasAbundance, addTo = "rasAbundance", cols = "black", title = "")
#'
#' #########################################
#' # save iterations to a stack to make animated GIF
#' ########################################
#' tmpStack <- raster::rasterTmpFile("stackToAnimate")
#' out <- spread3(rasAbundance = rasAbundance,
#'         rasQuality = rasQuality,
#'         advectionDir = advectionDir,
#'         advectionMag = advectionMag,
#'         meanDist = 2600, verbose = 2,
#'         plot.it = 0, saveStack = tmpStack)
#'
#'  out2 <- raster::stack(tmpStack)
#'  gifName <- file.path(tempdir(), "animation.gif")
#'  animation::saveGIF(interval = 0.1,
#'                     movie.name = gifName,
#'                     expr = {
#'                       for(i in seq(numLayers(out2)))
#'                         plot(out2[[i]])
#'                      })
#'
#'
#' @importFrom CircStats deg rad
#' @importFrom fpCompare %>=% %>>%
#' @importFrom raster xyFromCell
#' @importFrom SpaDES.tools spread2
#' @importFrom raster pointDistance
#' @importFrom quickPlot Plot clearPlot
#' @importFrom data.table := setattr
#' @importFrom stats pexp
#' @export
spread3 <- function(start, rasQuality, rasAbundance, advectionDir,
                    advectionMag, meanDist, plot.it = 2,
                    minNumAgents = 50, verbose = getOption("LandR.verbose", 0),
                    saveStack = NULL) {

  testEquivalentMetadata(rasAbundance, rasQuality)

  if (is(advectionDir, "Raster")) {
    testEquivalentMetadata(rasAbundance, advectionDir)
    advectionDir <- advectionDir[]
  } else if (length(advectionDir) != 1) {
    if (length(advectionDir) != ncell(rasAbundance)) {
      stop("advectionDir must be length 1, length ncell(rasAbundance), or a Raster with ",
           "identical metadata as rasAbundance")
    }
  }
  if (is(advectionMag, "Raster")) {
    testEquivalentMetadata(rasAbundance, advectionMag)
    advectionMag <- advectionMag[]
  } else if (length(advectionMag) != 1) {
    if (length(advectionMag) != ncell(rasAbundance)) {
      stop("advectionMag must be length 1, length ncell(rasAbundance), or a Raster with ",
           "identical metadata as rasAbundance")
    }
  }

  if (any(abs(advectionDir) > 2 * pi)) {
    messAngles <- "degrees"
    advectionDir <- CircStats::rad(advectionDir)
  } else {
    messAngles <- "radians"
  }
  message("assuming that advectionDir is in geographic ", messAngles,
          "(i.e., North is 0)")

  if (missing(start))
    start <- which(!is.na(rasAbundance[]) & rasAbundance[] > 0)

  start <- spread2(rasQuality, start, iterations = 0, returnDistances = TRUE,
          returnDirections = TRUE, returnFrom = TRUE, asRaster = FALSE)
  start[, `:=`(abundActive = rasAbundance[][start$pixels],
               abundSettled = 0)]
  abundanceDispersing <- sum(start$abundActive)
  plotMultiplier <- mean(start$abundActive) /
    ((meanDist * 10 / res(rasQuality)[1]))
  rasIterations <- raster(rasQuality)
  rasIterations[] <- NA
  rasIterations[start$pixels] <- 0

  while (abundanceDispersing > minNumAgents) {
    b <- spread2(landscape = rasQuality, start = start,
                 spreadProb = 1, iterations = 1, asRaster = FALSE,
                 returnDistances = TRUE, returnFrom = TRUE,
                 returnDirections = TRUE,
                 circle = TRUE, allowOverlap = 3)
    #b <- b[!duplicated(b, by = c("initialPixels", "pixels"))]
    spreadState <- attr(b, "spreadState")
    # faster than assessing with a which()
    active <- spreadState$whActive
    inactive <- spreadState$whInactive


    iteration <- spreadState$totalIterations
    if (verbose > 1) message("Iteration ", iteration)
    if (isTRUE(plot.it > 1)) {
      rasIterations[b[active]$pixels] <- iteration
      Plot(rasIterations, new = iteration == 1,
           legendRange = c(0, meanDist / (res(rasQuality)[1] / 12)))
    }

    fromPts <- xyFromCell(rasQuality, b[active]$from)
    toPts <- xyFromCell(rasQuality, b[active]$pixels)
    dists <- pointDistance(p1 = fromPts, p2 = toPts, lonlat = FALSE)
    dirs <- b[active]$direction

    # Convert advection vector into length of dirs from pixels, if length is not 1
    advectionDirTmp <- if (length(advectionDir) > 1) {
      advectionDir[b[active]$pixels]
    } else {
      advectionDir
    }
    advectionMagTmp <- if (length(advectionMag) > 1) {
      advectionMag[b[active]$pixels]
    } else {
      advectionMag
    }
    xDist <- round(sin(advectionDirTmp) * advectionMagTmp + sin(dirs) * dists, 4)
    yDist <- round(cos(advectionDirTmp) * advectionMagTmp + cos(dirs) * dists, 4)

    # This calculates: "what fraction of the distance being moved is along the dirs axis"
    #   This means that negative mags is "along same axis, but in the opposite direction"
    #   which is dealt with next, see "opposite direction"
    b[active, mags := round(sin(dirs) * xDist + cos(dirs) * yDist, 3)]
    negs <- b[active]$mags < 0
    negs[is.na(negs)] <- FALSE
    #dirs2 <- dirs
    anyNegs <- any(negs)
    nonNA <- !is.na(b[active]$direction)
    if (any(anyNegs)) { # "opposite direction"
      dirs2 <- (dirs[negs] + pi) %% (2*pi)
      b[active[negs], mags := -mags]
      b[active[negs], newDirs := dirs2]
      b[active[nonNA], newMags := mags + mags[match(round(direction, 4), round(newDirs, 4))],
        by = "initialPixels"]
      nonNANewMags <- !is.na(b[active]$newMags)
      b[active[nonNANewMags], mags := newMags]
      nonNANewDirs <- !is.na(b[active]$newDirs)
      b[active[nonNANewDirs], mags := 0]
      set(b, NULL, c("newDirs", "newMags"), NULL)
    }
    b[active[nonNA], prop := round(mags / sum(mags), 3), by = c("from", "initialPixels")]
    if (FALSE) # almost
      b[active, -c("abundActive")][b[inactive, c("pixels", "initialPixels", "abundActive")],
                                   on = c("from" = "pixels", "initialPixels")]

    b[distance %>>% ((iteration - 2) * res(rasAbundance)[1]),
      srcAbundActive := abundActive[match(from, pixels)], by = "initialPixels"]


    # Expected number, based on advection and standard spread2
    b[active, abund := srcAbundActive * prop]

    b[active, lenRec := .N, by = c("pixels", "initialPixels")]
    b[active, lenSrc := min(2.5, .N), by = c("from", "initialPixels")]

    # Sum all within a receiving pixel,
    #    then collapse so only one row per receiving cell,
    #    it is a markov chain of order 1 only, except for some initial info
    b[active, `:=`(sumAbund = sum(abund, na.rm = TRUE),
                   indWithin = seq(.N),
                   indFull = .I), by = c('initialPixels', 'pixels')]
    b[active, meanNumNeighs := mean(lenSrc / lenRec) * mean(mags), by = c("pixels", "initialPixels")]
    keepRows <- which(b$indWithin == 1 | is.na(b$indWithin))
    b <- b[keepRows]
    active <- na.omit(match(active, b$indFull))
    #b[, indFull := seq(NROW(b))]

    b[active, sumAbund2 := sumAbund * meanNumNeighs/ mean(mags)]

    totalSumAbund <- sum(b[active]$sumAbund, na.rm = TRUE)
    totalSumAbund2 <- sum(b[active]$sumAbund2, na.rm = TRUE)
    multiplyAll <- totalSumAbund/totalSumAbund2

    b[active, sumAbund := sumAbund2 * multiplyAll]
    set(b, NULL, c("abund", "sumAbund2", "mags", "lenSrc", "lenRec",
                   "meanNumNeighs", "prop", "indWithin", "indFull"),
        NULL)

    # Some of those active will not stop: estimate here by kernel probability
    advectionMagTmp <- if (length(advectionMag) > 1) {
      advectionMag[b[active]$pixels]
    } else {
      advectionMag
    }
    b[active, abundSettled :=
        pexp(q = distance, rate = pi/(meanDist+advectionMagTmp)^1.5) * sumAbund] # kernel is 1 dimensional,
    # b[active, abundSettled :=
    #     dexp(x = distance, rate = 1/(meanDist+advectionMag)) * sumAbund] # kernel is 1 dimensional,
    # but spreading is dropping agents in 2 dimensions
             # It doesn't work to use ^2, I think because we are discretizing the landscape
             # from a continuous surface, so, the number of pixels with agents settled
             # is not actually the full square on a 1 dimensional line ... I might be wrong
    # Some of the estimated dropped will not drop because of quality
    #   First place to round to whole numbers
    b[active, abundSettled := abundSettled  * rasQuality[][pixels]]
    b[active, abundActive := sumAbund - abundSettled]
    b[active[active %in% which(b$abundActive < 1)], abundActive := 0]

    abundanceDispersing <- sum(b[active]$abundActive, na.rm = TRUE)
    if (verbose > 1) message("Number still dispersing ", abundanceDispersing)
    if (isTRUE(plot.it > 0)) {
      b2 <- b[, sum(abundSettled), by = "pixels"]
      rasAbundance[b2$pixels] <- ceiling(b2$V1)
      needNew <- FALSE
      if (max(ceiling(b2$V1), na.rm = TRUE) > plotMultiplier) {
        plotMultiplier <- plotMultiplier * 1.5
        needNew <- TRUE
      }
      Plot(rasAbundance, new = iteration == 1 || needNew,
           legendRange = c(0, plotMultiplier), title = "Abundance")
    }

    newInactive <- b[active]$abundActive == 0
    b[active[newInactive], state := "inactive"]
    spreadState$whActive <- active[!newInactive]
    spreadState$whInactive <- c(spreadState$whInactive, active[newInactive])
    setattr(b, "spreadState", spreadState)

    # clean up temporary columns
    set(b, NULL, c("sumAbund", "srcAbundActive"), NULL)

    start <- b
  }
  if (!is.null(saveStack)) {
    b[, distGrp := floor(distance/100)]
    ras <- raster(rasAbundance)
    out1 <- lapply(unique(b$distGrp), function(x)  {
      r <- raster(ras)
      x1 <- b[distGrp <= x, sum(abundSettled), by = "pixels"]
      r[x1$pixels] <- ceiling(x1$V1)
      r
    })
    writeRaster(raster::stack(out1), filename = saveStack, overwrite = TRUE)
    message("stack saved to ", saveStack)
  }

  return(start)
}

#' Return the (approximate) middle pixel on a raster
#'
#' This calculation is slightly different depending on whether
#' the \code{nrow(ras)} and \code{ncol(ras)} are even or odd.
#' It will return the exact middle pixel if these are odd, and
#' the pixel just left and/or above the middle pixel if either
#' dimension is even, respectively.
#' @param ras A \code{Raster}
#' @export
middlePixel <- function(ras) {
  if (nrow(ras) %% 2 == 1) {
    floor(ncell(ras) / 2)
  } else {
    floor(nrow(ras)/2) * ncol(ras) - floor(ncol(ras)/2)
  }
}

#' Test that metadata of 2 or more objects is the same
#'
#' Currently, only Raster class has a useful method. Defaults to
#' \code{all(sapply(list(...)[-1], function(x) identical(list(...)[1], x)))}
#' @param ... 2 or more of the same type of object to test
#'   for equivalent metadata
#' @export
testEquivalentMetadata <- function(...) {
  UseMethod("testEquivalentMetadata")
}


#' @export
testEquivalentMetadata.Raster <- function(...) {
  d <- list(...)
  res <- lapply(d[-1], function(x) {
    if (!identical(d[[1]]@extent, x@extent)) {
      stop("rasAbundance and rasQuality must have same extent")
    }
    if (!identical(d[[1]]@crs, x@crs)) {
      stop("rasAbundance and rasQuality must have same crs")
    }
    if (!identical(res(d[[1]]), res(x))) {
      stop("rasAbundance and rasQuality must have same res")
    }

  })
  return(invisible())
}
