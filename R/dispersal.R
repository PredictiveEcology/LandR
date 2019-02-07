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
#' @param advectionDir A single number in degrees from North = 0, indicating
#'   the direction of advective forcing (i.e., wind). Will soon allow a
#'   raster of advection vectors.
#' @param advectionMag A single number in distance units of the \code{rasQuality},
#'   e.g., meters, indicating the relative forcing that will occur for
#'   each 1 pixel of attempted spreading.
#' @param meanDist A single number indicating the mean distance parameter
#'    for a negative exponential distribution dispersal kernel (e.g., \code{dexp}).
#' @param plot.it Logical. If \code{TRUE}, there will be 2 plots that occur
#'    as iterations happen.
#' @param minNumAgents Single numeric indicating the minimum number of agents
#'    to consider all dispersing finished. Default is 50
#' @examples
#'
#' library(raster)
#' library(sf)
#' library(quickPlot)
#' #a <- randomStudyArea(size = 1e8)
#' maxDim <- 10000
#' ras <- raster::raster(extent(c(0, maxDim, 0, maxDim)), res = 100, vals = 0)
#' rasQuality <- raster(ras)
#' rasQuality[] <- 1
#' #ras <- raster::raster(extent(a), res = 100)
#' #mask <- fasterize::fasterize(st_as_sf(a), ras)
#' #rasQuality <- gaussMap(ras)
#' #crs(rasQuality) <- crs(a)
#' #rasQuality[is.na(mask)] <- NA
#' # rescale so min is 0.75 and max is 1
#' #rasQuality[] <- rasQuality[] / (maxValue(rasQuality * 4) ) + 3/4
#' #rasQuality[] <- 1
#' rasAbundance <- raster(rasQuality)
#' rasAbundance[] <- 0
#' startPixel <- middlePixel(rasAbundance)
#' startPixel <- sample(seq(ncell(rasAbundance)), 3)
#' rasAbundance[startPixel] <- 1000
#' advectionDir <- 90
#' advectionMag <- 4 * res(rasAbundance)[1]
#'
#' clearPlot()
#' spread3(rasAbundance = rasAbundance,
#'         rasQuality = rasQuality,
#'         advectionDir = advectionDir,
#'         advectionMag = advectionMag,
#'         meanDist = 600)
#'
#'
#' @importFrom CircStats deg rad
#' @importFrom fpCompare %>=%
#' @importFrom raster xyFromCell
#' @importFrom SpaDES.tools spread2 pointDistance
#' @importFrom quickPlot Plot clearPlot
#' @importFrom data.table := setattr [
spread3 <- function(start, rasQuality, rasAbundance, advectionDir,
                    advectionMag, kernel, meanDist, plot.it = TRUE,
                    minNumAgents = 50) {
  if (advectionDir > 2 * pi) {
    message("assuming that advectionDir is in geographic degrees")
    advectionDir <- CircStats::rad(advectionDir)
  }

  if (missing(start))
    start <- which(!is.na(rasAbundance[]) & rasAbundance[] > 0)

  start <- spread2(rasQuality, start, iterations = 0, returnDistances = TRUE,
          returnDirections = TRUE, returnFrom = TRUE, asRaster = FALSE)
  start[, `:=`(abundActive = rasAbundance[][start$pixels],
               indWithin = 1L,
               indFull = 1L)]
  abundanceDispersing <- sum(start$abundActive)
  rasIterations <- raster(rasQuality)
  rasIterations[] <- NA
  rasIterations[start$pixels] <- 0
  while (abundanceDispersing > minNumAgents) {
    b <- spread2(landscape = rasQuality, start = start,
                 spreadProb = 1, iterations = 1, asRaster = FALSE,
                 returnDistance = TRUE, returnFrom = TRUE,
                 returnDirections = TRUE,
                 circle = TRUE, allowOverlap = 3)
    #b <- b[!duplicated(b, by = c("initialPixels", "pixels"))]
    spreadState <- attr(b, "spreadState")
    # faster than assessing with a which()
    active <- spreadState$whActive
    inactive <- spreadState$whInactive


    iteration <- spreadState$totalIterations
    if (isTRUE(plot.it)) {
      rasIterations[b[active]$pixels] <- iteration
      Plot(rasIterations, new = TRUE)
    }

    fromPts <- xyFromCell(rasQuality, b[active]$from)
    toPts <- xyFromCell(rasQuality, b[active]$pixels)
    dists <- pointDistance(p1 = fromPts, p2 = toPts, lonlat = FALSE)
    dirs <- b[active]$direction
    #dists <- b[active]$distance
    xDist <- round(sin(advectionDir) * advectionMag + sin(dirs) * dists, 4)
    yDist <- round(cos(advectionDir) * advectionMag + cos(dirs) * dists, 4)

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
                   #maxAbund = max(abund, na.rm = TRUE),
                   indWithin = seq(.N),
                   indFull = .I), by = c('initialPixels', 'pixels')]
    b[active, meanNumNeighs := mean(lenSrc / lenRec) * mean(mags), by = c("pixels", "initialPixels")]
    b[active, sumAbund2 := sumAbund * meanNumNeighs/ mean(mags)]
    toSubtract <- sum(b[active]$sumAbund2 - b[active]$sumAbund, na.rm = TRUE)
    totalSumAbund <- sum(b[active]$sumAbund, na.rm = TRUE)
    totalSumAbund2 <- sum(b[active]$sumAbund2, na.rm = TRUE)
    multiplyAll <- totalSumAbund/totalSumAbund2
    b[active, sumAbund := sumAbund2 * multiplyAll]
    b[, abund := NULL]
    #b[is.na(ind), ind := 1]
    #rmFromActive <- which(b$ind > 1)
    keepRows <- which(b$indWithin == 1)
    b2 <- copy(b)
    b <- b[keepRows]
    active <- na.omit(match(active, b$indFull))
    b[, indFull := seq(NROW(b))]
    #active <- active[na.omit(match(keepRows, active))]

    # Some of those active will not stop: estimate here by kernel probability
    b[active, abundSettled := pexp(q = distance, rate = 1/meanDist) * sumAbund]
    # Some of the estimated dropped will not drop because of quality
    b[active, abundSettled := floor(abundSettled  * rasQuality[][pixels])]
    b[active, abundActive := sumAbund - abundSettled]

    abundanceDispersing <- sum(b[active]$abundActive, na.rm = TRUE)
    if (isTRUE(plot.it)) {
      b1 <- copy(b)
      b2 <- b1[, sum(abundSettled), by = "pixels"]
      rasAbundance[b2$pixels] <- ceiling(b2$V1)
      Plot(rasAbundance, new = TRUE)
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

  return(start)
}

#' Return the (approximate) middle pixel on a raster
#'
#' This calculation is different depending on whether
#' the \code{nrow(ras)} and \code{ncol(ras)} are even or odd.
#' @export
middlePixel <- function(ras) {
  if (nrow(ras) %% 2 == 1) {
    floor(ncell(ras) / 2)
  } else {
    floor(nrow(ras)/2) * ncol(ras) - floor(ncol(ras)/2)
  }
}
