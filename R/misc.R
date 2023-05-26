utils::globalVariables(c(
  ".", ":=", "highProb", "highShadetol", "lightProb", "lowProb", "lowShadetol",
  "shadetolerance", "siteShade", "year"
))

.extractFunction <- getFromNamespace(".extractFunction", "reproducible")

#' Assign light probability
#'
#' @template sufficientLight
#'
#' @param newCohortData  a modified version of `cohortData` that contains new cohorts.
#'
#' @param interpolate Logical. Activates interpolation of probabilities of establishment between
#'   any two values of shade tolerance in the sufficient light table, allowing species shade tolerance
#'   trait values to take any decimal value between 1 and 5 (inclusively). If false, species shade tolerances
#'   can only take integer values between 1 and 5  (inclusively).
#' @template doAssertion
#'
#' @return  `newCohortData` with a `lightProb` column
#'
#' @export
assignLightProb <- function(sufficientLight, newCohortData, interpolate = TRUE,
                            doAssertion = getOption("LandR.assertions", TRUE)) {

  ## 2022-06-22 AMC removed the assertion and always coerce to data.frame as needed
  if (!is(sufficientLight, "data.frame") | is(sufficientLight, "data.table")) {
    warning("Coercing object 'sufficientLight' from ", is(sufficientLight)[1], " to data.frame.")
    sufficientLight <- as.data.frame(sufficientLight)
    stopifnot(is(sufficientLight, "data.frame"))
  }

  ## for each line, get the survival probability from sufficientLight table note
  ## that sufficentLight is a table of survival probs for each tolerance level
  ## (row) by and shade level (column) siteShade + 2 is necessary to skip the
  ## first column
  if (interpolate) {
    ## calculate floors/ceilings of shade tolerance and corresponding probs.
    tempDT <- newCohortData[, list(shadetolerance = shadetolerance,
                                   siteShade = siteShade,
                                   lowShadetol = floor(shadetolerance),
                                   highShadetol = ceiling(shadetolerance))]
    tempDT[, `:=`(lowProb = sufficientLight[cbind(lowShadetol, siteShade + 2)],
                  highProb = sufficientLight[cbind(highShadetol, siteShade + 2)])]

    ## interpolate between floor/ceiling to find prob, then add to newCohortData
    tempDT[, lightProb := .interpolateLightProb(x = shadetolerance, x0 = lowShadetol, x1 = highShadetol,
                                                y0 = lowProb, y1 = highProb)]
    newCohortData[, lightProb := tempDT$lightProb]
  } else {
    ## are there any decimals in shade tolerance trait values?
    if (!all(newCohortData$shadetolerance == round(newCohortData$shadetolerance)))
      stop(paste("Species shade tolerance values (in sim$species) have decimals,",
                 "but interpolation of germination probabilities between",
                 "shade tolerance categories in sim$sufficientLight is FALSE.\n",
                 "Set 'interpolate = TRUE', or provide integer shadetolerance values."))

    newCohortData[, lightProb := sufficientLight[cbind(shadetolerance, siteShade + 2)]]
  }
}

#' Find interpolated value of light probability
#'
#' @param x the species shade tolerance trait value for which we want to find
#'   the interpolated probability.
#'
#' @param x0 the `floor` of `x` corresponding to a class of shade
#'   tolerance in the `sufficientLight` table.
#'
#' @param x1 the `ceiling` of `x` corresponding to a class of shade
#'   tolerance in the `sufficientLight` table.
#'
#' @param y0 the probability of germination in the `sufficientLight` table
#'   corresponding to `x0`.
#'
#' @param y1 the probability of germination in the `sufficientLight` table
#'   corresponding to `x1`.
#'
#' @return  `vector` of the interpolated value

.interpolateLightProb <- function(x, x0, x1, y0, y1) {
  ## if floor and ceiling are the same and equal to x, shade tolerances (x) are integers
  ## and the probability does not need to be interpolated
  if (all((x1 - x0) == 0) & all(x == x1)) {
    y <- y0
  } else {
    y <- ((y1 - y0)/(x1 - x0)) * (x - x0) + y0
    y[is.nan(y)] <- y0[is.nan(y)]  ## as before, when shade tol. is an integer
  }
  y
}

#' Convert numeric values to rounded integers
#'
#' Essentially a wrapper around `round`, rather than `truncate`, which is what `as.integer`
#' does. Internally, this is simply `as.integer(floor(x + 0.5))`.
#'
#' @note Values ending in `.5` will be rounded up, whether positive or negative.
#' This is different than `round`.
#'
#' @param x A numeric vector
#'
#' @return An integer vector of length `x`, rounded to zero decimal places
#'   prior to `as.integer`
#'
#' @export
#' @examples
#' x <- seq(-2, 2, 0.25)
#' data.frame(dbl = x, int = asInteger(x))
asInteger <- function(x)
  as.integer(floor(x + 0.5))

#' Resample
#'
#' Imports the non-exported function `SpaDES.tools:::resample`.
#'
#' @importFrom utils getFromNamespace
#' @keywords internal
#' @rdname resample
#' @seealso [SpaDES.tools::resample()]
.resample <- getFromNamespace("resample", "SpaDES.tools")

#' Test whether disturbance should be scheduled
#'
#' @param disturbanceLayer a `RasterLayer` object
#' @param currentYear time of simulation
#'
#' @return Logical indicating whether to schedule a disturbance event
#'
#' @export
#' @examples
#' \dontrun{
#'   doEvent <- scheduleDisturbance(sim$rstCurrentBurn, time(sim), disturbanceType = "Burn")
#' }
scheduleDisturbance <- function(disturbanceLayer, currentYear) {
  if (is.null(disturbanceLayer) ||
      is.null(disturbanceLayer@data@attributes$Year) ||
      disturbanceLayer@data@attributes$Year != currentYear) {
    TRUE
  } else {
    FALSE
  }
}

#' Log-transformed values, with a floor (> 0)
#'
#' Avoid `-Inf` problems when `x == 0` by setting a non-zero floor value for `x`.
#' This is preferred over using some `log(x + d)` transformation, as the choice of `d` is
#' arbitrary, and will affect model fit.
#'
#' @param x     Numeric.
#' @param floor Minimum age value boundary. Default `0.3`.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   x <- sample(0:400, 1e7, replace = TRUE)
#'   floor <- 0.3
#'   logFloor <- log(floor)
#'   microbenchmark::microbenchmark(
#'     log(pmax(floor, x)),
#'     pmax(log(floor), log(x)),
#'     {x[x < floor] <- floor; log(x)},
#'     {y <- log(x); y[is.infinite(y)] <- logFloor} ## fastest; TODO: implement? not a bottleneck
#'   )
#' }
.logFloor <- function(x, floor = 0.3) {
  log(pmax(floor, x))
}

## TODO: rename the function "logTrunc"? implement a ceiling?

#' Search/return for current data layers
#'
#' @param dataLayers a named list/stack of `SpatRasters` with
#'   data layers to be searched (i.e., data layers for multiple simulation
#'   periods). Names *must* be suffixed with `"*_<year>"`, where
#'   "<year>" is the *first year* of the data period to use. For instance,
#'   if layers available are `"temperature_2021"` and `"temperature_2041"`
#'   and `currentYear` is 2025, `"temperature_2021"` will be used.
#'   If `currentYear` is 2050, `"temperature_2041"` will be used.
#'   If `currentYear` is 2020, the function will error as there are no
#'   climate layers available for periods before 2021.
#' @param currentYear the year for which we want find climate data for.
#'
#' @return a filtered list of data layers

currentDataLayers <- function(dataLayers, currentYear) {
  availYears <- as.integer(gsub("[^0-9]", "", names(dataLayers)))
  availYears <- unique(availYears[!is.na(availYears)])
  searchYears <- unique(sort(c(availYears, Inf)))

  if (currentYear < min(availYears)) {
    stop(paste0("The present simulation year, ", currentYear, ", is lower than the first climate",
                " projection year avalaible, ", min(availYears), ".\n  Please provide a climate",
                " projection layer that can be used for this period,\n  or consider increasing",
                " the simulation start year."))
  }
  climateYear <- as.character(cut(currentYear, breaks = searchYears, labels = availYears,
                                  right = FALSE))

  ## make a vector of years and select layers for a
  yearsOnly <- sub(".*_", "", names(dataLayers))

  return(dataLayers[yearsOnly == climateYear])
}

