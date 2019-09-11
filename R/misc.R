if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", "lightProb", "shadetolerance", "siteShade", "year"))
}

#' Assign light probability
#'
#' @param sufficientLight a data.table containing probability of establishment, given a site's light conditions (X0-X5) for
#'    each level of a species shade tolerance (1-5)
#' @param newCohortData  a modified version of /code{cohortData} that contains new cohorts.
#'
#' @return  /code{newCohortData} with a lightProb column
#'
#' @export
assignLightProb <- function(sufficientLight, newCohortData) {
  ## for each line, get the survival probability from sufficientLight table note
  ## that sufficentLight is a table of survival probs for each tolerance level
  ## (row) by and shade level (column) siteShade + 2 is necessary to skip the
  ## first column
  newCohortData[, lightProb := sufficientLight[cbind(shadetolerance, siteShade + 2)]]
}

#' Convert numeric values to rounded integers
#'
#' Essentially a wrapper around \code{round}, rather than \code{truncate}, which is what \code{as.integer}
#' does. Internally, this is simply \code{as.integer(floor(x + 0.5))}.
#'
#' @note Values ending in \code{.5} will be rounded up, whether positive or negative.
#' This is different than \code{round}.
#'
#' @param x A numeric vector
#'
#' @return An integer vector of length \code{x}, rounded to zero decimal places
#'   prior to \code{as.integer}
#'
#' @export
#' @examples
#' x <- seq(-2, 2, 0.25)
#' data.frame(dbl = x, int = asInteger(x))
asInteger <- function(x)
  as.integer(floor(x + 0.5))

#' Resample
#'
#' Imports the non-exported function \code{SpaDES.tools:::resample}.
#'
#' @importFrom utils getFromNamespace
#' @keywords internal
#' @rdname resample
#' @seealso \code{\link[SpaDES.tools]{resample}}
.resample <- getFromNamespace("resample", "SpaDES.tools")

#' Test whether disturbance should be scheduled
#'
#' @param disturbanceLayer a \code{RasterLayer} object
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


#' Activate serotiny after a (fire) disturbance
#'
#' @param burnedPixelCohortData A /code{cohortData} table expanded to pixelIndex containing the only the cohorts in
#'    the pixels that were affected by the disturbance
#' @param postFirePixelCohortData  an empty /code{cohortData}-like table with columns "age", "B", "mortality",
#'    "aNPPAct", and "sumB" removed and "pixelIndex" added.
#' @param postFireRegenSummary a data.table summarizing for which species serotiny/resprouting were
#'    activated and in how many pixels, for each year. Only necessary if /code{calibrate = TRUE}.
#' @param species a \code{data.table} with species traits such as longevity, shade tolerance, etc.
#' @param sufficientLight a data.table containing probability of establishment, given a site's light conditions (X0-X5) for
#'    each level of a species shade tolerance (1-5)
#' @param speciesEcoregion a \code{data.table} with \code{speciesEcoregion} values
#' @param currentTime integer. The current simulation time obtained with /code{time(sim)}
#' @param treedFirePixelTableSinceLastDisp a vector of pixels that burnt and were forested in the previous time step.
#' @param calibrate logical. Determines whether to output /code{postFirePixelCohortData}. Defaults to FALSE
#'
#' @return  A list of objects. /code{postFirePixelCohortData}, a data.table with the cohorts that will undergo serotiny;
#'    /code{serotinyPixel} a vector of pixels where serotiny was activated;
#'    /code{postFireRegenSummary} the updated postFireRegenSummary, if /code{calibrate = TRUE}
#'
#' @export
#'

doSerotiny <- function(burnedPixelCohortData, postFirePixelCohortData,
                       postFireRegenSummary = NULL, species, sufficientLight,
                       speciesEcoregion, currentTime, treedFirePixelTableSinceLastDisp,
                       calibrate = FALSE) {
  ## checks
  if (calibrate & is.null(postFireRegenSummary)) {
    stop("missing postFireRegenSummary table for doSerotiny")
  }

  ## subset spp with serotiny
  tempspecies <- species[postfireregen == "serotiny", .(speciesCode, postfireregen)]

  ## join tables to make a serotiny table
  serotinyPixelCohortData <- burnedPixelCohortData[tempspecies, nomatch = 0][, postfireregen := NULL]

  if (NROW(serotinyPixelCohortData)) {
    ## assess potential serotiny reg: add sexual maturity to the table and compare w/ age
    ## as long as one cohort is sexually mature, serotiny is activated
    serotinyPixelCohortData <- serotinyPixelCohortData[species[, .(speciesCode, sexualmature)],
                                                       on = "speciesCode", nomatch = 0]
    #serotinyPixelCohortData <- setkey(serotinyPixelCohortData, speciesCode)[species[,.(speciesCode, sexualmature)],
    #                                                                          nomatch = 0]
    serotinyPixelCohortData <- serotinyPixelCohortData[age >= sexualmature] %>% # NOTE should be in mortalityFromDisturbance module or event
      unique(., by = c("pixelGroup", "speciesCode"))
    set(serotinyPixelCohortData, NULL, "sexualmature", NULL)

    ## select the pixels that have potential serotiny regeneration and assess them
    serotinyPixelTable <- treedFirePixelTableSinceLastDisp[pixelGroup %in% unique(serotinyPixelCohortData$pixelGroup)]

    ## from now on the regeneration process is assessed for each potential pixel
    #setkey(serotinyPixelTable, pixelGroup)
    #setkey(serotinyPixelCohortData, pixelGroup)
    serotinyPixelCohortData <- serotinyPixelTable[serotinyPixelCohortData, allow.cartesian = TRUE,
                                                  nomatch = 0, on = "pixelGroup"] ## join table to add pixels

    ## light check: add shade tolerance to table and set shade to 0 (100% mortality.)
    ## the get survival probs and subset survivors with runif
    serotinyPixelCohortData <- serotinyPixelCohortData[species[, .(speciesCode, shadetolerance)],
                                                       nomatch = 0, on = "speciesCode"]
    # serotinyPixelCohortData[, siteShade := 0]
    # serotinyPixelCohortData <- setkey(serotinyPixelCohortData, speciesCode)[species[,.(speciesCode, shadetolerance)],
    #                                                     nomatch = 0][, siteShade := 0]
    serotinyPixelCohortData <- assignLightProb(sufficientLight = sufficientLight,
                                               serotinyPixelCohortData)
    serotinyPixelCohortData <- serotinyPixelCohortData[lightProb %>>% runif(nrow(serotinyPixelCohortData), 0, 1)]  ## subset survivors
    set(serotinyPixelCohortData, NULL, c("shadetolerance", "siteShade", "lightProb"), NULL)   ## clean table again

    ## get establishment probs and subset species that establish with runif
    specieseco_current <- speciesEcoregion[year <= round(currentTime)]
    specieseco_current <- specieseco_current[year == max(specieseco_current$year),
                                             .(ecoregionGroup, speciesCode, establishprob)]
    serotinyPixelCohortData <- serotinyPixelCohortData[specieseco_current, on = c("ecoregionGroup", "speciesCode"), nomatch = 0]
    #serotinyPixelCohortData <- setkey(serotinyPixelCohortData, ecoregionGroup, speciesCode)[specieseco_current, nomatch = 0]  ## join table to add probs
    serotinyPixelCohortData <- serotinyPixelCohortData[runif(nrow(serotinyPixelCohortData), 0, 1) %<<% establishprob][, establishprob := NULL]

    ## only need one cohort per spp per pixel survives/establishes
    serotinyPixelCohortData <- unique(serotinyPixelCohortData, by = c("pixelIndex", "speciesCode"))

    if (NROW(serotinyPixelCohortData)) {
      ## rm age
      serotinyPixelCohortData <- serotinyPixelCohortData[,.(pixelGroup, ecoregionGroup, speciesCode, pixelIndex)] #
      serotinyPixelCohortData[, type := "serotiny"]
      if (calibrate) {
        serotinyRegenSummary <- serotinyPixelCohortData[,.(numberOfRegen = length(pixelIndex)), by = speciesCode]
        serotinyRegenSummary <- serotinyRegenSummary[,.(year = currentTime, regenMode = "Serotiny",
                                                        speciesCode, numberOfRegen)]
        serotinyRegenSummary <- setkey(serotinyRegenSummary, speciesCode)[species[,.(species, speciesCode)],
                                                                          nomatch = 0]
        serotinyRegenSummary[, ':='(speciesCode = species, species = NULL)]
        setnames(serotinyRegenSummary, "speciesCode", "species")
        postFireRegenSummary <- rbindlist(list(postFireRegenSummary, serotinyRegenSummary))
      } else {
        postFireRegenSummary <- NULL
      }
      serotinyPixel <- unique(serotinyPixelCohortData$pixelIndex) # save the pixel index for resprouting assessment use,
      # i.e., removing these pixel from assessing resprouting
      ## append table to postFirePixelCohortData
      postFirePixelCohortData <- rbindlist(list(postFirePixelCohortData, serotinyPixelCohortData), fill = TRUE)
    } else {
      serotinyPixel <- NULL
    }
  } else {
    serotinyPixel <- NULL
  }

  return(list(postFirePixelCohortData = postFirePixelCohortData,
              serotinyPixel = serotinyPixel,
              postFireRegenSummary = postFireRegenSummary))
}

#' Activate resprouting after a (fire) disturbance
#'
#' @param burnedPixelCohortData A /code{cohortData} table expanded to pixelIndex containing the only the cohorts in
#'    the pixels that were affected by the disturbance
#' @param postFirePixelCohortData  an empty /code{cohortData}-like table with columns "age", "B", "mortality",
#'    "aNPPAct", and "sumB" removed and "pixelIndex" added.
#' @param postFireRegenSummary a data.table summarizing for which species serotiny/resprouting were
#'    activated and in how many pixels, for each year. Only necessary if /code{calibrate = TRUE}.
#' @param species a \code{data.table} with species traits such as longevity, shade tolerance, etc.
#' @param sufficientLight a data.table containing probability of establishment, given a site's light conditions (X0-X5) for
#'    each level of a species shade tolerance (1-5)
#' @param speciesEcoregion a \code{data.table} with \code{speciesEcoregion} values
#' @param currentTime integer. The current simulation time obtained with /code{time(sim)}
#' @param treedFirePixelTableSinceLastDisp a vector of pixels that burnt and were forested in the previous time step.
#' @param calibrate logical. Determines whether to output /code{postFirePixelCohortData}. Defaults to FALSE
#'
#' @return  A list of objects. /code{postFirePixelCohortData}, a data.table with the cohorts that will undergo serotiny;
#'    /code{serotinyResproutSuccessPixels} a vector of pixels where serotiny and resprouting were activated;
#'    /code{postFireRegenSummary} the updated postFireRegenSummary, if /code{calibrate = TRUE}
#' @export

doResprouting <- function(burnedPixelCohortData, postFirePixelCohortData,
                          postFireRegenSummary = NULL, serotinyPixel,
                          treedFirePixelTableSinceLastDisp, currentTime,
                          species, sufficientLight, calibrate = FALSE) {
  ## checks
  if (calibrate & is.null(postFireRegenSummary)) {
    stop("missing postFireRegenSummary table for doResprouting")
  }

  ## make a table of pixels where resprouting occurs.
  if (is.null(serotinyPixel)) {
    resproutingPixelTable <- setkey(treedFirePixelTableSinceLastDisp, pixelGroup)
    # availableToResprout <- burnedPixelCohortData[0,]
    availableToResprout <- copy(burnedPixelCohortData)    ## Ceres - fix

  } else {
    # Replacing here -- ELiot -- THis was removing entire pixels that had successful serotiny -- now only species-pixel combos are removed
    ## should be done by pixel and species -- Eliot: it works ok now because there are no serotinous species that are resprouters
    full <- treedFirePixelTableSinceLastDisp[unique(burnedPixelCohortData, by = c("pixelGroup", "speciesCode")),
                                             on = "pixelGroup", allow.cartesian = TRUE] #

    # anti join to remove species-pixels that had successful serotiny/survivors
    # Ceres: i don't know if I agree with this...
    availableToResprout <- full[!postFirePixelCohortData, on = c("pixelIndex", "speciesCode")]
  }

  ## assess whether reprouting can occur in burnt pixels
  species_temp <- species[postfireregen == "resprout",
                          .(speciesCode, postfireregen,
                            resproutage_min, resproutage_max, resproutprob)]

  resproutingPixelCohortData <- availableToResprout[species_temp, nomatch = 0, on = "speciesCode"]
  resproutingPixelCohortData <- resproutingPixelCohortData[age >= resproutage_min & age <= resproutage_max]
  set(resproutingPixelCohortData, NULL, c("resproutage_min", "resproutage_max", "postfireregen", "age"), NULL)

  if (NROW(resproutingPixelCohortData)) {
    ## assess potential resprouting reg: add reprout probability, siteShade/tolerance to the table and assess who resprouts
    ## as long as one cohort can resprout, resprouting is activated
    #resproutingAssessCohortData <- unique(resproutingAssessCohortData, by = c("pixelGroup", "speciesCode"))
    #setkey(resproutingAssessCohortData, pixelGroup)

    ## make new table joing resprouters with burnt pixels
    #newPixelCohortData <- resproutingPixelTable[resproutingAssessCohortData, nomatch = 0, allow.cartesian = TRUE]

    ## light check: add shade tolerance to table and set shade to 0 (100% mortality.)
    ## the get survival probs and subset survivors with runif
    resproutingPixelCohortData <- resproutingPixelCohortData[species[, .(speciesCode, shadetolerance)],
                                                             nomatch = 0, on = "speciesCode"]
    # resproutingPixelCohortData[,siteShade := 0]
    # resproutingPixelCohortData <- setkey(resproutingPixelCohortData, speciesCode)[species[,.(speciesCode, shadetolerance)],
    #                                                     nomatch = 0][, siteShade := 0]
    resproutingPixelCohortData <- assignLightProb(sufficientLight = sufficientLight,
                                                  resproutingPixelCohortData)

    resproutingPixelCohortData <- resproutingPixelCohortData[lightProb %>>% runif(nrow(resproutingPixelCohortData), 0, 1)]
    resproutingPixelCohortData <- resproutingPixelCohortData[resproutprob %>>% runif(nrow(resproutingPixelCohortData), 0, 1)]

    resproutingPixelCohortData <- unique(resproutingPixelCohortData, by = c("pixelIndex", "speciesCode"))
    set(resproutingPixelCohortData, NULL, c("resproutprob", "shadetolerance", "siteShade", "lightProb"), NULL)

    # remove all columns that were used temporarily here
    if (NROW(resproutingPixelCohortData)) {
      resproutingPixelCohortData <- resproutingPixelCohortData[,.(pixelGroup, ecoregionGroup, speciesCode, pixelIndex)]#
      resproutingPixelCohortData[, type := "resprouting"]
      if (calibrate) {
        resproutRegenSummary <- resproutingPixelCohortData[,.(numberOfRegen = length(pixelIndex)), by = speciesCode]
        resproutRegenSummary <- resproutRegenSummary[,.(year = currentTime, regenMode = "Resprout",
                                                        speciesCode, numberOfRegen)]
        resproutRegenSummary <- setkey(resproutRegenSummary, speciesCode)[species[,.(species, speciesCode)],
                                                                          nomatch = 0]
        resproutRegenSummary[,':='(speciesCode = species, species = NULL)]
        setnames(resproutRegenSummary, "speciesCode", "species")
        postFireRegenSummary <- rbindlist(list(postFireRegenSummary, resproutRegenSummary))
      } else {
        postFireRegenSummary <- NULL
      }
      ## append resprouters to the table
      postFirePixelCohortData <- rbindlist(list(postFirePixelCohortData, resproutingPixelCohortData), fill = TRUE)
      postFirePixelCohortData[, type := factor(type)]
      serotinyResproutSuccessPixels <- c(serotinyPixel, unique(resproutingPixelCohortData$pixelIndex))
    } else {
      serotinyResproutSuccessPixels <- serotinyPixel
    }
  } else {
    serotinyResproutSuccessPixels <- serotinyPixel
  }

  return(list(postFirePixelCohortData = postFirePixelCohortData,
              serotinyResproutSuccessPixels = serotinyResproutSuccessPixels,
              postFireRegenSummary = postFireRegenSummary))
}
