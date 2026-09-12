utils::globalVariables(c(
  ".", ":=", "lightProb", "numberOfRegen",
  "resproutage_min", "resproutage_max", "sexualmature", "shadetolerance", "siteShade",
  "type", "year"
))

#' Activate serotiny after a (fire) disturbance
#'
#' @template burnedPixelCohortData
#' @template postFirePixelCohortData
#' @template postFireRegenSummary
#' @template species
#' @template sufficientLight
#' @template speciesEcoregion
#' @param currentTime integer. The current simulation time obtained with `time(sim)`
#' @param treedFirePixelTableSinceLastDisp a vector of pixels that burnt and were forested
#'     in the previous time step.
#' @template calibrate
#'
#' @return  A list of objects:
#'     `postFirePixelCohortData`, a `data.table` with the cohorts that undergo serotiny;
#'     `serotinyPixel`, a vector of pixels where serotiny was activated;
#'     `postFireRegenSummary`, the updated `postFireRegenSummary`, if `calibrate = TRUE`.
#'
#' @export
doSerotiny <- function(burnedPixelCohortData, postFirePixelCohortData,
                       postFireRegenSummary = NULL, species, sufficientLight,
                       speciesEcoregion, currentTime, treedFirePixelTableSinceLastDisp,
                       calibrate = FALSE) {
  ## checks
  if (calibrate && is.null(postFireRegenSummary)) {
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
      on = "speciesCode", nomatch = 0
    ]
    ## NOTE: should be in mortalityFromDisturbance module or event
    serotinyPixelCohortData <- serotinyPixelCohortData[age >= sexualmature] |>
      unique(by = c("pixelGroup", "speciesCode"))
    set(serotinyPixelCohortData, NULL, "sexualmature", NULL)

    ## select the pixels that have potential serotiny regeneration and assess them
    serotinyPixelTable <- treedFirePixelTableSinceLastDisp[
      pixelGroup %in% unique(serotinyPixelCohortData$pixelGroup)
    ]

    ## from now on the regeneration process is assessed for each potential pixel
    # setkey(serotinyPixelTable, pixelGroup)
    # setkey(serotinyPixelCohortData, pixelGroup)
    serotinyPixelCohortData <- serotinyPixelTable[serotinyPixelCohortData,
      allow.cartesian = TRUE,
      nomatch = 0, on = "pixelGroup"
    ] ## join table to add pixels

    ## light check: add shade tolerance to table and set shade to 0 (100% mortality).
    ## then get survival probs and subset survivors with runif
    serotinyPixelCohortData <- serotinyPixelCohortData[species[, .(speciesCode, shadetolerance)],
      nomatch = 0, on = "speciesCode"
    ]
    serotinyPixelCohortData <- assignLightProb(sufficientLight, serotinyPixelCohortData)
    serotinyPixelCohortData <- serotinyPixelCohortData[runif(nrow(serotinyPixelCohortData)) %<=% lightProb] ## subset survivors
    set(serotinyPixelCohortData, NULL, c("shadetolerance", "siteShade", "lightProb"), NULL) ## clean table again

    ## get establishment probs and subset species that establish with runif
    specieseco_current <- speciesEcoregion[year <= round(currentTime)]
    specieseco_current <- specieseco_current[
      year == max(specieseco_current$year),
      .(ecoregionGroup, speciesCode, establishprob)
    ]
    serotinyPixelCohortData <- serotinyPixelCohortData[specieseco_current,
      on = c("ecoregionGroup", "speciesCode"),
      nomatch = 0
    ]
    # serotinyPixelCohortData <- setkey(serotinyPixelCohortData, ecoregionGroup, speciesCode)[
    #   specieseco_current, nomatch = 0]  ## join table to add probs
    serotinyPixelCohortData <- serotinyPixelCohortData[runif(nrow(serotinyPixelCohortData)) %<=% establishprob][
      , establishprob := NULL
    ]

    ## only need one cohort per spp per pixel survives/establishes
    serotinyPixelCohortData <- unique(serotinyPixelCohortData, by = c("pixelIndex", "speciesCode"))

    if (NROW(serotinyPixelCohortData)) {
      serotinyPixelCohortData <- serotinyPixelCohortData[, .(
        pixelGroup, ecoregionGroup,
        speciesCode, pixelIndex
      )]
      serotinyPixelCohortData[, type := "serotiny"]
      if (calibrate) {
        serotinyRegenSummary <- serotinyPixelCohortData[, .(numberOfRegen = length(pixelIndex)),
          by = speciesCode
        ]
        serotinyRegenSummary <- serotinyRegenSummary[, .(
          year = currentTime, regenMode = "Serotiny",
          speciesCode, numberOfRegen
        )]
        serotinyRegenSummary <- setkey(serotinyRegenSummary, speciesCode)[
          species[, .(species, speciesCode)],
          nomatch = 0
        ]
        serotinyRegenSummary[, ":="(speciesCode = species, species = NULL)]
        setnames(serotinyRegenSummary, "speciesCode", "species")
        postFireRegenSummary <- rbindlist(list(postFireRegenSummary, serotinyRegenSummary))
      } else {
        postFireRegenSummary <- NULL
      }

      ## save the pixel index for resprouting assessment use,
      ## i.e., removing these pixel from assessing resprouting
      serotinyPixel <- unique(serotinyPixelCohortData$pixelIndex)

      ## append table to postFirePixelCohortData
      postFirePixelCohortData <- rbindlist(list(postFirePixelCohortData, serotinyPixelCohortData), fill = TRUE)
    } else {
      serotinyPixel <- NULL
    }
  } else {
    serotinyPixel <- NULL
  }

  return(list(
    postFirePixelCohortData = postFirePixelCohortData,
    serotinyPixel = serotinyPixel,
    postFireRegenSummary = postFireRegenSummary
  ))
}

#' Activate resprouting after a (fire) disturbance
#'
#' @template burnedPixelCohortData
#' @template postFirePixelCohortData
#' @template postFireRegenSummary
#' @param serotinyPixel a vector of pixels where serotiny was activated;
#' @template species
#' @template sufficientLight
#' @param currentTime integer. The current simulation time obtained with `time(sim)`
#' @param treedFirePixelTableSinceLastDisp a vector of pixels that burnt and were forested
#'        in the previous time step.
#' @template calibrate
#'
#' @return  A list of objects:
#'     `postFirePixelCohortData`, a `data.table` with the cohorts that undergo serotiny;
#'     `serotinyPixel`, a vector of pixels where serotiny was activated;
#'     `postFireRegenSummary`, the updated `postFireRegenSummary`, if `calibrate = TRUE`.
#'
#' @export
doResprouting <- function(burnedPixelCohortData, postFirePixelCohortData,
                          postFireRegenSummary = NULL, serotinyPixel,
                          treedFirePixelTableSinceLastDisp, currentTime,
                          species, sufficientLight, calibrate = FALSE) {
  ## checks
  if (calibrate && is.null(postFireRegenSummary)) {
    stop("missing postFireRegenSummary table for doResprouting")
  }

  ## make a table of pixels where resprouting occurs.
  if (is.null(serotinyPixel)) {
    resproutingPixelTable <- setkey(treedFirePixelTableSinceLastDisp, pixelGroup)
    # availableToResprout <- burnedPixelCohortData[0, ]
    availableToResprout <- copy(burnedPixelCohortData) ## Ceres - fix
  } else {
    ## Replacing here -- Eliot -- This was removing entire pixels that had successful serotiny
    ## -- now only species-pixel combos are removed.
    ## should be done by pixel and species -- Eliot: it works ok now because there are no serotinous
    ## species that are resprouters.
    full <- treedFirePixelTableSinceLastDisp[
      unique(burnedPixelCohortData,
        by = c("pixelGroup", "speciesCode")
      ),
      on = "pixelGroup", allow.cartesian = TRUE
    ]

    ## anti join to remove species-pixels that had successful serotiny/survivors
    ## Ceres: i don't know if I agree with this...
    availableToResprout <- full[!postFirePixelCohortData, on = c("pixelIndex", "speciesCode")]
  }

  ## assess whether reprouting can occur in burnt pixels
  species_temp <- species[
    postfireregen == "resprout",
    .(
      speciesCode, postfireregen,
      resproutage_min, resproutage_max, resproutprob
    )
  ]

  resproutingPixelCohortData <- availableToResprout[species_temp, nomatch = 0, on = "speciesCode"]
  resproutingPixelCohortData <- resproutingPixelCohortData[age >= resproutage_min & age <= resproutage_max]
  set(resproutingPixelCohortData, NULL, c("resproutage_min", "resproutage_max", "postfireregen", "age"), NULL)

  if (NROW(resproutingPixelCohortData)) {
    ## assess potential resprouting reg: add reprout probability, siteShade/tolerance to the table and assess who resprouts
    ## as long as one cohort can resprout, resprouting is activated
    # resproutingAssessCohortData <- unique(resproutingAssessCohortData, by = c("pixelGroup", "speciesCode"))
    # setkey(resproutingAssessCohortData, pixelGroup)

    ## make new table joing resprouters with burnt pixels
    # newPixelCohortData <- resproutingPixelTable[resproutingAssessCohortData, nomatch = 0, allow.cartesian = TRUE]

    ## light check: add shade tolerance to table and set shade to 0 (100% mortality.)
    ## the get survival probs and subset survivors with runif
    resproutingPixelCohortData <- resproutingPixelCohortData[
      species[, .(speciesCode, shadetolerance)],
      nomatch = 0, on = "speciesCode"
    ]
    # resproutingPixelCohortData[,siteShade := 0]    ## no longer part of resprouting
    # resproutingPixelCohortData <- setkey(resproutingPixelCohortData, speciesCode)[species[,.(speciesCode, shadetolerance)],
    #                                                     nomatch = 0][, siteShade := 0]
    resproutingPixelCohortData <- assignLightProb(
      sufficientLight = sufficientLight,
      resproutingPixelCohortData
    )

    resproutingPixelCohortData <- resproutingPixelCohortData[lightProb %>>% runif(nrow(resproutingPixelCohortData), 0, 1)]
    resproutingPixelCohortData <- resproutingPixelCohortData[resproutprob %>>% runif(nrow(resproutingPixelCohortData), 0, 1)]

    resproutingPixelCohortData <- unique(resproutingPixelCohortData, by = c("pixelIndex", "speciesCode"))
    set(resproutingPixelCohortData, NULL, c("resproutprob", "shadetolerance", "siteShade", "lightProb"), NULL)

    # remove all columns that were used temporarily here
    if (NROW(resproutingPixelCohortData)) {
      resproutingPixelCohortData <- resproutingPixelCohortData[, .(pixelGroup, ecoregionGroup, speciesCode, pixelIndex)]
      resproutingPixelCohortData[, type := "resprouting"]
      if (calibrate) {
        resproutRegenSummary <- resproutingPixelCohortData[
          , .(numberOfRegen = length(pixelIndex)),
          by = speciesCode
        ]
        resproutRegenSummary <- resproutRegenSummary[, .(
          year = currentTime, regenMode = "Resprout",
          speciesCode, numberOfRegen
        )]
        resproutRegenSummary <- setkey(resproutRegenSummary, speciesCode)[
          species[, .(species, speciesCode)],
          nomatch = 0
        ]
        resproutRegenSummary[, ":="(speciesCode = species, species = NULL)]
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

  return(list(
    postFirePixelCohortData = postFirePixelCohortData,
    serotinyResproutSuccessPixels = serotinyResproutSuccessPixels,
    postFireRegenSummary = postFireRegenSummary
  ))
}
