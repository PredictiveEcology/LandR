utils::globalVariables(c(
  ".", ":=", "B", "pixelGroup", "speciesCode",
  "prefireB", "postfireB", "severityB"
))

#' Calculate fire severity
#'
#' Calculates fire severity as the loss of pre-fire to
#'   post-fire biomass.
#'
#' @template cohortData
#' @template burnedPixelCohortData
#'
#' @note if \code{burnedPixelCohortData} does not have a \code{B}
#'    column, the fire is assumed to be stand replacing (i.e.
#'    we assume B to be 0 across all pixels/cohorts in
#'    \code{burnedPixelCohortData})
#'
#' @export
#' @return \code{data.table} with columns \code{pixelIndex},
#'   \code{pixelGroup} and  \code{severityB}
#'
#'
calcSeverityB <- function(cohortData, burnedPixelCohortData) {
  if (!"B" %in% names(burnedPixelCohortData)) {
    message("Assuming a stand replacing fire")
    burnedPixelCohortData[, B := 0]
  }

  severityData <- burnedPixelCohortData[, .(pixelIndex, pixelGroup)]

  ## add initial and post-fire B to severityData
  severityData <- cohortData[, .(speciesCode, B, pixelGroup)][severityData, on = "pixelGroup"]
  setnames(severityData, "B", "prefireB")

  severityData <- burnedPixelCohortData[, .(speciesCode, B, pixelIndex)][severityData, on = .(speciesCode, pixelIndex)]
  setnames(severityData, "B", "postfireB")

  ## sum B's across species, and drop species
  severityData[, `:=`(prefireB = sum(prefireB),
                      postfireB = sum(postfireB)), by = pixelIndex]
  set(severityData, j = "speciesCode", value = NULL)
  severityData <- unique(severityData)

  ## calculate severity in terms of biomass
  severityData[, severityB := prefireB - postfireB]

  ## keep only certain columns
  return(severityData[, .(pixelGroup, pixelIndex, severityB)])
}
