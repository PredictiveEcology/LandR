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
  ## start with post-fire B - depending on the module used B
  ## may be present or not
  severityData <- copy(burnedPixelCohortData)

  if (!"B" %in% names(severityData)) {
    message("Assuming a stand replacing fire")
    severityData[, B := 0]
  }

  ## calculate post-fire stand biomass, drop unnecessary columns
  severityData <- severityData[, list(postfireB = sum(B),
                                      pixelGroup = pixelGroup), by = pixelIndex] %>%
    unique(.)

  ## sum initial B to stand level and add to severityData - this expands cohortData's B to pixels.
  severityData <- cohortData[, list(prefireB = sum(B)), by = pixelGroup][severityData, on = "pixelGroup"]

  ## calculate severity in terms of biomass
  severityData[, severityB := prefireB - postfireB]

  ## drop prefireB and postfireB columns
  return(severityData[, .(pixelGroup, pixelIndex, severityB)])
}
