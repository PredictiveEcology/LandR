utils::globalVariables(c("Broadleaf", "cover", "cover2", "decid", "LandR", "totalBiomass"))

#' Partition biomass according to cover estimates
#'
#' This function will partition `totalBiomass` into each cohort.
#' It will discount deciduous cover, if `x < 1`.
#' @param x The ratio for deciduous cover:biomass, where conifer cover:biomass = 1
#' @param pixelCohortData A full `pixelCohortData` object (i.e., not `cohortData`)
#' @export
partitionBiomass <- function(x = 1, pixelCohortData) {
  if (!"decid" %in% colnames(pixelCohortData)) {
    sppEquivalencies_CA <- get(data("sppEquivalencies_CA", package = "LandR",
                                    envir = environment()), inherits = FALSE)

    colName <- equivalentNameColumn(as.character(unique(pixelCohortData$speciesCode)),
                                    sppEquivalencies_CA)
    decidSp <- equivalentName(sppEquivalencies_CA[Broadleaf == TRUE, LandR],
                              sppEquivalencies_CA, colName)
    decidSp <- decidSp[nzchar(decidSp)]
    pixelCohortData[, decid := speciesCode %in% decidSp]
  }

  pixelCohortData[, cover2 := cover * c(1, x)[decid + 1]]
  pixelCohortData[, cover2 := cover2 / sum(cover2), by = "pixelIndex"]
  pixelCohortData[, B := totalBiomass * cover2]

  set(pixelCohortData, NULL, c("decid", "cover2"), NULL)

  pixelCohortData
}
