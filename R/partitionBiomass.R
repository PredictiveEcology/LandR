utils::globalVariables(c("Broadleaf", "cover", "cover2", "decid", "LandR", "totalBiomass"))

#' Partition biomass according to cover estimates
#'
#' This function will partition \code{totalBiomass} into each cohort.
#' It will discount deciduous cover, if \code{x < 1}.
#' @param x The ratio for deciduous cover:biomass, where conifer cover:biomass = 1
#' @param pixelCohortData A full \code{pixelCohortData} object (i.e., not \code{cohortData})
#' @export
partitionBiomass <- function(x = 1, pixelCohortData) {
  if (!"decid" %in% colnames(pixelCohortData)) {
    colName <- equivalentNameColumn(as.character(unique(pixelCohortData$speciesCode)), LandR::sppEquivalencies_CA)
    decidSp <- equivalentName(LandR::sppEquivalencies_CA[Broadleaf == TRUE, LandR],
                              LandR::sppEquivalencies_CA, colName)
    decidSp <- decidSp[nzchar(decidSp)]
    pixelCohortData[, decid := speciesCode %in% decidSp]
  }

  pixelCohortData[, cover2 := cover * c(1, x)[decid + 1]]
  pixelCohortData[, cover2 := cover2 / sum(cover2), by = "pixelIndex"]
  pixelCohortData[, B := totalBiomass * cover2]
  pixelCohortData
}
