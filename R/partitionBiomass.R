#' Partition biomass according to cover estimates
#'
#' This function will partition totalBiomass into each cohort. It will discount
#' deciduous cover, if x is < 1
#' @param x The ratio for deciduous cover:biomass, where conifer cover:biomass = 1
#' @param pixelCohortData A full pixelCohortData object (i.e., not cohortData)
#' @param decidSp A character vector of species that represent deciduous. Defaults
#'   to the equivalent names of \code{c("Popu_Tre", "Betu_Pap")}
#' @export
partitionBiomass <- function(x = 1, pixelCohortData, decidSp) {
  if (!"decid" %in% colnames(pixelCohortData)) {
    if (missing(decidSp)) {
      colName <- equivalentNameColumn(as.character(unique(pixelCohortData$speciesCode)), LandR::sppEquivalencies_CA)
      decidSp <- equivalentName(c("Popu_Tre", "Popu_Bal", "Betu_Pap"), LandR::sppEquivalencies_CA, colName)
    }
    pixelCohortData[, decid := speciesCode %in% decidSp]
  }

  pixelCohortData[, cover2 := cover * c(1,x)[decid + 1]]
  pixelCohortData[, cover2 := cover2/sum(cover2), by = "pixelIndex"]
  pixelCohortData[, B := totalBiomass*cover2]
  pixelCohortData

}
