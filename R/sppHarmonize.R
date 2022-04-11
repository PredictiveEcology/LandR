#' Harmonize the 3 components that bring species into Biomass_**
#'
#' This function will attempt to harmonize many potential issues/conflicts that
#' may arise under different combinations of supplied objects to the three
#' arguments. See manual for details.
#'
#' @param sppEquiv A data.table or NULL representing the sppEquiv table (see manual).
#'   If NULL, then it will use the `LandR::sppEquivalencies_CA`
#' @param sppNameVector A character vector of species to use. These species must all
#'   be from one naming convention, i.e., from one column in the sppEquiv.
#' @param sppEquivCol A character string normally provided from the P(sim)$sppEquivCol
#'   (see manual). If `NA`, the default, then this will try to determine which
#'   column the `sppNameVector` used and use that. If `sppNameVector` is NULL, then
#'   it will default to `"Boreal"`.
#' @return Returns a named list with the same names as the arguments. These should
#'   likely be assigned to the `sim` object in the module following this function call.
#' @export
#' @examples
#' sppOuts <- sppHarmonize(sim$sppEquiv, sim$sppNameVector, P(sim)$sppEquivCol)
#' sim$sppEquiv <- sppOuts$sppEquiv
#' sim$sppNameVector <- sppOuts$sppNameVector
#' P(sim)$sppEquivCol <- sppOuts$sppEquivCol
#'
sppHarmonize <- function(sppEquiv, sppNameVector, sppEquivCol) {
  if (is.null(sppEquiv) && is.null(sppNameVector)) {
    sppNameConvention <- "KNN"
    sppNameVector <-  Cache(LandR::speciesInStudyArea, studyArea)$speciesList
  }

  if (is.null(sppEquiv)) {
    # can't use suppliedElsewhere because the sequence below needs this;
    # can't specify where = "user" either because it could be supplied by another module prior to this
    sppEquiv <- LandR::sppEquivalencies_CA
  }

  if (is.na(sppEquivCol)) {
    if (is.null(sppNameVector)) {
      sppEquivCol <- "Boreal"
    } else {
      sppNameConvention <- LandR::equivalentNameColumn(sppNameVector, sppEquiv)
      sppEquivCol <- sppNameConvention
    }
  }

  if (is.null(sppNameVector)) {
    sppNameVector <- sort(unique(sppEquiv[[sppEquivCol]]))
  } else {
    sppNameConvention <- LandR::equivalentNameColumn(sppNameVector, sppEquiv)
  }

  if (!exists("sppNameConvention", inherits = FALSE)) {
    sppNameConvention <- sppEquivCol
  }
  sppMissing <- sppNameVector %in% sppEquiv[[sppNameConvention]]
  sppMissing <- sppNameVector[!(sppMissing)]

  if (length(sppMissing)) {
    stop("Different naming convention detected for species (", paste(sppMissing, collapse = ", "), ")",
         " in sim$sppNameVector. Please ensure all species use a single column in spp$sppEquiv ",
         "or LandR::sppEquivalencies_CA")
  }


  sppNameVectorInSppEquivCol <- if (!identical(sppEquivCol, sppNameConvention)) {
    LandR::equivalentName(sppNameVector, df = sppEquiv,
                          column = sppEquivCol, searchColumn = sppNameConvention)
  } else {
    sppNameVector
  }

  sppMissing <- sppNameVectorInSppEquivCol %in% sppEquiv[[sppEquivCol]]
  sppMissing <- sppNameVector[!(sppMissing)]

  if (length(sppMissing)) {
    stop("Species missing (", paste(sppMissing, collapse = ", "), ")",
         " from sim$sppEquiv in P(sim)$sppEquivCol (", sppEquivCol, ")")
  }

  sppEquiv <- sppEquiv[get(sppEquivCol) %in% sppNameVectorInSppEquivCol]

  sppEquiv <- sppEquiv[!"", on = sppEquivCol]
  sppEquiv <- na.omit(sppEquiv, sppEquivCol)

  return(list(sppEquiv = sppEquiv, sppNameVector = sppNameVector, sppEquivCol = sppEquivCol))
}
