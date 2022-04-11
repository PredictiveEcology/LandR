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
