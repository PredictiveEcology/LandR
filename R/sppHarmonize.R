#' Harmonize the three components that bring species into `Biomass_**` modules
#'
#' This function will attempt to harmonize many potential issues/conflicts that
#' may arise under different combinations of supplied objects to the three
#' arguments. See manual for details.
#'
#' @template sppEquiv
#'
#' @template sppNameVector
#'
#' @param sppEquivCol A character string normally provided from the `P(sim)$sppEquivCol`
#'   (see manual). If `NA`, the default, then this will try to determine which
#'   column the `sppNameVector` used and use that. If `sppNameVector` is NULL, then
#'   it will default to `"Boreal"`.
#' @param dPath Passed to `speciesInStudyArea` (which then passes to `preProcess`)
#'
#' @template sppColorVect
#'
#' @template vegLeadingProportion
#'
#' @template studyArea
#'
#' @return Returns a named list with the same names as the arguments. These should
#'   likely be assigned to the `sim` object in the module following this function call.
#' @export
#' @examples
#' ## not run. usage example within module
#' # sppOuts <- sppHarmonize(sim$sppEquiv, sim$sppNameVector, P(sim)$sppEquivCol)
#' # sim$sppEquiv <- sppOuts$sppEquiv
#' # sim$sppNameVector <- sppOuts$sppNameVector
#' # P(sim)$sppEquivCol <- sppOuts$sppEquivCol
#'
sppHarmonize <- function(sppEquiv, sppNameVector, sppEquivCol, sppColorVect,
                         vegLeadingProportion = 0, studyArea,
                         dPath = getOption("reproducible.destinationPath")) {
  if (is.null(sppEquiv) && is.null(sppNameVector)) {
    sppNameConvention <- "KNN"
    if (!missing(studyArea)) {
      sppNameVector <-  Cache(LandR::speciesInStudyArea, studyArea, dPath = dPath)$speciesList
    } else {
      stop("Please pass studyArea polygon if not passing sppEquiv/sppNameVector")
    }
  }

  # Catch both cases of sppEquiv is NULL or is an object with NROW 0 e.g., NULL data.table
  if (is.null(sppEquiv))
    sppEquiv <- data.table()
  if (NROW(sppEquiv) == 0) {
    ## note that this step MUST come after the previous
    sppEquiv <- LandR::sppEquivalencies_CA
  }

  if (is.na(sppEquivCol)) {
    if (is.null(sppNameVector)) {
      sppEquivCol <- "Boreal"
    } else {
      sppNameConvention <- LandR::equivalentNameColumn(sppNameVector, sppEquiv)
      sppEquivCol <- sppNameConvention
    }
  } else {
    if (sppEquivCol %in% names(sppEquiv)) {
      message(paste("Using species listed under", sppEquivCol,
                    "column in the 'sppEquiv' table"))
    } else {
      stop("'sppEquivCol' is not a column in 'sppEquiv' or 'LandR::sppEquivalencies_CA'.",
           " Please provide conforming 'sppEquivCol' and 'sppEquiv',",
           " or a column that exists in 'LandR::sppEquivalencies_CA'")
    }
  }

  if (is.null(sppNameVector)) {
    sppNameVector <- sort(unique(sppEquiv[[sppEquivCol]]))
    sppNameVector <- sppNameVector[nzchar(sppNameVector)] ## drop empty strings
  } else {
    sppNameConvention <- LandR::equivalentNameColumn(sppNameVector, sppEquiv)
  }

  if (!exists("sppNameConvention", inherits = FALSE)) {
    sppNameConvention <- sppEquivCol
  }
  sppMissing <- !(sppNameVector %in% sppEquiv[[sppNameConvention]])
  sppMissing <- sppNameVector[sppMissing]

  if (length(sppMissing)) {
    stop("Species in 'sppNameVector' (", paste(unique(sppMissing), collapse = ", "), ") ",
         "not found in 'sppEquiv', or 'LandR::sppEquivalencies_CA'. This may be ",
         "due to different species in each set, or same species with different ",
         "naming conventions. Please ensure all species correspond to a column ",
         "in 'sppEquiv' or 'LandR::sppEquivalencies_CA'.")
  }


  sppNameVectorInSppEquivCol <- if (!identical(sppEquivCol, sppNameConvention)) {
    LandR::equivalentName(sppNameVector, df = sppEquiv,
                          column = sppEquivCol, searchColumn = sppNameConvention)
  } else {
    sppNameVector
  }

  sppMissing <- !(sppNameVectorInSppEquivCol %in% sppEquiv[[sppEquivCol]])
  sppMissing <- sppNameVector[sppMissing]

  if (length(sppMissing)) {
    stop("Species missing (", paste(sppMissing, collapse = ", "), ")",
         " from 'sppEquiv' in 'sppEquivCol' (", sppEquivCol, ")")
  }

  sppEquiv <- sppEquiv[get(sppEquivCol) %in% sppNameVectorInSppEquivCol]

  sppEquiv <- sppEquiv[!"", on = sppEquivCol]
  sppEquiv <- na.omit(sppEquiv, sppEquivCol)

  ## sanity checks/warnings. the code above allows several species
  ## from sppNameVector to match one spp in sppEquivCol (e.g. merged spp)
  if (isFALSE(all(sppNameVector %in% sppEquiv[[sppEquivCol]]))) {
    warning("Some species in 'sppNameVector' are not in 'sppEquiv[[sppEquivCol]]'.",
            "This may be because species are being merged. Please verify.")
  }

  ## add default colors for species used in model
  if (is.null(sppColorVect)) {
    sppColorVect <- sppColors(sppEquiv, sppEquivCol, newVals = "Mixed", palette = "Accent")
    message("No 'sppColorVect' provided; making one with colour palette: Accent")
  } else {
    if (any(grepl("Mixed", names(sppColorVect)))) {
      if (length(sppColorVect) != (length(unique(sppEquiv[[sppEquivCol]])) + 1) ||
          length(sppColorVect) != (length(sppNameVector) + 1)) {
        stop("Length of 'sppColorVect' differs from number species in final 'sppEquiv'.",
             " Check species provided in 'sppColorVect', 'sppNameVector' and/or ",
             "'sppEquiv[[sppEquivCol]]'")
      }
    } else {
      if (length(sppColorVect) != length(unique(sppEquiv[[sppEquivCol]])) ||
          length(sppColorVect) != length(sppNameVector)) {
        stop("Length of 'sppColorVect' differs from number species in final 'sppEquiv'.",
             " Check species provided in 'sppColorVect', 'sppNameVector' and/or ",
             "'sppEquiv[[sppEquivCol]]'")
      }
    }
    if (!all(names(sppColorVect)[names(sppColorVect) != "Mixed"] %in% sppEquiv[[sppEquivCol]])) {
      message("'sppColorVect' names do not match 'sppEquivCol' (", sppEquivCol, ")",
              " and will be converted if possible.")
      tempNames <-  LandR::equivalentName(names(sppColorVect), df = sppEquiv,
                                          column = sppEquivCol)
      if (any(is.na(tempNames))) {
        missingNames <- names(sppColorVect)[is.na(tempNames)]
        stop("Could not convert 'sppColorVect' names. The following are missing",
             " from 'sppEquiv' in 'sppEquivCol' (", sppEquivCol, "):",
             paste(missingNames, collapse = ", "))
      }
      names(sppColorVect) <- tempNames
    }
  }

  ## vegLeadingProportion can be NULL below:
  if (isTRUE(vegLeadingProportion > 0) && isTRUE(is.na(sppColorVect["Mixed"]))) {
    stop("'vegLeadingProportion'  is > 0 but there is no 'Mixed' color in 'sppColorVect'. ",
         "Please supply 'sppColorVect' with a 'Mixed' color or set 'vegLeadingProportion' to zero.")
  }

  return(list(sppEquiv = sppEquiv, sppNameVector = sppNameVector,
              sppEquivCol = sppEquivCol, sppColorVect = sppColorVect))
}
