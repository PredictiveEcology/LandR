if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(":=", ".SD", "Area", "col1", "growthcurve",
                           "leafLignin", "leaflongevity", "mortalityshape",
                           "seeddistance_eff", "seeddistance_max",
                           "species", "species1", "species2", "wooddecayrate"))
}

#' Default LANDIS-II project repo url
#'
#' @keywords internal
landisIIrepo <- paste0("https://raw.githubusercontent.com/LANDIS-II-Foundation/",
                      "Extensions-Succession/master/biomass-succession-archive/",
                      "trunk/tests/v6.0-2.0/")

#' Download and prepare a species traits table for use with LBMR module
#'
#' TODO: add detailed description
#'
#' @note This one is tailored to Canadian forests (?)
#'
#' @param url If NULL (the default), uses one from D. Cyr's LANDIS-II files:
#' \url{https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv}).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to \code{Cache}.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom magrittr %>%
#' @importFrom reproducible asPath Cache prepInputs
#' @rdname speciesTable
getSpeciesTable <- function(url = NULL, dPath = tempdir(), cacheTags = NULL) {
  if (is.null(url))
    url <- paste0("https://raw.githubusercontent.com/",
                  "dcyr/LANDIS-II_IA_generalUseFiles/",
                  "master/speciesTraits.csv")

  speciesTable <- Cache(prepInputs, "speciesTraits.csv",
                        destinationPath = asPath(dPath),
                        url = url,
                        fun = "utils::read.csv",
                        header = TRUE, stringsAsFactors = FALSE,
                        userTags = c(cacheTags, "speciesTable")) %>%
    data.table()

  return(speciesTable)
}

#' @param speciesTable  A raw species traits table
#'
#' @param speciesLayers TODO: description needed
#'
#' @param sppEquiv TODO: description needed
#'
#' @param sppEquivCol TODO: description needed
#'
#' @return A \code{data.table} with columns ... TODO
#'
#' @export
#' @rdname speciesTable
prepSpeciesTable <- function(speciesTable, speciesLayers, sppEquiv = NULL, sppEquivCol = "LandR") {

  if (is.null(sppEquiv))
    sppEquiv <- data.table(utils::data("sppEquivalencies_CA",
                                                 package = "pemisc",
                                                 envir = environment()))

  names(speciesTable) <- c(
    "species",
    "Area",
    "longevity",
    "sexualmature",
    "shadetolerance",
    "firetolerance",
    "seeddistance_eff",
    "seeddistance_max",
    "resproutprob",
    "resproutage_min",
    "resproutage_max",
    "postfireregen",
    "leaflongevity",
    "wooddecayrate",
    "mortalityshape",
    "growthcurve",
    "leafLignin",
    "hardsoft"
  )

  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]),]
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  speciesTable <- speciesTable[species %in% equivalentName(sppNameVector, sppEquiv,
                                                           "LANDIS_traits", multi = TRUE) &
                                 Area %in% c("BSW", "BP", "MC")]

  speciesTable[, species := equivalentName(speciesTable$species, sppEquiv, sppEquivCol)]
  speciesTable <- speciesTable[, lapply(.SD, function(x) {
    if (is.numeric(x)) min(x, na.rm = TRUE) else x[1]
  }), by = "species"]

  return(speciesTable)
}

#' Download and prepare a species traits table for use with \code{LandR_BiomassGMOrig} module
#'
#' TODO: add detailed description
#'
#' @note This one is tailored to Canadian forests (?)
#'
#' @param url If NULL (the default), uses one from the LANDIS-II project:
#' \url{https://github.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession_test.txt"}).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to \code{Cache}.
#'
#' @export
#' @importFrom data.table data.table setcolorder
#' @importFrom reproducible asPath Cache prepInputs
#'
#' @return A \code{data.table} with columns ... TODO
#'
#' @export
#' @rdname prepInputsSpecies
prepInputsSpecies <- function(url = NULL, dPath, cacheTags = NULL) {
  if (is.null(url))
    url <- paste0(landisIIrepo, "species.txt")

  mainInput <- prepInputsMainInput(url = NULL, dPath, cacheTags) ## uses default URL

  maxcol <- 13#max(count.fields(file.path(dPath, "species.txt"), sep = ""))
  species <- Cache(prepInputs,
                   url = url,
                   targetFile = "species.txt",
                   destinationPath = dPath,
                   fun = "utils::read.table",
                   fill = TRUE, row.names = NULL, #purge = 7,
                   sep = "",
                   header = FALSE,
                   blank.lines.skip = TRUE,
                   col.names = c(paste("col", 1:maxcol, sep = "")),
                   stringsAsFactors = FALSE,
                   overwrite = TRUE)
  species <- data.table(species[, 1:11])
  species <- species[col1 != "LandisData",]
  species <- species[col1 != ">>",]
  colNames <- c("species", "longevity", "sexualmature", "shadetolerance",
                "firetolerance", "seeddistance_eff", "seeddistance_max",
                "resproutprob", "resproutage_min", "resproutage_max",
                "postfireregen")
  names(species) <- colNames
  species[, ':='(seeddistance_eff = gsub(",", "", seeddistance_eff),
                 seeddistance_max = gsub(",", "", seeddistance_max))]
  # change all columns to integer
  species <- species[, lapply(.SD, as.integer), .SDcols = names(species)[-c(1, NCOL(species))],
                     by = "species,postfireregen"]
  setcolorder(species, colNames)

  # get additional species traits
  speciesAddon <- mainInput
  startRow <- which(speciesAddon$col1 == "SpeciesParameters")
  speciesAddon <- speciesAddon[(startRow + 1):(startRow + nrow(species)), 1:6, with = FALSE]
  names(speciesAddon) <- c("species", "leaflongevity", "wooddecayrate",
                           "mortalityshape", "growthcurve", "leafLignin")
  speciesAddon[, ':='(leaflongevity = as.numeric(leaflongevity),
                      wooddecayrate = as.numeric(wooddecayrate),
                      mortalityshape = as.numeric(mortalityshape),
                      growthcurve = as.numeric(growthcurve),
                      leafLignin = as.numeric(leafLignin))]

  species <- setkey(species, species)[setkey(speciesAddon, species), nomatch = 0]

  ## TODO: use species equivalency table here
  ## rename species for compatibility across modules (Genu_spe)
  species$species1 <- as.character(substring(species$species, 1, 4))
  species$species2 <- as.character(substring(species$species, 5, 7))
  species[, ':='(species = paste0(toupper(substring(species1, 1, 1)),
                                  substring(species1, 2, 4), "_",
                                  species2))]

  species[, ':='(species1 = NULL, species2 = NULL)]

  return(species)
}

#' @export
#' @rdname prepInputsSpecies
prepInputsMainInput <- function(url = NULL, dPath, cacheTags) {
  if (is.null(url))
    url <- paste0("https://raw.githubusercontent.com/LANDIS-II-Foundation/",
                  "Extensions-Succession/master/biomass-succession-archive/",
                  "trunk/tests/v6.0-2.0/biomass-succession_test.txt")

  maxcol <- 7L
  mainInput <- Cache(prepInputs,
                     url = url,
                     targetFile = "biomass-succession_test.txt",
                     destinationPath = dPath,
                     userTags = cacheTags,
                     fun = "utils::read.table",
                     fill = TRUE,  #purge = 7,
                     sep = "",
                     header = FALSE,
                     col.names = c(paste("col", 1:maxcol, sep = "")),
                     blank.lines.skip = TRUE,
                     stringsAsFactors = FALSE)

  mainInput <- data.table(mainInput)
  mainInput <- mainInput[col1 != ">>",]

  return(mainInput)
}

#' Prepare ecoregion table
#'
#' Get the dummy ecoregion table from LANDIS-II examples.
#'
#' @param url If NULL (the default), uses one from the LANDIS-II project:
#' \url{https://github.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/ecoregion.txt"}).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to \code{Cache}.
#'
#' @return A \code{data.table}
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom reproducible Cache prepInputs
#' @importFrom utils count.fields
#' @rdname prepInputsEcoregion
prepInputsEcoregion <- function(url = NULL, dPath, cacheTags = NULL) {
  if (is.null(url))
    url <- paste0(landisIIrepo, "ecoregions.txt")

  maxcol <- 5 #max(count.fields(file.path(dPath, "ecoregions.txt"), sep = ""))
  ecoregion <- Cache(prepInputs,
                     url = url,
                     targetFile = "ecoregions.txt",
                     destinationPath = dPath,
                     fun = "utils::read.table",
                     fill = TRUE,
                     sep = "",
                     # purge = 7,
                     header = FALSE,
                     blank.lines.skip = TRUE,
                     stringsAsFactors = FALSE,
                     userTags = cacheTags)
  maxcol <- max(count.fields(file.path(dPath, "ecoregions.txt"), sep = ""))
  colnames(ecoregion) <- c(paste("col", 1:maxcol, sep = ""))
  ecoregion <- data.table(ecoregion)
  ecoregion <- ecoregion[col1 != "LandisData",]
  ecoregion <- ecoregion[col1 != ">>",]
  names(ecoregion)[1:4] <- c("active", "mapcode", "ecoregion", "description")
  ecoregion$mapcode <- as.integer(ecoregion$mapcode)

  return(ecoregion)
}

#' Prepare species ecoregion table
#'
#' Get the dummy ecoregion table from LANDIS-II examples.
#'
#' @param url If NULL (the default), uses one from the LANDIS-II project:
#' \url{https://github.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession-dynamic-inputs_test.txt"}).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to \code{Cache}.
#'
#' @return A \code{data.table}
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom reproducible Cache prepInputs
#' @importFrom utils count.fields
#' @rdname prepInputsSpeciesEcoregion
prepInputsSpeciesEcoregion <- function(url = NULL, dPath, cacheTags = NULL) {
  if (is.null(url))
    url <- paste0(landisIIrepo, "biomass-succession-dynamic-inputs_test.txt")

  speciesEcoregion <- Cache(prepInputs,
                            url = url,
                            fun = "utils::read.table",
                            destinationPath = dPath,
                            targetFile = "biomass-succession-dynamic-inputs_test.txt",
                            fill = TRUE,
                            sep = "",
                            header = FALSE,
                            blank.lines.skip = TRUE,
                            stringsAsFactors = FALSE,
                            userTags = cacheTags)
  maxcol <- max(count.fields(file.path(dPath, "biomass-succession-dynamic-inputs_test.txt"),
                             sep = ""))
  colnames(speciesEcoregion) <- paste("col", 1:maxcol, sep = "")
  speciesEcoregion <- data.table(speciesEcoregion)
  speciesEcoregion <- speciesEcoregion[col1 != "LandisData",]
  speciesEcoregion <- speciesEcoregion[col1 != ">>",]
  keepColNames <- c("year", "ecoregion", "species", "establishprob", "maxANPP", "maxB")
  names(speciesEcoregion)[1:6] <- keepColNames
  speciesEcoregion <- speciesEcoregion[, keepColNames, with = FALSE]
  integerCols <- c("year", "establishprob", "maxANPP", "maxB")
  speciesEcoregion[, (integerCols) := lapply(.SD, as.integer), .SDcols = integerCols]

  ## TODO: use species equivalency table here
  ## rename species for compatibility across modules (Genu_spe)
  speciesEcoregion$species1 <- as.character(substring(speciesEcoregion$species, 1, 4))
  speciesEcoregion$species2 <- as.character(substring(speciesEcoregion$species, 5, 7))
  speciesEcoregion[, ':='(species = paste0(toupper(substring(species1, 1, 1)),
                                           substring(species1, 2, 4), "_", species2))]

  speciesEcoregion[, ':='(species1 = NULL, species2 = NULL)]

  return(speciesEcoregion)
}
