utils::globalVariables(c(
  ":=", ".SD", "Area", "col1", "firetolerance", "growthcurve", "hardsoft",
  "leafLignin", "leaflongevity", "mortalityshape",
  "seeddistance_eff", "seeddistance_max", "species", "species1", "species2",
  "wooddecayrate"
))

#' Default LANDIS-II project repo url
#'
#' @keywords internal
landisIIrepo <- paste0("https://raw.githubusercontent.com/LANDIS-II-Foundation/",
                       "Extensions-Succession/master/biomass-succession-archive/",
                       "trunk/tests/v6.0-2.0/")

#' Download and prepare a species traits table for use with `Biomass_core` module
#'
#' `prepSpeciesTable`
#'
#' @note This one is tailored to Canadian forests (?)
#'
#' @param url If NULL (the default), uses one from D. Cyr's LANDIS-II files:
#' <https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv>).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to `Cache`.
#'
#' @export
#' @importFrom data.table data.table
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
                        fun = "data.table::fread",
                        header = TRUE, stringsAsFactors = FALSE,
                        userTags = c(cacheTags, "speciesTable"))

  return(speciesTable)
}

#' Species Table Column Names
#'
#' @keywords internal
.speciesTableColNames <- c("species", "Area", "longevity", "sexualmature", "shadetolerance",
                           "firetolerance", "seeddistance_eff", "seeddistance_max", "resproutprob",
                           "resproutage_min", "resproutage_max", "postfireregen", "leaflongevity",
                           "wooddecayrate", "mortalityshape", "growthcurve", "leafLignin",
                           "hardsoft")

#' @template speciesTable
#'
#' @param speciesLayers Deprecated.
#' @param areas A character vector of areas to use. Can be one or more of
#'   `c("Acadian", "AM", "NorthShore", "BP", "BSE", "BSW", "LSJ", "MC", "PM", "WestON")`.
#'   If it is more than one, this function will take the minimum value, within a species.
#'   These are short versions of the Canada Ecoprovinces.
#'   Currently defaults to `c("BSW", "BP", "MC")` for historical reasons.
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @return A `data.table` with columns ... TODO
#'
#' @export
#' @rdname speciesTable
prepSpeciesTable <- function(speciesTable, speciesLayers = NULL,
                             sppEquiv = NULL, sppEquivCol = "LandR",
                             areas = c("BSW", "BP", "MC")) {
  if (!"Area" %in% names(speciesTable)) {
    stop("Please add an 'Area' column of ecoprovinces to 'sim$speciesTable'")
  }

  if (is.null(sppEquiv))
    sppEquiv <- data.table(utils::data("sppEquivalencies_CA", package = "LandR", envir = environment()))

  names(speciesTable) <- .speciesTableColNames

  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  speciesTable <- speciesTable[species %in% equivalentName(sppNameVector, sppEquiv, "LANDIS_traits", multi = TRUE) &
                                 Area %in% areas]

  speciesTable[, species := equivalentName(speciesTable$species, sppEquiv, sppEquivCol)]
  speciesTable <- speciesTable[, lapply(.SD, function(x) {
    if (is.numeric(x)) min(x, na.rm = TRUE) else x[1]
  }), by = "species"]

  ## use integers (instead of numerics) where possible; these are asserted in Biomass_core
  speciesTable[, `:=`(Area = as.factor(Area),
                      growthcurve = as.numeric(growthcurve),
                      shadetolerance = as.numeric(shadetolerance),
                      hardsoft = as.factor(hardsoft),
                      seeddistance_eff = asInteger(seeddistance_eff),
                      seeddistance_max = asInteger(seeddistance_max),
                      resproutage_min = asInteger(resproutage_min),
                      resproutage_max = asInteger(resproutage_max),
                      mortalityshape = asInteger(mortalityshape),
                      postfireregen = as.factor(postfireregen))]

  return(speciesTable[])
}

#' Change species table of parameters/traits
#'
#' Changes longevity and shade tolerance values in the species table.
#' Longevity values are changed to follow Burton & Cumming (1995) for the following species:
#' *Abies balsamea*, *Abies lasiocarpa*, *Betula papyrifera*, *Larix laricina*,
#' *Larix occidentalis*, *Picea engelmannii*, *Picea glauca*, *Picea mariana*,
#' *Pinus banksiana*, *Pinus contorta*, *Pinus resinosa*, *Pinus strobus*,
#' *Populus balsamifera v. balsamifera*, *Populus tremuloides*, *Pseudotsuga menziesii var. glauca*,
#' *Pseudotsuga menziesii*, *Thuja plicata*, *Tsuga heterophylla*,
#' *Tsuga mertensiana x heterophylla*, and only for the  Boreal Shield West (BSW), Boreal Plains (BP)
#'  and Montane Cordillera (MC) `speciesTable$Area`s.
#' Note that BSW and BP areas correspond more closely to the region considered in Table 2 of
#' Burton & Cumming (1995), while MC will correspond to both tables.
#'
#' Of the above species, shade tolerance values are changed for *Abies spp*, *Picea spp*,
#' and *Tsuga spp.* to reflect western boreal shade tolerances better.
#'
#' When different longetivity/shade tolerance trait values exist for a given species, the minimum
#' value across `Area`'s (BSW, BP, MC) is kept.
#'
#' ATTENTION: if none of species in `species` are from BSW, BP or MC area this function will not
#' change any values.
#'
#' All other species/Area trait values follow Dominic Cyr and Yan Boulanger's trait values available at:
#' (<https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv>).
#'
#'
#' @template species
#'
#' @template speciesTable
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @return An updated species `data.table`
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom crayon red
#'
#' @rdname speciesTableUpdate
speciesTableUpdate <- function(species, speciesTable, sppEquiv, sppEquivCol) {
  ## if "Area"is a column in the (final) traits table, then check and warn the user for area
  ## mismatches
  if (!"Area" %in% names(species)) {
    stop("Can't find 'Area' column in 'sim$species'")

    test <- !any(unique(species$Area) %in% c("BSW", "BP", "MC"))
    if (test) {
      message(red("Areas in 'species$Area' do not match any of 'BSW', 'BP' or 'MC',",
                  "\nno changes made to 'sim$species'."))
      return(species)
    }
  }

  if (is.null(sppEquiv))
    sppEquiv <- data.table(utils::data("sppEquivalencies_CA", package = "LandR", envir = environment()))

  if (is.null(sppEquivCol))
    stop("Please provide sppEquivCol")

  names(speciesTable) <- .speciesTableColNames

  ## make temporary table that will have new parameters for Boreal spp.
  ## longevity values from Burton & Cumming (1995)
  speciesTableShort <- speciesTable[Area %in% c("BSW", "BP", "MC"), .(species, longevity, shadetolerance)]
  speciesTableShort[species == "ABIE.BAL", longevity := 200] #default 150
  speciesTableShort[species == "ABIE.LAS", longevity := 240] #default 250
  speciesTableShort[species == "BETU.PAP", longevity := 140] #default 150
  speciesTableShort[species == "LARI.LAR", longevity := 350] #default 150
  speciesTableShort[species == "LARI.OCC", longevity := 450] #default 900!
  speciesTableShort[species == "PICE.ENG", longevity := 460] #default 450
  speciesTableShort[species == "PICE.GLA", longevity := 400] #default 250
  speciesTableShort[species == "PICE.MAR", longevity := 250] #default 200
  speciesTableShort[species == "PINU.BAN", longevity := 150]  #default 150 - no change
  speciesTableShort[species == "PINU.CON.LAT", longevity := 335] #default 300
  speciesTableShort[species == "PINU.PON", longevity := 575] #default 500
  speciesTableShort[species == "POPU.BAL", longevity := 200] #default 130
  speciesTableShort[species == "POPU.TRE", longevity := 200] #default 150
  speciesTableShort[species == "PSEU.MEN", longevity := 525] ## default 600, only in MC area, corresponding to var. glauca
  speciesTableShort[species == "THUJ.PLI", longevity := 1500] ##default 700, 1500 may be incorrect for MC
  speciesTableShort[species == "TSUG.HET", longevity := 500] #default 475
  speciesTableShort[species == "TSUG.MER", longevity := 800] #default 700

  speciesTableShort[species == "ABIE.BAL", shadetolerance := 3] #default 5
  speciesTableShort[species == "ABIE.LAS", shadetolerance := 3] #default 4
  speciesTableShort[species == "PICE.ENG", shadetolerance := 3] #default 4
  speciesTableShort[species == "PICE.GLA", shadetolerance := 2] #default 3
  speciesTableShort[species == "PICE.MAR", shadetolerance := 3] #default 4
  speciesTableShort[species == "TSUG.HET", shadetolerance := 4] #default 5
  speciesTableShort[species == "TSUG.MER", shadetolerance := 3] #default 4

  ## subset, rename and "merge" species by using the minimum value
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]
  sppNameVector <- species$species
  speciesTableShort <- speciesTableShort[species %in% equivalentName(sppNameVector, sppEquiv,
                                                                     "LANDIS_traits", multi = TRUE)]
  speciesTableShort[, species := equivalentName(speciesTableShort$species, sppEquiv, sppEquivCol)]
  speciesTableShort <- speciesTableShort[, .(longevity = min(longevity),
                                             shadetolerance = min(shadetolerance)), by = "species"]

  ## join to deal with eventual non-matching species ordering
  ## subset species table to common species, then add missing species lines
  ## (which did not have traits changed above)
  cols <- setdiff(names(species), c("longevity", "shadetolerance"))
  speciesTemp <- species[, ..cols]
  speciesTemp <- speciesTableShort[speciesTemp, on = "species", nomatch = 0]
  species <- rbind(species[!species %in% speciesTemp$species], speciesTemp)[order(species)]

  ## make sure updated columns have the correct class
  species[, `:=`(
    Area = as.factor(Area),
    firetolerance = asInteger(firetolerance),
    growthcurve = as.numeric(growthcurve),
    longevity = asInteger(longevity),
    leaflongevity = asInteger(leaflongevity),
    leafLignin = as.numeric(leafLignin),
    mortalityshape = asInteger(mortalityshape),
    postfireregen = as.factor(postfireregen),
    resproutage_max = asInteger(resproutage_max),
    resproutage_min = asInteger(resproutage_min),
    resproutprob = as.numeric(resproutprob),
    seeddistance_eff = asInteger(seeddistance_eff),
    seeddistance_max = asInteger(seeddistance_max),
    sexualmature = asInteger(sexualmature),
    shadetolerance = as.numeric(shadetolerance),
    species = as.character(species),
    wooddecayrate = as.numeric(wooddecayrate)
  )]

  return(species)
}

#' Download and prepare a species traits table for use with `Biomass_core` module
#'
#' TODO: add detailed description
#'
#' @note This one is tailored to Canadian forests (?)
#'
#' @param url If NULL (the default), uses one from the LANDIS-II project:
#' <https://github.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession_test.txt">).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to `Cache`.
#'
#' @export
#' @importFrom data.table data.table setcolorder
#' @importFrom reproducible asPath Cache prepInputs
#'
#' @return A `data.table` with columns ... TODO
#'
#' @export
#' @rdname prepInputsSpecies
prepInputsSpecies <- function(url = NULL, dPath, cacheTags = NULL) {
  if (is.null(url))
    url <- paste0(landisIIrepo, "species.txt")

  mainInput <- prepInputsMainInput(url = NULL, dPath, cacheTags) ## uses default URL

  maxcol <- 13 #max(count.fields(file.path(dPath, "species.txt"), sep = ""))
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
  species <- species[col1 != "LandisData", ]
  species <- species[col1 != ">>", ]
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
prepInputsMainInput <- function(url = NULL, dPath = tempdir(), cacheTags = NULL) {
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
#' <https://github.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/ecoregion.txt">).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to `Cache`.
#'
#' @return A `data.table`
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
#' <https://github.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/biomass-succession-dynamic-inputs_test.txt">).
#'
#' @param dPath The destination path.
#'
#' @param cacheTags User tags to pass to `Cache`.
#'
#' @return A `data.table`
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
