utils::globalVariables(c(
  "allPres", "allPresFac", "KNN", "pixel", "variable", "NTEMS_Species_Code"
))

####kNN ####

maskTo <- utils::getFromNamespace("maskTo", "reproducible")
projectTo <- utils::getFromNamespace("projectTo", "reproducible")

#' Make a species factor raster from KNN
#'
#' This will download all KNN layers in Forests of Canada, and make
#' a factor raster at resolution provided by `res` (larger is faster).
#'
#' @param year Default (and only implemented) is 2011. This will download the 2011 KNN data layers
#'
#' @param dPath A character string indicating where to download all the KNN layers
#'
#' @param res The resolution (one dimension, in m) for the resulting raster
#'
#' @param minPctCover An integer indicating what percent cover a species must have
#'   in a pixel to be considered present in that pixel.
#'
#' @return A `SpatRaster` object with 2 layers: `"speciesPresent"` is a factor, with
#' a legend (i.e., it is numbers on a map, that correspond to a legend) and
#' `"numberSpecies"` which represents the number of species in each pixel.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("googledrive", quietly = TRUE)) {
#'   # Make the dataset
#'   speciesPresent <- speciesPresentFromKNN(dPath = "~/data/KNN")
#'
#'   # To upload this:
#'   speciesPresentRas <- terra::rast(speciesPresent)[[1]]
#'   fn <- "SpeciesPresentInCanadianForests.tif"
#'   writeRaster(speciesPresentRas, file = fn)
#'   zipFn <- gsub(".tif", ".zip", fn)
#'   zip(files = dir(pattern = fn), zipFn)
#'   out <- googledrive::drive_put(zipFn)
#'
#'   ## Note: previous file (cropped to boreal forest)
#'   ## available at "1Oj78jJBeha5L6XDBBdWDAfimgNjYc9UD"
#'   driveID <- "1J8fN7clZeqjd7yhiDWi13uoCBL8OensF"
#'
#'   ## Get species list
#'   sa <- LandR::randomStudyArea(size = 1e11)
#'   species <- LandR::speciesInStudyArea(sa)
#' }
#' }
#'
#' @export
#' @rdname speciesPresent
speciesPresentFromKNN <- function(year = 2011, dPath = asPath("."), res = 2000, minPctCover = 10) {
  studyAreaED <- Cache(
    prepInputs,
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    destinationPath = dPath,
    # fun = quote(SA_ERIntersect(x = targetFilePath, studyArea)),
    overwrite = FALSE
  )

  opts <- options("reproducible.useTerra" = TRUE)
  on.exit(options(opts), add = TRUE)
  studyAreaER <- Cache(
    prepInputs,
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecodistrict_shp.zip",
    destinationPath = dPath,
    fun = "terra::vect",
    overwrite = TRUE
  )

  templateCRS <- reproducible::prepInputs(
    url = paste0(
      "https://www12.statcan.gc.ca/census-recensement/2021/",
      "geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip"
    ),
    destinationPath = dPath
  )
  sa <- projectTo(studyAreaER, terra::crs(templateCRS))

  allForestedStk <- Cache(loadAndAggregateKNN, dPath, res, sa)
  allForestedStk <- round(allForestedStk, 0)
  allForestedStk[allForestedStk <= minPctCover] <- 0

  numSp <- sum(allForestedStk > 0)

  mat <- terra::values(allForestedStk)
  dt <- as.data.table(mat)
  dt[, pixel := seq_len(.N)]
  dt2 <- melt(dt, measure.vars = setdiff(colnames(dt), "pixel"), na.rm = TRUE, id.vars = "pixel")
  dt2 <- dt2[value != 0]
  setorderv(dt2, c("pixel", "variable"))
  dt3 <- dt2[, list(allPres = paste(variable, collapse = "__")), by = "pixel"]
  dt3[, allPresFac := factor(allPres)]

  # Create a new empty rast
  speciesPres <- terra::rast(allForestedStk[[1]])
  # fill it with the integer values
  speciesPres[dt3$pixel] <- as.integer(dt3$allPresFac)
  names <- unique(dt3$allPresFac)
  numerics <- as.integer(names)
  # assign the levels
  levels(speciesPres) <- data.frame(ID = numerics, category = names)


  return(c(speciesPres, numSp))
}


####NTEMS ####

#
#' Make a species factor raster based on NTEMS Data
#'
#' This will download NTEMS dominant species layer for 2011 for forests of Canada, and make
#' a factor raster at resolution provided by `res` (larger is faster).
#'
#' @param dPath A character string indicating where to download all the KNN layers
#'
#' @param res The resolution (one dimension, in m) for the resulting raster
#'
#' @param year Default (and only implemented) is 2011. This will download the 2011 NTEMS species layer
#'
#' @template rasterToMatch
#'
#' @template studyArea
#'
#' @param ... Additional arguments passed to [reproducible::Cache()]
#'            and [equivalentName()]. Also valid: `outputPath`, and `studyAreaName`.
#'
#' @return A `SpatRaster` object with 2 layers: `"speciesPresent"` is a factor, with
#' a legend (i.e., it is numbers on a map, that correspond to a legend) and
#' `"numberSpecies"` which represents the number of species in each pixel.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("googledrive", quietly = TRUE)) {
#'   # Make the dataset
#'   speciesPresent <- speciesPresentFromNTEMS(dPath = "~/data/NTEMS")
#'
#'   # To upload this:
#'   speciesPresentRas <- terra::rast(speciesPresent)[[1]]
#'   fn <- "SpeciesPresentInCanadianForests_NTEMS.tif"
#'   writeRaster(speciesPresentRas, file = fn)
#'   zipFn <- gsub(".tif", ".zip", fn)
#'   zip(files = dir(pattern = fn), zipFn)
#'   out <- googledrive::drive_put(zipFn)
#'
#'   ## Get species list
#'   sa <- LandR::randomStudyArea(size = 1e11)
#'   species <- LandR::speciesInStudyArea(sa, dataSource = "NTEMS")
#' }
#' }
#'
#' @export
#' @rdname speciesPresent
speciesPresentFromNTEMS <- function(dPath = asPath("."), res = 2400, year = 2011,
                                    rasterToMatch = NULL, studyArea = NULL, ...) {
  dots <- list(...)

  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop(
      "Suggested package 'terra' not installed.\n",
      "Install it using `install.packages('terra')`."
    )
  }
  oPath <- if (!is.null(dots$outputPath)) dots$outputPath else dPath

  if ("cachePath" %in% names(dots)) {
    cachePath <- dots$cachePath
  } else {
    cachePath <- getOption("reproducible.cachePath")
  }

  sppEquiv <- get(data("sppEquivalencies_CA",
    package = "LandR",
    envir = environment()
  ), inherits = FALSE)

  opts <- options("reproducible.useTerra" = TRUE)
  on.exit(options(opts), add = TRUE)

  SppURL <- paste0("https://opendata.nfis.org/downloads/forest_change/CA_Tree_Species_Classification_", year, ".zip")
  SppRast <- prepInputs(
    url = SppURL,
    destinationPath = dPath
  )
  #For some reason the raster is not trimmed - this reduces the pixel count by 30%
  SppRast <- terra::trim(SppRast)
  NAflag(SppRast) <- 0 #so we don't make an NA raster when we use aggregate/segregate

  aggName <- file.path(dPath, paste0("CA_Tree_Species_NTEMS_", year, "_120m.tif"))

  NAflag(SppRast) <- 0

  SppRast <- terra::aggregate(SppRast, fact = 4, fun = "modal", na.rm = TRUE,
                              filename = aggName, overwrite = TRUE)
  uniqueVals <- as.data.table(terra::unique(SppRast)) #for later

  segFile <- file.path(dPath, paste0("CA_Tree_Species", year, "_120m_seg.tif"))

  SppRast <- terra::segregate(SppRast, classes = uniqueVals$NTEMS_Species_Code,
                              keep = TRUE, other = 0, filename = segFile,
                              overwrite = TRUE)

  newAggName <- file.path(dPath, paste0("CA_Tree_Species", year, "_", res, ".tif"))

  SppRast <- terra::aggregate(SppRast, ceiling(res/120), fun = "max", overwrite = TRUE, na.rm = TRUE,
                              filename = newAggName)

  sppEquiv <- sppEquiv[, .SD, .SDcol = c("NTEMS_Species_Code", "LandR")] # matching NTEMS spp code to sppEquivCol


  setnames(uniqueVals, new = "NTEMS_Species_Code")

  uniqueVals <- sppEquiv[uniqueVals, on = c("NTEMS_Species_Code")] # pulling all species from NTEMS layer
  uniqueVals <- unique(uniqueVals[, .(NTEMS_Species_Code, LandR)])

  #note that Pinu_alb and Pice_abi both dropped
  ActualNTEMS <- as.integer(names(SppRast))
  uniqueVals <- uniqueVals[NTEMS_Species_Code %in% ActualNTEMS,]
  names(SppRast) <- c(uniqueVals$LandR)

  #one last time because aggregate or segregate reintroduces zeroes
  NAflag(SppRast) <- 0

  numSp <- sum(SppRast > 0, na.rm = TRUE)

  mat <- terra::values(SppRast)
  dt <- as.data.table(mat)
  dt[, pixel := seq_len(.N)]
  dt2 <- melt(dt, measure.vars = setdiff(colnames(dt), "pixel"), na.rm = TRUE, id.vars = "pixel")
  dt2 <- dt2[value != 0]
  setorderv(dt2, c("pixel", "variable"))
  dt3 <- dt2[, list(allPres = paste(variable, collapse = "__")), by = "pixel"]
  dt3[, allPresFac := factor(allPres)]

  # Create a new empty rast
  speciesPres <- terra::rast(SppRast[[1]])
  # fill it with the integer values
  speciesPres[dt3$pixel] <- as.integer(dt3$allPresFac)
  names <- unique(dt3$allPresFac)
  numerics <- as.integer(names)
  # assign the levels
  levels(speciesPres) <- data.frame(ID = numerics, category = names)


  return(c(speciesPres, numSp))
}

#' Get species list in a given study area for a forest in Canada
#'
#' `speciesInStudyArea` defaults to use a url of a dataset uploaded to Google Drive that is
#' from Canadian Forests, but a different factor raster can be passed e.g.,
#' from `speciesPresentFromKNN`.
#'
#' @template studyArea
#'
#' @param url A url to get a `speciesPresence` raster e.g., from `peciesPresentFromKNN`
#'
#' @param speciesPresentRas A factor raster where the character string is a string of
#'   species names, separated by 2 underscores, sorted alphabetically. Can be produced
#'   with `speciesPresentFromKNN`
#'
#' @param sppEquivCol An optional column from `LandR::sppEquivalencies_CA`.
#'   If passed the KNN species will be returned according to this naming convention.
#'
#' @param dataSource Character. Either KNN or NTEMS. Defaults to KNN to obtain species from layer
#'   created using species cover data from KNN. Also able to obtain species from NTEMS
#'   Dominant Species Layer (though less species will likely be included).
#'
#' @param dPath Passed to `destinationPath` in `preProcess`.
#'
#' @return A named list of length 2: `speciesRas` is a factor `RasterLayer`
#' and `speciesList` is a character string containing the unique, sorted
#' species on the `speciesRas`, for convenience.
#'
#' @export
speciesInStudyArea <- function(studyArea, url = NULL, speciesPresentRas = NULL, sppEquivCol = NULL,
                               dataSource = "KNN",
                               dPath = getOption("reproducible.destinationPath")) {
  if (!(dataSource %in% c("KNN", "NTEMS"))) {
    stop("Data Source must be either KNN or NTEMS")
  }
  if (dataSource == "KNN") {
    if (is.null(speciesPresentRas)) {
      if (is.null(url)) {
        url <- "https://drive.google.com/file/d/1J8fN7clZeqjd7yhiDWi13uoCBL8OensF"
      }
      speciesPres <- preProcess(url = url, destinationPath = dPath)
      speciesPresRas <- rasterRead(speciesPres$targetFilePath)
    } else {
      speciesPresRas <- speciesPresentRas
    }

    bb <- postProcess(x = speciesPresRas, studyArea = studyArea)

    rasLevs <- as.data.table(levels(bb))
    # if (is(speciesPresRas, "RasterLayer")) {
    #   bb <- raster::deratify(bb)
    # }
    IDcol <- names(rasLevs)[1]
    speciesCommunities <- na.omit(rasLevs[rasLevs[[IDcol]] %in% as.vector(bb[[1]])]$category)
    species <- as.character(speciesCommunities)
    species <- unique(unlist(strsplit(species, "__")))

    if (!is.null(sppEquivCol) & is.null(speciesPresentRas)) {
      sppEquiv <- LandR::sppEquivalencies_CA
      species <- unique(sppEquiv[KNN %in% species, .SD, ][[sppEquivCol]])
      species <- species[!species == ""]
    }
  }
  if (dataSource == "NTEMS") {
    warning("the NTEMS data includes only dominant species and thus using this dataset will likely result in fewer species returned")
    if (is.null(speciesPresentRas)) {
      if (is.null(url)) {
        url <- "https://drive.google.com/file/d/1AytH5-4FZ3XpOFOiVFqUCLrMtavNmHuj"
      }
      speciesPres <- preProcess(url = url, destinationPath = dPath)
      speciesPresRas <- rasterRead(speciesPres$targetFilePath)
    } else {
      speciesPresRas <- speciesPresentRas
    }

    bb <- postProcess(x = speciesPresRas, studyArea = studyArea)
    rasLevs <- as.data.table(levels(speciesPresRas[[1]]))

    speciesCommunities <- na.omit(rasLevs[rasLevs[["value"]] %in% as.vector(bb[[1]])]$category)

    # else
    if (!is.null(sppEquivCol)) {
      sppEquiv <- LandR::sppEquivalencies_CA
      species <- unique(sppEquiv[LandR %in% speciesCommunities, ][[sppEquivCol]])
      species <- species[!species == ""]
    } else {
      species <- speciesCommunities
    }
  }

  return(list(speciesRas = bb, speciesList = species))
}


#' @keywords internal
forestOutline <- function(x) {
  x1 <- terra::vect(x)
  bf2 <- terra::simplifyGeom(x1, tolerance = 5000)
  bf3 <- terra::makeValid(bf2)
  bf4 <- terra::aggregate(bf3)
  bf5 <- terra::makeValid(bf4)
  bf6 <- terra::buffer(bf5, 6000)
  bf7 <- terra::aggregate(bf6)
  bf8 <- terra::buffer(bf7, -6000)
}

# ## TODO: randomized URL changes
# borealForestURL <- "https://d278fo2rk9arr5.cloudfront.net/downloads/boreal.zip"

#' @keywords internal
SA_ERIntersect <- function(x, studyArea) {
  x <- sf::st_read(x)
  sa_sf <- sf::st_as_sf(studyArea)
  ecoregions <- sf::st_transform(x, sf::st_crs(sa_sf))
  studyAreaER <- sf::st_intersects(ecoregions, sa_sf, sparse = FALSE)
  terra::vect(sf::as_Spatial(ecoregions[studyAreaER, ]))
}

#' @keywords internal
loadAndAggregateKNN <- function(dPath, res, sa) {
  ll <- loadkNNSpeciesLayers(dPath, sppEquiv = LandR::sppEquivalencies_CA, sppEquivCol = "KNN")
  llCoarse <- terra::aggregate(ll, res / 250)
  postProcess(llCoarse, cropTo = sa, maskTo = sa, method = "near")
}
