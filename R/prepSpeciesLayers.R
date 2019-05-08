if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("AGE", "CC", "GID", "keepSpecies", "layerName", "pct", "value"))
}

#' Load CASFRI data
#'
#' TODO: description needed
#'
#' @param CASFRIRas TODO: description needed
#' @param attrFile TODO: description needed
#' @param headerFile TODO: description needed
#' @param sppEquiv TODO: description needed
#' @param sppEquivCol TODO: description needed
#' @param type Character string. Either \code{"cover"} or \code{"age"}.
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom data.table data.table fread melt set setkey
#' @importFrom reproducible asPath Cache
loadCASFRI <- function(CASFRIRas, attrFile, headerFile, sppEquiv, sppEquivCol,
                       type = c("cover", "age")) {
  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  sppNameVectorCASFRI <- equivalentName(sppNameVector, sppEquiv,  column = "CASFRI", multi = TRUE)

  # CASFRI stuff
  CASFRIheader <- fread(headerFile, skip = 14, nrows = 49, header = FALSE, sep = "", fill = TRUE)
  header <- apply(CASFRIheader, 1, function(x) sub(pattern = "(\t+| ).*$", "", x))
  CASFRIheader <- header[nchar(header) != 0]

  wantedColumns <- grep(CASFRIheader, pattern = "^SPECIES|^GID|^AGE")

  CASFRIattr <- fread(asPath(attrFile), select = wantedColumns)

  setnames(CASFRIattr, CASFRIheader[wantedColumns])
  setkey(CASFRIattr, "GID")

  NAVals <- c("XXXX MISS", "UNDEF", "XXXX ERRC")
  numSpeciesColumns <- length(grep("SPECIES_PER", names(CASFRIattr), value = TRUE))
  if (type[1] == "cover") {
    for (i in seq(numSpeciesColumns)) {
      set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_", i)]] %in% NAVals),
          paste0("SPECIES_", i), NA_character_)
      set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_", i)]] %in% NAVals),
          paste0("SPECIES_", i), NA_character_)
    }
    for (i in 1:1) {
      message("remove CASFRI entries with <15 cover as dominant species,",
              " i.e., these pixels are deemed untreed")
      CASFRIattr <- CASFRIattr[which(CASFRIattr[[paste0("SPECIES_PER_", i)]] > 15), ]
    }
    message("set CASFRI entries with <15 cover in 2nd-5th dominance class to NA")
    for (i in 2:5) {
      set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_", i)]] <= 15),
          paste0("SPECIES_", i), NA_character_)
    }

    CASFRIattrLong <- melt(CASFRIattr, id.vars = c("GID"),
                           measure.vars = paste0("SPECIES_", 1:5))
    CA2 <- melt(CASFRIattr, id.vars = c("GID"),
                measure.vars = c(paste0("SPECIES_PER_", 1:5)))
    CASFRIattrLong[, pct := CA2$value]
    rm(CA2)
    CASFRIattrLong <- na.omit(CASFRIattrLong)
    CASFRIattrLong <- CASFRIattrLong[value %in% sppNameVectorCASFRI]
  } else {
    CASFRIattrLong <- CASFRIattr[, .(GID, AGE)]
    CASFRIattrLong <- CASFRIattrLong[!is.na(AGE) & AGE > -1]
  }

  CASFRIdt <- data.table(GID = CASFRIRas[], rastInd = 1:ncell(CASFRIRas))
  CASFRIdt <- CASFRIdt[!is.na(GID)]
  #CASFRIdt <- CASFRIdt[isNA == FALSE]
  setkey(CASFRIdt, GID)
  #set(CASFRIdt, NULL, "isNA", NULL)

  return(list(#keepSpecies = keepSpecies,
    CASFRIattrLong = CASFRIattrLong,
    CASFRIdt = CASFRIdt))
}

#' CASFRItoSpRasts
#'
#' TODO: description and title needed
#'
#' @param CASFRIRas TODO: description needed
#' @param CASFRIattrLong TODO: description needed
#' @param CASFRIdt TODO: description needed
#' @param sppEquiv TODO: description needed
#' @param sppEquivCol TODO: description needed
#' @param destinationPath TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom data.table setkey
#' @importFrom magrittr %>%
#' @importFrom reproducible asPath Cache
#' @importFrom raster crs crs<- raster setValues stack writeRaster
CASFRItoSpRasts <- function(CASFRIRas, CASFRIattrLong, CASFRIdt,
                            sppEquiv, sppEquivCol, destinationPath) {
  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This
  sppListMergesCASFRI <-lapply(sppNameVector, function(x)
    equivalentName(x, sppEquiv,  column = "CASFRI", multi = TRUE)
  )

  ## create list and template raster
  spRasts <- list()
  spRas <- raster(CASFRIRas) %>% setValues(., NA_integer_)

  ## NOT SURE IF THESE LINES ABOUT NA are relevant -- Eliot Dec 7
  ## selected spp absent from CASFRI data
  NA_Sp <- which(is.na(sppListMergesCASFRI))#setdiff(speciesLandR, unique(keepSpecies$spGroup))

  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp))
    warning(cat("Not all selected species are in loadedCASFRI. Check if this is correct:\n",
                paste(paste0(keepSpecies$CASFRI[NA_Sp], collapse = ", "), "absent\n")))

  ## empty rasters for NA_sp
  for (sp in NA_Sp) {
    message("  running ", sp, ". Assigning NA, because absent from CASFRI")
    spRasts[[sp]] <- spRas
    spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                           filename = asPath(file.path(destinationPath,
                                                       paste0("CASFRI", sp, ".tif"))),
                           overwrite = TRUE, datatype = "INT2U")
  }

  sppTODO <- unique(names(sppListMergesCASFRI))

  for (sp in sppTODO) {
    spCASFRI <- sppListMergesCASFRI[[sp]]
    spRasts[[sp]] <- spRas
    message("starting ", sp)
    if (length(spCASFRI) > 1)
      message("  Merging ", paste(spCASFRI, collapse = ", "), "; becoming: ", sp)
    aa2 <- CASFRIattrLong[
      value %in% spCASFRI][
        , min(100L, sum(pct)), by = GID]
    setkey(aa2, GID)
    cc <- aa2[CASFRIdt] %>% na.omit()
    rm(aa2)
    spRasts[[sp]][cc$rastInd] <- cc$V1
    message("  ", sp, " writing to disk")

    startCRS <- crs(spRasts[[sp]])
    spRasts[[sp]] <- writeRaster(spRasts[[sp]],
                                 filename = asPath(file.path(destinationPath,
                                                             paste0("CASFRI", sp, ".tif"))),
                                 datatype = "INT1U", overwrite = TRUE)

    if (is(spRasts[[sp]], "Raster")) {
      # Rasters need to have their disk-backed value assigned, but not shapefiles
      # This is a bug in writeRaster was spotted with crs of rastTmp became
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
      # should have stayed at
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0
      if (!identical(startCRS, crs(spRasts[[sp]])))
        crs(spRasts[[sp]]) <- startCRS
    }
    message("  ", sp, " done")
  }

  raster::stack(spRasts)
}

#' Prepare species layers
#'
#' TODO: description needed
#'
#' @param destinationPath TODO: description needed
#' @param outputPath TODO: description needed
#' @param url TODO: description needed; if \code{NULL}, the default, use the default source url
#' @param studyArea TODO: description needed
#' @param rasterToMatch TODO: description needed
#' @param sppEquiv TODO: description needed
#' @param sppEquivCol TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom reproducible asPath Cache prepInputs
#' @rdname prepSpeciesLayers
prepSpeciesLayers_CASFRI <- function(destinationPath, outputPath,
                                     url = NULL,
                                     studyArea, rasterToMatch,
                                     sppEquiv,
                                     sppEquivCol) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h/view?usp=sharing"

  CASFRItiffFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs.tif"))
  CASFRIattrFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs_attributes3.csv"))
  CASFRIheaderFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs_README.txt"))

  message("  Loading CASFRI layers...")
  CASFRIRas <- Cache(prepInputs,
                     #targetFile = asPath("Landweb_CASFRI_GIDs.tif"),
                     targetFile = basename(CASFRItiffFile),
                     archive = asPath("CASFRI for Landweb.zip"),
                     url = url,
                     alsoExtract = c(CASFRItiffFile, CASFRIattrFile, CASFRIheaderFile),
                     destinationPath = destinationPath,
                     fun = "raster::raster",
                     studyArea = studyArea,
                     rasterToMatch = rasterToMatch,
                     method = "bilinear", ## ignore warning re: ngb (#5)
                     datatype = "INT4U",
                     filename2 = NULL, #TRUE,
                     overwrite = TRUE,
                     userTags =  c("CASFRIRas", "stable"))

  message("Load CASFRI data and headers, and convert to long format, and define species groups")

  #Cache
  loadedCASFRI <- Cache(loadCASFRI,
                        CASFRIRas = CASFRIRas,
                        attrFile = CASFRIattrFile,
                        headerFile = CASFRIheaderFile, ## TODO: this isn't used internally
                        sppEquiv = sppEquiv,
                        sppEquivCol = sppEquivCol,
                        type = "cover"#,
                        #userTags = c("function:loadCASFRI", "BigDataTable",
                        #"speciesLayers", "KNN")
  )

  message("Make stack from CASFRI data and headers")
  CASFRISpStack <- CASFRItoSpRasts(CASFRIRas = CASFRIRas,
                                   sppEquiv = sppEquiv,
                                   sppEquivCol = sppEquivCol,
                                   CASFRIattrLong = loadedCASFRI$CASFRIattrLong,
                                   CASFRIdt = loadedCASFRI$CASFRIdt,
                                   destinationPath = outputPath)

  return(CASFRISpStack)
}

#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_KNN <- function(destinationPath, outputPath,
                                  url = NULL,
                                  studyArea, rasterToMatch,
                                  sppEquiv,
                                  sppEquivCol) {
  if (is.null(url))
    url <- "http://tree.pfc.forestry.ca/kNN-Species.tar"

  loadkNNSpeciesLayers(
    dPath = destinationPath,
    rasterToMatch = rasterToMatch,
    studyArea = studyArea,
    sppEquiv = sppEquiv,
    knnNamesCol = "KNN",
    sppEquivCol = sppEquivCol,
    thresh = 10,
    url = url,
    userTags = c("speciesLayers", "KNN"))
}

#' makePickellStack
#'
#' TODO: description and title needed
#'
#' @param PickellRaster TODO: description needed
#' @param sppEquiv TODO: description needed
#' @param sppEquivCol TODO: description needed
#' @param destinationPath TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom reproducible asPath Cache
#' @importFrom raster raster rasterOptions setValues stack
makePickellStack <- function(PickellRaster, sppEquiv, sppEquivCol, destinationPath) {
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the speciesEquivalency table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  PickellSpp <- c("Pice_mar", "Pice_gla", "Pinu_con", "Popu_tre") ## pick LandR as standard
  names(PickellSpp) <- PickellSpp

  # Pick the full LandR dataset, which should be broad. We will change to sppEquivCol below
  sppOfInterest <- equivalentName(sppNameVector, sppEquiv, "LandR", multi = TRUE)
  sppInPickell <- lapply(PickellSpp, function(sp)
    equivalentName(sp, sppEquiv, "LandR", multi = TRUE)
  )

  # Check that each of the layers that Pickell did are actually desired in speciesEquivalency
  needPickell <- vapply(sppInPickell, function(sp) {
    any(sp %in% sppOfInterest)
  }, logical(1))

  # These are the ones in Pickell data set that we want according to speciesEquivalency
  PickellSpp <- equivalentName(PickellSpp[needPickell], sppEquiv, sppEquivCol)

  ## bring to memory and replace water, non veg by NAs
  PickellRaster[] <- PickellRaster[]
  PickellRaster[PickellRaster[] %in% c(230, 220, 255)] <- NA_integer_

  ## create list and template raster
  spRasts <- list()
  spRas <- raster(PickellRaster) %>% setValues(., NA_integer_)

  rasterOptions(maxmemory = 1e9)

  ## converting existing species codes into percentages
  for (sp in PickellSpp) {
    message("  converting Pickell's codes to pct cover raster, for ", sp)

    if (!is.na(equivalentName("Pice_gla", sppEquiv, sppEquivCol))) {
      if (sp == equivalentName("Pice_gla", sppEquiv, sppEquivCol)) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(41, 42, 43)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(44)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(14, 34)] <- 40
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U")
      }
    }

    if (!is.na(equivalentName("Pice_mar", sppEquiv, sppEquivCol))) {
      if (sp == equivalentName("Pice_mar", sppEquiv, sppEquivCol)) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(23, 26)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(22)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(32, 42)] <- 40
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U")
      }
    }

    if (any(!is.na(equivalentName("Pinu_ban", sppEquiv, sppEquivCol)),
            !is.na(equivalentName("Pinu_con", sppEquiv, sppEquivCol)),
            !is.na(equivalentName("Pinu_spp", sppEquiv, sppEquivCol)))) {
      if (sp %in% c(equivalentName("Pinu_ban", sppEquiv, sppEquivCol),
                    equivalentName("Pinu_con", sppEquiv, sppEquivCol),
                    equivalentName("Pinu_sp", sppEquiv, sppEquivCol))) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(31, 32, 34)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(33)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(23, 43)] <- 40
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U")
      }
    }

    if (!is.na(equivalentName("Popu_tre", sppEquiv, sppEquivCol))) {
      if (sp == equivalentName("Popu_tre", sppEquiv, sppEquivCol)) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(14)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(11)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(31, 41)] <- 40
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT2U")
      }
    }
  }

  ## species in Pickell's data
  raster::stack(spRasts)
}

#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_Pickell <- function(destinationPath, outputPath,
                                      url = NULL,
                                      studyArea, rasterToMatch,
                                      sppEquiv,
                                      sppEquivCol) {

  if (is.null(url))
    url <- "https://drive.google.com/open?id=1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4"

  speciesLayers <- Cache(prepInputs,
                         targetFile = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.dat"),
                         url = url,
                         archive = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.zip"),
                         alsoExtract = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.hdr"),
                         destinationPath = destinationPath,
                         fun = "raster::raster",
                         studyArea = studyArea,
                         rasterToMatch = rasterToMatch,
                         method = "bilinear", ## ignore warning re: ngb (#5)
                         datatype = "INT2U",
                         filename2 = NULL,
                         overwrite = TRUE,
                         userTags = c("speciesLayers", "KNN", "Pickell", "stable"))

  makePickellStack(PickellRaster = speciesLayers,
                   sppEquiv = sppEquiv,
                   sppEquivCol = sppEquivCol,
                   destinationPath = destinationPath)
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom map mapAdd maps
#' @importFrom raster maxValue minValue stack
#' @rdname prepSpeciesLayers
prepSpeciesLayers_ForestInventory <- function(destinationPath, outputPath,
                                              url = NULL,
                                              studyArea, rasterToMatch,
                                              sppEquiv,
                                              sppEquivCol) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1JnKeXrw0U9LmrZpixCDooIm62qiv4_G1/view?usp=sharing"

  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This includes LandType because it will use that at the bottom of this function to
  #  remove NAs
  CClayerNames <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "LandType")
  CClayerNamesFiles <- paste0(gsub(" ", "", CClayerNames), "1.tif")
  options(map.useParallel = FALSE) ## TODO: pass additional arg to function
  ml <- mapAdd(rasterToMatch, isRasterToMatch = TRUE, layerName = "rasterToMatch",
               #useSAcrs = TRUE, #poly = TRUE,
               #      columnNameForLabels = "NSN",
               filename2 = NULL)

  ml <- mapAdd(studyArea, map = ml, isStudyArea = TRUE, layerName = "studyArea",
               useSAcrs = TRUE, #poly = TRUE,
               #      columnNameForLabels = "NSN",
               filename2 = NULL)

  ml <- mapAdd(map = ml, url = url, layerName = CClayerNames, CC = TRUE,
               destinationPath = destinationPath,
               targetFile = CClayerNamesFiles, filename2 = NULL,
               alsoExtract = "similar", leaflet = FALSE, method = "ngb")

  ccs <- ml@metadata[CC == TRUE & !(layerName == "LandType"), ]
  CCs <- maps(ml, layerName = ccs$layerName)
  CCstack <- raster::stack(CCs)
  CCstackNames <- names(CCstack)

  assertthat::assert_that(identical(raster::minValue(CCstack), rep(0L, nlayers(CCstack))))
  assertthat::assert_that(identical(raster::maxValue(CCstack), rep(10L, nlayers(CCstack))))

  CCstack[CCstack[] < 0] <- 0  ## turns stack into brick, so need to restack later
  CCstack[CCstack[] > 10] <- 10
  CCstack <- CCstack * 10 # convert back to percent
  NA_ids <- which(is.na(ml$LandType[]) |  # outside of studyArea polygon
                    ml$LandType[] == 1)   # 1 is cities -- NA it here -- will be filled in with another veg layer if available (e.g. Pickell)
  message("  Setting NA, 1 in LandType to NA in speciesLayers in ForestInventory data")
  CCstack[NA_ids] <- NA
  names(CCstack) <- equivalentName(CCstackNames, sppEquiv, sppEquivCol)

  stack(CCstack)
}
