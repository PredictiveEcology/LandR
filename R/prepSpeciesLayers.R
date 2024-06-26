utils::globalVariables(c(
  "AGE", "CC", "GID", "id", "keepSpecies", "layerName", "name", "pct", "value"
))

#' Load CASFRI data
#'
#' TODO: description needed
#'
#' @param CASFRIRas TODO: description needed
#' @param attrFile TODO: description needed
#' @param headerFile TODO: description needed
#' @template sppEquiv
#' @template sppEquivCol
#' @param type Character string. Either `"cover"` or `"age"`.
#'
#' @return TODO: description needed
#'
#' @export
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

  CASFRIdt <- data.table(GID = as.vector(values(CASFRIRas)), rastInd = 1:ncell(CASFRIRas))
  CASFRIdt <- CASFRIdt[!is.na(GID)]
  #CASFRIdt <- CASFRIdt[isNA == FALSE]
  setkey(CASFRIdt, GID)
  #set(CASFRIdt, NULL, "isNA", NULL)

  return(list(CASFRIattrLong = CASFRIattrLong, CASFRIdt = CASFRIdt))
}

#' `CASFRItoSpRasts`
#'
#' TODO: description and title needed
#'
#' @param CASFRIRas TODO: description needed
#' @param CASFRIattrLong TODO: description needed
#' @param CASFRIdt TODO: description needed
#' @template sppEquiv
#' @template sppEquivCol
#' @template destinationPath
#'
#' @return TODO: description needed
#'
#' @export
CASFRItoSpRasts <- function(CASFRIRas, CASFRIattrLong, CASFRIdt,
                            sppEquiv, sppEquivCol, destinationPath) {
  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This
  sppListMergesCASFRI <- lapply(sppNameVector, function(x)
    equivalentName(x, sppEquiv,  column = "CASFRI", multi = TRUE)
  )

  ## create list and template raster
  spRasts <- list()
  spRas <- CASFRIRas
  spRas[] <- NA_integer_

  ## NOT SURE IF THESE LINES ABOUT NA are relevant -- Eliot Dec 7
  ## selected spp absent from CASFRI data
  NA_Sp <- which(is.na(sppListMergesCASFRI))#setdiff(speciesLandR, unique(keepSpecies$spGroup))

  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp))
    warning("Not all selected species are in loadedCASFRI. Check if this is correct:\n",
            paste(paste0(keepSpecies$CASFRI[NA_Sp], collapse = ", "), "absent\n"))

  ## empty rasters for NA_sp
  for (sp in NA_Sp) {
    message("  running ", sp, ". Assigning NA, because absent from CASFRI")
    spRasts[[sp]] <- spRas
    NAval <- 65535L
    spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                           filename = asPath(file.path(destinationPath,
                                                       paste0("CASFRI_", sp, ".tif"))),
                           overwrite = TRUE, datatype = "INT2U", NAflag = NAval)
    ## NAvals need to be converted back to NAs
    spRasts[[sp]] <- .NAvalueFlag(spRasts[[sp]], NAval)

  }

  sppTODO <- unique(names(sppListMergesCASFRI))

  for (sp in sppTODO) {
    spCASFRI <- sppListMergesCASFRI[[sp]]
    spRasts[[sp]] <- spRas
    message("starting ", sp)
    if (length(spCASFRI) > 1)
      message("  Merging ", paste(spCASFRI, collapse = ", "), "; becoming: ", sp)
    aa2 <- CASFRIattrLong[value %in% spCASFRI][, min(100L, sum(pct)), by = GID]
    setkey(aa2, GID)
    cc <- aa2[CASFRIdt] |>
      na.omit()
    rm(aa2)
    spRasts[[sp]][cc$rastInd] <- cc$V1
    message("  ", sp, " writing to disk")

    startCRS <- crs(spRasts[[sp]])
    NAval <- 255L
    spRasts[[sp]] <- writeRaster(spRasts[[sp]],
                                 filename = asPath(file.path(destinationPath,
                                                             paste0("CASFRI_", sp, ".tif"))),
                                 datatype = "INT1U", overwrite = TRUE, NAflag = NAval)
    ## NAvals need to be converted back to NAs
    spRasts[[sp]] <- .NAvalueFlag(spRasts[[sp]], NAval)

    if (is(spRasts[[sp]], "Raster") || is(spRasts[[sp]], "SpatRaster")) {
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

  ## if the original raster was `RasterLayer`, those in the list will be also.
  spRasts <- .stack(spRasts)

  return(spRasts)
}

#' Prepare species layers
#'
#' TODO: description needed
#'
#' @template destinationPath
#' @param outputPath TODO: description needed
#' @param url if `NULL`, the default, use the default source url
#' @template studyArea
#' @template rasterToMatch
#' @template sppEquiv
#' @template sppEquivCol
#' @param thresh threshold \% cover used to defined the species as "present" in the study area.
#'    If at least one pixel has `cover >= thresh` , the species is considered "present".
#'    Otherwise the raster is excluded from the output. Defaults to 10.
#' @param ... other arguments, used for compatibility with other `prepSpeciesLayers` functions.
#'
#' @return TODO: description needed
#'
#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_KNN <- function(destinationPath, outputPath,
                                  url = NULL,
                                  studyArea, rasterToMatch,
                                  sppEquiv,
                                  sppEquivCol,
                                  thresh = 10, ...) {
  stopifnot(requireNamespace("RCurl", quietly = TRUE))

  dots <- list(...)

  if ("year" %in% names(dots)) {
    year <- dots[["year"]]
  } else {
    year <- 2001
  }

  if (is.null(url)) {
    url <- paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                  "canada-forests-attributes_attributs-forests-canada/", year,
                  "-attributes_attributs-", year, "/")
  }

  shared_drive_url <- NULL
  if (!RCurl::url.exists(url)) {  ## ping website and use gdrive if not available
    if (requireNamespace("googledrive", quietly = TRUE)) {
      driveFolder <- paste0("kNNForestAttributes_", year)
      shared_drive_url <- "https://drive.google.com/drive/folders/0AJE09VklbHOuUk9PVA"
      # url <- googledrive::with_drive_quiet(
      #   googledrive::drive_link(
      #     googledrive::drive_ls(
      #       driveFolder,
      #       shared_drive = googledrive::as_id(shared_drive_url)
      #     )
      #   )
      # )

      driveDT <- as.data.table(googledrive::drive_ls(googledrive::as_id(shared_drive_url)))
      url <- googledrive::with_drive_quiet(
        googledrive::drive_link(driveDT[name == driveFolder, id])
      )
    }
  }

  loadkNNSpeciesLayers(
    dPath = destinationPath,
    knnNamesCol = "KNN",
    outputPath = outputPath,
    rasterToMatch = rasterToMatch,
    studyArea = studyArea,
    studyAreaName = dots$studyAreaName,
    sppEquiv = sppEquiv,
    sppEquivCol = sppEquivCol,
    thresh = thresh,
    url = url,
    year = year,
    shared_drive_url = shared_drive_url,
    userTags = c("speciesLayers", "KNN")
  )
}

#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_CASFRI <- function(destinationPath, outputPath,
                                     url = NULL,
                                     studyArea, rasterToMatch,
                                     sppEquiv,
                                     sppEquivCol, ...) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h"

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
                     studyArea = studyArea,
                     rasterToMatch = rasterToMatch,
                     method = "bilinear", ## ignore warning re: ngb (#5)
                     datatype = "INT4U",
                     filename2 = NULL,
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
prepSpeciesLayers_Pickell <- function(destinationPath, outputPath,
                                      url = NULL,
                                      studyArea, rasterToMatch,
                                      sppEquiv,
                                      sppEquivCol, ...) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4"

  speciesLayers <- Cache(prepInputs,
                         targetFile = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.dat"),
                         url = url,
                         archive = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.zip"),
                         alsoExtract = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.hdr"),
                         destinationPath = destinationPath,
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
                   destinationPath = outputPath)
}

#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_ForestInventory <- function(destinationPath, outputPath,
                                              url = NULL,
                                              studyArea, rasterToMatch,
                                              sppEquiv,
                                              sppEquivCol, ...) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1JnKeXrw0U9LmrZpixCDooIm62qiv4_G1"

  ## TODO: add terra compatible methods.

  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This includes LandType because it will use that at the bottom of this function to
  #  remove NAs
  CClayerNames <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "LandType")
  CClayerNamesFiles <- paste0(gsub(" ", "", CClayerNames), "1.tif")

  lr <- lapply(CClayerNamesFiles, prepInputs, studyArea = studyArea, rasterToMatch = rasterToMatch,
               url = url, alsoExtract = "similar", method = "ngb",
               destinationPath = destinationPath, filename2 = NULL)
  rs <- raster::stack(lr)
  names(rs) <- CClayerNames

  CCstack <- dropLayer(rs, which(grepl("LandType", CClayerNames)))
  CCstackNames <- names(CCstack)

  if (!all(min(CCstack[], na.rm = TRUE) >= 0)) stop("problem with min. of CCstack (< 0)")
  if (!all(max(CCstack[], na.rm = TRUE) <= 10)) stop("problem with max. of CCstack (> 10)")

  CCstack <- CCstack * 10 ## convert back to percent
  ## NA means outside of studyArea polygon; 1 is cities, Set to NA here:
  ## will be filled in with another veg layer if available (e.g. Pickell)
  NA_ids <- which(is.na(rs$LandType[]) | rs$LandType[] == 1)
  message("  Setting NA, 1 in LandType to NA in speciesLayers in ForestInventory data")
  aa <- try(CCstack[NA_ids] <- NA, silent = TRUE)
  ## unclear why line above sometimes fails: 'Error in value[j, ] : incorrect number of dimensions'
  if (is(aa, "try-error")) {
    l <- unstack(CCstack)
    CCstack <- lapply(l, function(x) {
      x[NA_ids] <- NA
      x
    })
  }
  names(CCstack) <- equivalentName(CCstackNames, sppEquiv, sppEquivCol)

  stack(CCstack)
}

#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_MBFRI <- function(destinationPath, outputPath,
                                    url = NULL,
                                    studyArea, rasterToMatch,
                                    sppEquiv,
                                    sppEquivCol, ...) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1KTqNBntNrEsDL6jk-5bchsBOcraDqNHe"

  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  ## NOTE: This includes LandType because we will use it to remove NAs
  CClayerNames <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "Landtype") ## file is 'Landtype'
  CClayerNames2 <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "LandType") ## needs 'LandType'
  CClayerNamesFiles <- paste0("MB_", gsub(" ", "", CClayerNames), "2016_NRV.tif")

  lr <- lapply(CClayerNamesFiles, prepInputs, studyArea = studyArea, rasterToMatch = rasterToMatch,
               url = url, alsoExtract = "similar", method = "ngb",
               destinationPath = destinationPath, filename2 = NULL)
  rs <- stack(lr)
  names(rs) <- CClayerNames2

  CCstack <- dropLayer(rs, which(grepl("LandType", CClayerNames2)))

  if (!all(min(CCstack[], na.rm = TRUE) >= 0)) stop("problem with min. of CCstack (< 0)")
  if (!all(max(CCstack[], na.rm = TRUE) <= 10)) stop("problem with max. of CCstack (> 10)")

  CCstack <- CCstack * 10 # convert back to percent
  ## NA means outside of studyArea polygon; 1 is cities, Set to NA here:
  ## will be filled in with another veg layer if available (e.g. Pickell)
  NA_ids <- which(is.na(rs$LandType[]) | rs$LandType[] == 1)
  message("  Setting NA, 1 in LandType to NA in speciesLayers in MB FRI data")
  aa <- try(CCstack[NA_ids] <- NA, silent = TRUE)
  ## unclear why line above sometimes fails: 'Error in value[j, ] : incorrect number of dimensions'
  if (is(aa, "try-error")) {
    l <- unstack(CCstack)
    CCstack <- lapply(l, function(x) {
      x[NA_ids] <- NA
      x
    })
  }
  names(CCstack) <- equivalentName(CClayerNames2, sppEquiv, sppEquivCol)

  stack(CCstack)
}

#' @export
#' @rdname prepSpeciesLayers
prepSpeciesLayers_ONFRI <- function(destinationPath, outputPath,
                                    url = NULL,
                                    studyArea, rasterToMatch,
                                    sppEquiv,
                                    sppEquivCol, ...) {

  ## TODO: this is sneaky and annoying (studyAreaName is part of outputPath)
  if (grepl("AOU", dirname(outputPath))) {
    sA <- "ceon"
    sAN <- "AOU"
  } else if (grepl("ROF", dirname(outputPath))) {
    sA <- "rof"
    sAN <- "ROF"
  }

  if (grepl("res125", dirname(outputPath))) {
    res <- 125L
  } else {
    res <- 250L
  }

  if (is.null(url)) {
    if (sA == "ceon") {
      url <- "https://drive.google.com/file/d/1iJq1wv06FYTvnkBXKr5HKa84ZF7QbalS"
      url2 <- "https://drive.google.com/file/d/1eg9yhkAKDsQ8VO5Nx4QjBg4yiB0qyqng"
    } else if (sA == "rof") {
      if (res == 125) {
        url <- "https://drive.google.com/file/d/12C2a3GpnwIz5j2EFERZQxCFYyZ2tWqFm"
        url2 <- "https://drive.google.com/file/d/1JouBj0iJOPB1qQeXkRRePMN6MZSX_R_q"
      } else if (res == 250) {
        url <- "https://drive.google.com/file/d/1MCQhlzwVc7KhNNPfn7uK5BlI6OgUarFS"
        url2 <- "https://drive.google.com/file/d/1-2XSrSp_WrZCnqUhHTaj0rQpzOcSLrfS"
      }
    }
  }

  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  if (is.null(sppEquiv[["ONFRI"]]))
    stop("Column 'ONFRI' not found in species equivalent table (sppEquiv).")

  FRIlayerNames <- unique(sppEquiv[["ONFRI"]])
  FRIlayerNamesFiles <- paste0(FRIlayerNames, "_fri_", sA, "_", res, "m.tif")
  FRIlccname <- paste0("lcc_fri_", sA, "_", res, "m.tif")

  sppLayers <- rast(lapply(FRIlayerNamesFiles, function(f) {
    prepInputs(url = url, studyArea = studyArea, rasterToMatch = rasterToMatch,
               destinationPath = destinationPath, targetFile = f, filename2 = NULL,
               alsoExtract = NA, method = "near")
  }))
  names(sppLayers) <- FRIlayerNames

  if (!all(minmax(sppLayers)[1, ] >= 0)) stop("problem with min. of species layers stack (< 0)")
  if (!all(minmax(sppLayers)[2, ] <= 100)) stop("problem with max. of species layers stack (> 100)")

  ## merge species layers (currently only Popu; TODO: pine?)
  idsPopu <- grep("Popu", FRIlayerNames)
  mergedPopu <- app(c(sppLayers[[idsPopu]]), sum, na.rm = TRUE)
  sppLayers <- sppLayers[[!names(sppLayers) %in% names(sppLayers)[idsPopu]]]
  sppLayers[["Popu_sp"]] <- mergedPopu ## NOTE: addLayer sporadically fails to add layer, w/o warning

  idThuj <- grep("Thuj_spp", names(sppLayers))
  names(sppLayers[[idThuj]]) <- "Thuj_sp"

  sppLayers
}

#' @template destinationPath
#' @param outputPath path to output directory
#' @export
#' @rdname LandR-deprecated
prepSpeciesLayers_KNN2011 <- function(destinationPath, outputPath, url = NULL, studyArea,
                                      rasterToMatch, sppEquiv, sppEquivCol, thresh = 10, ...) {
  .Deprecated("loadkNNSpeciesLayers",
              msg = paste("prepSpeciesLayers_KNN2011 is deprecated.",
                          "Please use 'loadkNNSpeciesLayers' and supply URL/year to validation layers."))

  loadkNNSpeciesLayers(
    dPath = destinationPath,
    rasterToMatch = rasterToMatch,
    studyArea = studyArea,
    sppEquiv = sppEquiv,
    year = 2011,
    knnNamesCol = "KNN",
    sppEquivCol = sppEquivCol,
    thresh = thresh,
    url = url,
    ## dots start here:
    outputPath = outputPath,
    userTags = c("speciesLayers", "KNN"),
    ...
  )
}

#' `makePickellStack`
#'
#' TODO: description and title needed
#'
#' @param PickellRaster TODO: description needed
#' @template sppEquiv
#' @template sppEquivCol
#' @template destinationPath
#'
#' @return TODO: description needed
#'
#' @export
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
  PickellRaster[] <- as.vector(PickellRaster[])
  PickellRaster[PickellRaster[] %in% c(230, 220, 255)] <- NA_integer_

  ## create list and template raster
  spRasts <- list()
  spRas <- PickellRaster
  spRas[] <- NA_integer_

  rasterOptions(maxmemory = 1e9)
  terraOptions(memmax = 1L)

  ## converting existing species codes into percentages
  for (sp in PickellSpp) {
    message("  converting Pickell's codes to pct cover raster, for ", sp)

    if (!is.na(equivalentName("Pice_gla", sppEquiv, sppEquivCol))) {
      if (sp == equivalentName("Pice_gla", sppEquiv, sppEquivCol)) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(41, 42, 43)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(44)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(14, 34)] <- 40
        NAval <- 255L
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell_", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U", NAflag = NAval)
        ## NAvals need to be converted back to NAs
        spRasts[[sp]] <- .NAvalueFlag(spRasts[[sp]], NAval)
      }
    }

    if (!is.na(equivalentName("Pice_mar", sppEquiv, sppEquivCol))) {
      if (sp == equivalentName("Pice_mar", sppEquiv, sppEquivCol)) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(23, 26)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(22)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(32, 42)] <- 40

        NAval <- 255L
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell_", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U", NAflag = NAval)
        ## NAvals need to be converted back to NAs
        spRasts[[sp]] <- .NAvalueFlag(spRasts[[sp]], NAval)
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

        NAval <- 255L
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell_", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U", NAflag = NAval)
        ## NAvals need to be converted back to NAs
        spRasts[[sp]] <- .NAvalueFlag(spRasts[[sp]], NAval)
      }
    }

    if (!is.na(equivalentName("Popu_tre", sppEquiv, sppEquivCol))) {
      if (sp == equivalentName("Popu_tre", sppEquiv, sppEquivCol)) {
        spRasts[[sp]] <- spRas
        spRasts[[sp]][PickellRaster[] %in% c(14)] <- 60
        spRasts[[sp]][PickellRaster[] %in% c(11)] <- 80
        spRasts[[sp]][PickellRaster[] %in% c(31, 41)] <- 40

        NAval <- 65535L
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell_", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT2U", NAflag = NAval)
        ## NAvals need to be converted back to NAs
        spRasts[[sp]] <- .NAvalueFlag(spRasts[[sp]], NAval)
      }
    }
  }

  ## species in Pickell's data
  .stack(spRasts)
}


#' Convert `NA` values in `speciesLayers` to zeros
#'
#' Pixels that are `NA` but are inside `rasterToMatch` may need to be converted to `0`,
#' as they can could still potentially be forested
#'
#' @template speciesLayers
#' @template rasterToMatch
#'
#' @return the `speciesLayers` with `0` in pixels that had `NA`
#'
#' @export
NAcover2zero <- function(speciesLayers, rasterToMatch) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop("Suggested package 'terra' not installed.\n",
         "Install it using `install.packages('terra')`.")
  }

  tempRas <- rasterToMatch
  tempRas[!is.na(as.vector(tempRas[]))] <- 0
  namesLayers <- names(speciesLayers)

  message("...making sure empty pixels inside study area have 0 cover, instead of NAs ...")
  # Changed to terra Nov 17 by Eliot --> this was many minutes with raster::cover --> 3 seconds with terra
  if (is(rasterToMatch, "Raster")) {
    speciesLayers <- terra::rast(speciesLayers)
    tempRas <- terra::rast(tempRas)
  }
  speciesLayers <- terra::cover(speciesLayers, tempRas)

  if (is(rasterToMatch, "Raster")) {
    speciesLayers <- raster::stack(speciesLayers)
  }

  names(speciesLayers) <- namesLayers
  message("   ...done")

  return(speciesLayers)
}
