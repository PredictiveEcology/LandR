#' Prepare species layers
#'
#'  TODO: description needed
#'
#' @param destinationPath TODO: description needed
#' @param outputPath TODO: description needed
#' @param url TODO: description needed; if \code{NULL}, the default, use the default source url
#' @param studyArea TODO: description needed
#' @param rasterToMatch TODO: description needed
#' @param sppEquiv TODO: description needed
#' @param sppEquivCol TODO: description needed
#'
#' @return
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
  CASFRIheaderFile <- asPath(file.path(destinationPath,"Landweb_CASFRI_GIDs_README.txt"))

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

  message('Make stack from CASFRI data and headers')
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
#' @importFrom map mapAdd maps
#' @importFrom raster stack
#' @rdname prepSpeciesLayers
prepSpeciesLayers_ForestInventory <- function(destinationPath, outputPath,
                                              url = NULL,
                                              studyArea, rasterToMatch,
                                              sppEquiv,
                                              sppEquivCol) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1JnKeXrw0U9LmrZpixCDooIm62qiv4_G1/view?usp=sharing"

  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]),]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This
  sppListMergesCASFRI <-lapply(sppNameVector, function(x)
    equivalentName(x, sppEquiv,  column = "CASFRI", multi = TRUE)
  )

  # This includes LandType because it will use that at the bottom of this function to
  #  remove NAs
  CClayerNames <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "LandType")
  CClayerNamesWDots <- gsub(" ", ".", CClayerNames)
  CClayerNamesLandR <- equivalentName(CClayerNamesWDots, sppEquiv, sppEquivCol, multi = TRUE)
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
  CCstack[CCstack[] < 0] <- 0
  CCstack[CCstack[] > 10] <- 10
  CCstack <- CCstack * 10 # convert back to percent
  NA_ids <- which(is.na(ml$LandType[]) |  # outside of studyArea polygon
                    ml$LandType[] == 1)   # 1 is cities -- NA it here -- will be filled in with another veg layer if available (e.g. Pickell)
  message("  Setting NA, 1 in LandType to NA in speciesLayers in ForestInventory data")
  CCstack[NA_ids] <- NA

  names(CCstack) <- equivalentName(names(CCstack), sppEquiv, sppEquivCol)

  CCstack
}
