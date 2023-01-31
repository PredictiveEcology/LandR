utils::globalVariables(c(
  "CC", "id", "layerName", "name", "value"
))

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
#' @importFrom reproducible asPath Cache prepInputs
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
                   destinationPath = outputPath)
}

#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom map mapAdd maps
#' @importFrom raster crop extend maxValue minValue origin origin<- stack unstack
#' @rdname prepSpeciesLayers
prepSpeciesLayers_ForestInventory <- function(destinationPath, outputPath,
                                              url = NULL,
                                              studyArea, rasterToMatch,
                                              sppEquiv,
                                              sppEquivCol, ...) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1JnKeXrw0U9LmrZpixCDooIm62qiv4_G1"

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
               filename2 = NULL)

  ml <- mapAdd(studyArea, map = ml, isStudyArea = TRUE, layerName = "studyArea",
               useSAcrs = TRUE, filename2 = NULL)

  ## TODO: avoid this workaround for stack not liking rasters with different origins and extents
  lr <- lapply(CClayerNamesFiles, prepInputs, url = url, alsoExtract = "similar",
               destinationPath = destinationPath, fun = "raster::raster")
  le <- lapply(lr, extent)
  dts <- lapply(lr, function(r) {
    e <- extent(r)
    data.table(xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4])
  })
  dt <- rbindlist(dts)
  new_ext <- extent(max(dt$xmin), min(dt$xmax), max(dt$ymin), min(dt$ymax))
  lr2 <- lapply(lr, `origin<-`, value = origin(lr[[1]]))
  lr3 <- lapply(lr2, function(r) {
    crop(r, y = new_ext, filename = tempfile(fileext = ".tif"))
  })
  lr4 <- lapply(lr3, function(r) {
    extend(r, y = new_ext, filename = tempfile(fileext = ".tif"))
  })
  ## stack(lr4) ## confirm it works

  ## loop/apply doesn't work here???
  ml <- mapAdd(lr4[[1]], map = ml, layerName = CClayerNames[1], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(lr4[[2]], map = ml, layerName = CClayerNames[2], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(lr4[[3]], map = ml, layerName = CClayerNames[3], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(lr4[[4]], map = ml, layerName = CClayerNames[4], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(lr4[[5]], map = ml, layerName = CClayerNames[5], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(lr4[[6]], map = ml, layerName = CClayerNames[6], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ## END WORKAROUND

  # ml <- mapAdd(map = ml, url = url, layerName = CClayerNames, CC = TRUE,
  #              destinationPath = destinationPath,
  #              targetFile = CClayerNamesFiles, filename2 = NULL,
  #              alsoExtract = "similar", leaflet = FALSE, method = "ngb") ## TODO: fix error due to different extent

  ccs <- ml@metadata[CC == TRUE & !(layerName == "LandType"), ]
  CCs <- maps(ml, layerName = ccs$layerName)
  CCstack <- raster::stack(CCs)
  CCstackNames <- names(CCstack)

  if (!all(raster::minValue(CCstack) >= 0)) stop("problem with minValue of CCstack (< 0)")
  if (!all(raster::maxValue(CCstack) <= 10)) stop("problem with maxValue of CCstack (> 10)")

  CCstack <- CCstack * 10 # convert back to percent
  NA_ids <- which(is.na(ml$LandType[]) |  # outside of studyArea polygon
                    ml$LandType[] == 1)   # 1 is cities -- NA it here -- will be filled in with another veg layer if available (e.g. Pickell)
  message("  Setting NA, 1 in LandType to NA in speciesLayers in ForestInventory data")
  aa <- try(CCstack[NA_ids] <- NA, silent = TRUE)
  ## unclear why line above sometimes fails: 'Error in value[j, ] : incorrect number of dimensions'
  if (is(aa, "try-error")) {
    l <- unstack(CCstack)
    CCstack <- lapply(l, function(x) {x[NA_ids] <- NA; x})
  }
  names(CCstack) <- equivalentName(CCstackNames, sppEquiv, sppEquivCol)

  stack(CCstack)
}

#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom map mapAdd maps
#' @importFrom raster crop extend maxValue minValue origin origin<- stack unstack
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

  # This includes LandType because it will use that at the bottom of this function to
  #  remove NAs
  CClayerNames <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "Landtype") ## file is 'Landtype'
  CClayerNames2 <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "LandType") ## needs 'LandType'
  CClayerNamesFiles <- paste0("MB_", gsub(" ", "", CClayerNames), "2016_NRV.tif")
  options(map.useParallel = FALSE) ## TODO: pass additional arg to function
  ml <- mapAdd(rasterToMatch, isRasterToMatch = TRUE, layerName = "rasterToMatch",
               filename2 = NULL)

  ml <- mapAdd(studyArea, map = ml, isStudyArea = TRUE, layerName = "studyArea",
               useSAcrs = TRUE, filename2 = NULL)

  lr <- lapply(CClayerNamesFiles, prepInputs, url = url, alsoExtract = "similar",
               destinationPath = destinationPath, fun = "raster::raster")

  rs <- stack(lr)
  names(rs) <- CClayerNames2

  ## loop/apply doesn't work here???
  ml <- mapAdd(rs[[1]], map = ml, layerName = CClayerNames2[1], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(rs[[2]], map = ml, layerName = CClayerNames2[2], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(rs[[3]], map = ml, layerName = CClayerNames2[3], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(rs[[4]], map = ml, layerName = CClayerNames2[4], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(rs[[5]], map = ml, layerName = CClayerNames2[5], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")
  ml <- mapAdd(rs[[6]], map = ml, layerName = CClayerNames2[6], CC = TRUE,
               destinationPath = destinationPath,
               filename2 = tempfile(fileext = ".tif"), leaflet = FALSE, method = "ngb")

  ccs <- ml@metadata[CC == TRUE & !(layerName == "LandType"), ]
  CCs <- maps(ml, layerName = ccs$layerName)
  CCstack <- raster::stack(CCs)
  CCstackNames <- names(CCstack)

  if (!all(raster::minValue(CCstack) >= 0)) stop("problem with minValue of CCstack (< 0)")
  if (!all(raster::maxValue(CCstack) <= 10)) stop("problem with maxValue of CCstack (> 10)")

  CCstack <- CCstack * 10 # convert back to percent
  NA_ids <- which(is.na(ml$LandType[]) |  # outside of studyArea polygon
                    ml$LandType[] == 1)   # 1 is cities -- NA it here -- will be filled in with another veg layer if available (e.g. Pickell)
  message("  Setting NA, 1 in LandType to NA in speciesLayers in MB FRI data")
  aa <- try(CCstack[NA_ids] <- NA, silent = TRUE)
  ## unclear why line above sometimes fails: 'Error in value[j, ] : incorrect number of dimensions'
  if (is(aa, "try-error")) {
    l <- unstack(CCstack)
    CCstack <- lapply(l, function(x) {x[NA_ids] <- NA; x})
  }
  names(CCstack) <- equivalentName(CCstackNames, sppEquiv, sppEquivCol)

  stack(CCstack)
}

#' @export
#' @importFrom map mapAdd maps
#' @importFrom raster addLayer dropLayer maxValue minValue stack unstack
#' @rdname prepSpeciesLayers
prepSpeciesLayers_ONFRI <- function(destinationPath, outputPath,
                                    url = NULL,
                                    studyArea, rasterToMatch,
                                    sppEquiv,
                                    sppEquivCol, ...) {

  ## TODO: this is sneaky and annoying (runName is part of outputPath)
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
  options(map.useParallel = FALSE) ## TODO: pass additional arg to function
  ml <- mapAdd(rasterToMatch, isRasterToMatch = TRUE, layerName = "rasterToMatch", filename2 = NULL)

  ml <- mapAdd(studyArea, map = ml, isStudyArea = TRUE, layerName = "studyArea",
               useSAcrs = TRUE, filename2 = NULL)

  ml <- mapAdd(map = ml, url = url, layerName = FRIlayerNames, CC = TRUE,
               destinationPath = destinationPath,
               targetFile = FRIlayerNamesFiles, filename2 = NULL,
               alsoExtract = NA, leaflet = FALSE, method = "ngb")

  ml <- mapAdd(map = ml, url = url2, layerName = "LCC_FRI", CC = TRUE,
               destinationPath = destinationPath,
               targetFile = FRIlccname, filename2 = NULL,
               overwrite = TRUE, ## TODO: workaraund bug in prepInputs redownloading each time
               alsoExtract = NA, leaflet = FALSE, method = "ngb")

  ccs <- ml@metadata[CC == TRUE & !(layerName == "LCC_FRI"), ]
  CCs <- maps(ml, layerName = ccs$layerName)
  CCstack <- raster::stack(CCs)
  CCstackNames <- names(CCstack)

  if (!all(raster::minValue(CCstack) >= 0)) stop("problem with minValue of CCstack (< 0)")
  if (!all(raster::maxValue(CCstack) <= 100)) stop("problem with maxValue of CCstack (> 100)")

  ## merge species layers (currently only Popu; TODO: pine?)
  idsPopu <- grep("Popu", CCstackNames)
  mergedPopu <- calc(stack(CCstack[[idsPopu]]), sum, na.rm = TRUE)

  CCstack <- dropLayer(CCstack, idsPopu)
  CCstack[["Popu_sp"]] <- mergedPopu ## NOTE: addLayer sporadically fails to add layer, w/o warning

  idThuj <- grep("Thuj_spp", names(CCstack))
  names(CCstack[[idThuj]]) <- "Thuj_sp"

  stack(CCstack) ## ensure it's still a stack
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
#' @importFrom reproducible asPath Cache
#' @importFrom raster NAvalue<- raster rasterOptions setValues stack
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
  spRas <- raster(PickellRaster) |>
    setValues(NA_integer_)

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
        NAval <- 255L
        spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                               filename = asPath(file.path(destinationPath,
                                                           paste0("Pickell_", sp, ".tif"))),
                               overwrite = TRUE, datatype = "INT1U", NAflag = NAval)
        ## NAvals need to be converted back to NAs
        NAvalue(spRasts[[sp]]) <- NAval
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
        NAvalue(spRasts[[sp]]) <- NAval
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
        NAvalue(spRasts[[sp]]) <- NAval
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
        NAvalue(spRasts[[sp]]) <- NAval
      }
    }
  }

  ## species in Pickell's data
  raster::stack(spRasts)
}


#' Convert NA's in speciesLayers to zeros
#'
#' Pixels that have NAs but are inside `rasterToMatch`
#' may need to be converted to 0s as they can could still potentially
#' be forested
#'
#' @template speciesLayers
#' @template rasterToMatch
#'
#' @return the `speciesLayers` with 0s in pixels that had NAs
#'
#' @export
#' @importFrom raster stack
NAcover2zero <- function(speciesLayers, rasterToMatch) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    ## since terra is dependency of raster, it should already be installed, but just in case...
    stop("Suggested package 'terra' not installed.\n",
         "Install it using `install.packages('terra')`.")
  }

  tempRas <- rasterToMatch
  tempRas[!is.na(tempRas[])] <- 0
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
