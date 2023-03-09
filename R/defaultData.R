utils::globalVariables(c(
  "TSLF", "latitude", "longitude", "Degree_Of_", "Permafrost", "thermokrst"
))

## ------------------------------------
## SOURCING FUNCTIONS FOR DEFAULT DATA
## ------------------------------------

## LandR may not be the best pacakge for these funcitons, but they will live
## here for now.

#' Default environmental variables
#'
#' Downloads a set of default environmental variable rasters to model
#' maxB, SEP and thermokarst (amongst others).
#'
#' @param vars character. Variables to download.
#' @param studyArea a spatial polygon  (`sf`, `sp` or `SpatVector`) used to crop
#'   and mask data layers. Passed to `reproducible::Cache` and `reproducible::prepInputs`
#' @param rasterToMatch a raster layer (`RasterLayer` or `SpatRaster`). Only used to `rasterize`
#'   the default permafrost layer and fire perimeters.
#' @param userTags passed to `reproducible::Cache` and `reproducible::prepInputs`
#' @param destinationPath passed to `reproducible::Cache` and `reproducible::prepInputs`
#'
#' @details By default this function downloads four bioclimatic variables (mean annual temperature
#'   MAT, precipitation of warmest quarter, PPT_sm, precipitation of coldest quarter, PPT_wt, annual
#'   climate moisture index, CMI) obtained from AdaptWest for the reference period 1980-2010
#'   (https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6/normals/Normal_1981_2010_bioclim.zip),
#'   elevation (from https://adaptwest.databasin.org/pages/adaptwest-climatena/) - at 1 Km resolution -
#'   a binary map of wetlands (from Wulder et al. 2018), permafrost peatland complex cover and amount
#'   of thermokarst (from Gibson et al. 2021) -- which are combined to produce a compound variable called
#'   'thermXperm' -- and a map of time since fire (considering fire perimeters from the Canadian Wildland
#'   Fire Information System.
#'
#' @return named list of `SpatRasters`, with names following `vars`.
#'
#' @importFrom reproducible Cache prepInputs postProcess normPath
#' @importFrom sf st_as_sf st_write st_crs st_crs<- st_transform st_intersection
#' @importFrom terra aggregate mask rasterize rast levels values values<-
#' @export
defaultEnvirData <- function(vars = c("MAT", "PPT_wt", "PPT_sm", "CMI", "elevation",
                                      "timeSinceFire", "wetlands", "permafrost", "thermokarst",
                                      "thermXperm"),
                             studyArea, userTags, destinationPath, rasterToMatch) {
  vars <- match.arg(vars, several.ok = TRUE)

  outs <- list()

  if ("timeSinceFire" %in% vars) {
    sa <- if (is(studyArea, "sf")) {
      aggregate(studyArea, list(rep(1, nrow(studyArea))),
                FUN = function(x) x)
    } else {
      aggregate(studyArea)
    }
    firePerimeters <- Cache(
      prepInputsFireYear(fireField = "YEAR",
                         rasterToMatch = rasterToMatch,
                         earliestYear = 1900, ## go as far back as possible
                         url = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip",
                         destinationPath = destinationPath,
                         studyArea = sa,
                         maskWRTM = TRUE,
                         overwrite = TRUE,
                         userTags = c(userTags, "firePerimeters"),
                         omitArgs = c("userTags")),
      userTags = c(userTags, "firePerimeters")
    )

    ## from Baltzer et al 2021 synthesis paper
    ## we're only interested in TSLF (time since last fire) and the coordinates
    if (!inherits(studyArea, c("sf", "sfc"))) {
      SAsf <- st_as_sf(studyArea)
    } else {
      SAsf <- studyArea
    }

    fireData <- Cache(prepInputs,
                      targetFile = "TaigaPlains.xlsx",
                      url = "https://docs.google.com/spreadsheets/d/1D1nrKkygfjhO53K5f44DcF9m6mhTgTx5/edit?usp=sharing&ouid=101605422842478560921&rtpof=true&sd=true",
                      fun = "readxl::read_xlsx",
                      destinationPath = destinationPath,
                      userTags = c(userTags, "fireDataNWT"),
                      omitArgs = c("userTags")) |>
      as.data.table()
    fireData <- st_as_sf(fireData[, .(TSLF, latitude, longitude)], coords = c("longitude", "latitude"))
    st_write(obj = fireData, dsn = file.path(destinationPath, "TSLF.shp"), delete_layer = TRUE)
    st_crs(fireData) <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"  ## latlong
    fireData <- st_transform(fireData, crs = st_crs(SAsf))

    ## are there any points inside the study Area?
    test <- st_intersection(fireData, SAsf)
    if (!nrow(test)) {
      message("No fire data intersects with study area. ",
              "\nAll available data for the Taiga Plains will be used")
    } else {
      if (nrow(test) < 10) {
        message("There are less than ten records of time since ladt fire in the study area.",
                "\nAll available data for the Taiga Plains will be used")
      } else {
        fireData <- test
      }
    }

    ## exclude fires post-2015
    vals <- values(firePerimeters)
    vals[vals > 2015 ,] <- NA
    values(firePerimeters) <- vals

    ## time since fire raster:
    timeSinceFire <- 2015 - firePerimeters

    ## need to assume a year for NA values (otherwise biomod2 can't predict)
    ## using the mid point between earliest fire and 1900
    vals <- values(timeSinceFire)
    vals[is.na(vals),] <- mean(fireData$TSLF)
    values(timeSinceFire) <- vals
    timeSinceFire <- mask(timeSinceFire, rasterToMatch)
    outs[["timeSinceFire"]] <- timeSinceFire
  }

  if ("elevation" %in% vars) {
    outs[["elevation"]] <- Cache(prepInputs,
                               targetFile = "elevation.tif",
                               url = "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6/elevation.tif",
                               destinationPath = destinationPath,
                               alsoExtract = NA,
                               overwrite = TRUE,
                               # rasterToMatch = rasterToMatch,
                               studyArea = studyArea,   ## keep coarse res.
                               useSAcrs = FALSE,
                               userTags = c(userTags, "elevation"),
                               omitArgs = c("userTags"))
  }
  ## Mean annual temperature (mean of all of the monthly mean temperatures)
  bioClimURL <- paste0("https://s3-us-west-2.amazonaws.com/www.cacpd.org/",
                       "CMIP6/normals/Normal_1981_2010_bioclim.zip")
  if ("MAT" %in% vars) {
    outs[["MAT"]] <- Cache(prepInputs,
                         targetFile = "Normal_1981_2010_MAT.tif",
                         archive = "data/Normal_1981_2010_bioclim.zip",
                         url = bioClimURL,
                         destinationPath = destinationPath,
                         alsoExtract = NA,
                         overwrite = TRUE,
                         # rasterToMatch = rasterToMatch,
                         studyArea = studyArea,   ## keep coarse res.
                         useSAcrs = FALSE,
                         userTags = c(userTags, "MAT"),
                         omitArgs = c("userTags"))
  }

  ## Precipitation of warmest quarter (total precipitation over warmest quarter of the year)
  ## assuming warmest quarter to be June-August
  if ("PPT_sm" %in% vars) {
    outs[["PPT_sm"]] <- Cache(prepInputs,
                            targetFile = "Normal_1981_2010_PPT_sm.tif",
                            archive = "data/Normal_1981_2010_bioclim.zip",
                            url = bioClimURL,
                            destinationPath = destinationPath,
                            alsoExtract = NA,
                            overwrite = TRUE,
                            # rasterToMatch = rasterToMatch,
                            studyArea = studyArea,   ## keep coarse res.
                            useSAcrs = FALSE,
                            userTags = c(userTags, "PPT_sm"),
                            omitArgs = c("userTags"))
  }

  ## Precipitation of coldest quarter (total precipitation over the coldest quarter of the year)
  ## assuming coldest quarter to be December-February
  if ("PPT_wt" %in% vars) {
    outs[["PPT_wt"]] <- Cache(prepInputs,
                            targetFile = "Normal_1981_2010_PPT_wt.tif",
                            archive = "data/Normal_1981_2010_bioclim.zip",
                            url = bioClimURL,
                            destinationPath = destinationPath,
                            alsoExtract = NA,
                            overwrite = TRUE,
                            # rasterToMatch = rasterToMatch,
                            studyArea = studyArea,   ## keep coarse res.
                            useSAcrs = FALSE,
                            userTags = c(userTags, "PPT_wt"),
                            omitArgs = c("userTags"))
  }

  ## Annual climate moisture index (CMI = P PET, where P is annual precipitation and PET is annual potential evapotranspiration)
  ## this is Hogg's CMI
  if ("CMI" %in% vars) {
    outs[["CMI"]] <- Cache(prepInputs,
                         targetFile = "Normal_1981_2010_CMI.tif",
                         archive = "data/Normal_1981_2010_bioclim.zip",
                         url = bioClimURL,
                         destinationPath = destinationPath,
                         alsoExtract = NA,
                         overwrite = TRUE,
                         # rasterToMatch = rasterToMatch,
                         studyArea = studyArea,   ## keep coarse res.
                         useSAcrs = FALSE,
                         userTags = c(userTags, "CMI"),
                         omitArgs = c("userTags"))
  }

  ## --------------------------------------------------------
  ## Wetland data
  ## binary raster of wetlands 2000-2016 (wetland presence = wetlands were present >=80% of the years)
  if ("wetlands" %in% vars) {
    wetlands <- Cache(prepInputs,
                      targetFile = "CA_wetlands_post2000.tif",
                      archive = "CA_wetlands_post2000.zip",
                      url = "https://opendata.nfis.org/downloads/forest_change/CA_wetlands_post2000.zip",
                      destinationPath = destinationPath,
                      # rasterToMatch = rasterToMatch,
                      studyArea = studyArea,   ## keep coarse res.
                      useSAcrs = FALSE,
                      method = "near",
                      overwrite = TRUE,
                      userTags = c(userTags, "wetlands"),
                      omitArgs = c("userTags"))
    ## replace NAs inside study area with 0s. otherwise we loose data for SDMs
    sa <- postProcessTerra(studyArea, projectTo = wetlands)
    wetlands[is.na(wetlands[])] <- 0
    outs[["wetlands"]] <- mask(wetlands, sa)
    rm(sa)
  }

  ## --------------------------------------------------------
  ## Permafrost and thermokarst data
  ## not using raster::shapefile results in a lot of alignment issues!
  ## the URL supplied cannot be used to retrieve the data:
  if (any(c("permafrost", "thermokarst", "thermXperm") %in% vars)) {
    cmd <- paste("curl --header 'Host: app.nwtgeoscience.ca'",
                 "--user-agent 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:99.0) Gecko/20100101 Firefox/99.0'",
                 "--header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8'",
                 "--header 'Accept-Language: en-CA,en-US;q=0.7,en;q=0.3'",
                 "--header 'Content-Type: application/x-www-form-urlencoded'",
                 "--header 'Origin: https://app.nwtgeoscience.ca'",
                 "--referer 'https://app.nwtgeoscience.ca/Journal.aspx?refnum=2020-010'",
                 "--cookie 'ASP.NET_SessionId=vzqpm3xbcsv0oorjkxjd2cy3'",
                 "--header 'Upgrade-Insecure-Requests: 1'",
                 "--header 'Sec-Fetch-Dest: document'",
                 "--header 'Sec-Fetch-Mode: navigate'",
                 "--header 'Sec-Fetch-Site: same-origin'",
                 "--header 'Sec-Fetch-User: ?1'",
                 "--request POST",
                 "--data-urlencode",
                 paste0("'__VIEWSTATE=VgcB84SNFQjOzui1LrkB/baCjwp+82NVTO63xsRzN3tIbPwB8",
                        "PoWe79ARC2a847RWR4FqYMEn6c16kqkMF5E6R8b6CCNoK3+E+e+uLDuTU5yqjU",
                        "0pPxhn3Djn/jJ6+VEl/6NkdFIieo/HDBwZYq9cp2Gj6R7OlR4D37gxjO7m9JbvI",
                        "WDJHA6X0qa8Z/+dR/lm8c9edrHUl76DXVAB+ja6wVt5ry8clL/U5bgWzHZG4hcC",
                        "D2L5Wh+ScFUQD++qTrCxloNafdxCMwfBDB1t4/pBz0NTKUz+zRiS7hSY7kl/1dM",
                        "88UISQgkEZTkCtdTNCvTnBq3okXggHDmTU9zZhAP9hfnsJbymnhqioh+Fh50kEd",
                        "s9IN1v2SVKVqoQvqy8bAiaBQw2ej7YwmJDeYiGcTS93Ig6Q1QF4xZ2dgKbpnsT2",
                        "KhB3RCtqQiCjJOPvZ48aH4w1RGISj+nV/vDDGkkq0/HvKvHEUa9NxbzsCkXNvFz",
                        "QteyQcPC48cG8WDgVVgY2gwmnGNM14eohVXtb58IlTyQ3uQgpw7lsc2nwyMKqp0",
                        "RoIqD7d2F7H62MaPXQhU3iTntBTOWCzYfUaF1tFHoMnOReM2suFS7z0J0X2c3E9",
                        "8hDV5vjrgixkqNsvQdifeW0rf4JQ68a0I/yePF9G7U2Rt+r2aneXd5eE8dLB+vY",
                        "K0Lbjh3JGAlxteN0r+73W3KNTZEAh/yu7+zbyebgSKnwcV6B4fld1utvgusVhHM",
                        "IR9KlIluB4NVB/s01CbuLtQ/mAvs5zJn7yQ9pkyB8V0grZdC0Nbj/aoZyDd6jrB",
                        "GbbgYMVnXz13aDodQgQKWbAoqMDqA0cgw6bXffY6UTh7KxXtWVYgnxpJ47CyzK5",
                        "cyiPO/RsqcBYCkz3ETVLutx1Rr04x/UymqR6xyoRehjGM+l9N4B3GkmDy0B0D09",
                        "4t7bcn7jxQLh9DJC7VlzKlWbrgE4HJTCV+dCBTeSAUC+KIGKjxWXoF04GzHCcEbN",
                        "i/57wPJ5s='"),
                 "--data-urlencode '__VIEWSTATEGENERATOR=4A21451A'",
                 "--data-urlencode",
                 paste0("'__EVENTVALIDATION=kSI1vOPq7X5qViNicd5Yp7wqlvVE0N/Qqn4+2sgSPE/gk",
                        "pUsNlHYEm37WmFa44YryKY7rd+/1sxYaWwb6gVVymojca9ybYAYibQbUVGiOvHNm",
                        "aUkJn7LEJW8NIGn0vvp'"),
                 "--data-urlencode 'btnDownload=Download' 'https://app.nwtgeoscience.ca/Journal.aspx?refnum=2020-010'",
                 "--output 'DOI_2020-010.zip'")

    fileBash <- normPath(file.path(destinationPath, "curlCommand.sh"))
    cat(paste("cd", destinationPath), file = fileBash, sep = "\n")
    cat(cmd, file = fileBash, append = TRUE, sep = "\n")

    out <- tryCatch(shell(fileBash), error = function(e) e)
    if (is(out, "error")) {
      warning("Could not download DOI_2020-010.zip (permafrost data).",
              " Please download manually from 'https://doi.org/10.46887/2020-010' and put zip file in ",
              destinationPath)
    }
    file.remove(fileBash)

    permafrost <- Cache(prepInputs,
                        targetFile = "Final_Organics_TaigaPlains.shp",
                        alsoExtract = "similar",
                        archive = "DOI_2020-010.zip",
                        destinationPath = destinationPath,
                        studyArea = studyArea,
                        useSAcrs = FALSE,
                        overwrite = TRUE,
                        fun = "terra::vect",
                        # useCache = "overwrite",
                        userTags = c(userTags, "permafrost"),
                        omitArgs = c("userTags"))

    ## "rasterize" amount of permafrost
    ## make a points object to extract values in native CRS
    pixIDs <- which(!is.na(as.vector(rasterToMatch[])))
    pixIDsCoords <- as.data.frame(xyFromCell(rasterToMatch, cell = pixIDs))
    pixIDsPoints <- vect(pixIDsCoords, geom = c("x", "y"), crs = crs(rasterToMatch))
    pixIDs <- data.table(ID = 1:length(pixIDs), pixelIndex = pixIDs)   ## these will be the points IDs

    permafrostData <- genericExtract(x = permafrost, y = pixIDsPoints, field = "Permafrost")
    ## points in between polys could touch more than one, take min value
    ## for a conservative estimate of permafrost
    if (any(duplicated(permafrostData$ID))) {
      suppressWarnings({
        permafrostData <- permafrostData[, list(Permafrost = min(Permafrost, na.rm = TRUE)),
                                         by = ID]
      })
      permafrostData[Permafrost == Inf, Permafrost := NA]
    }
    permafrostData <- pixIDs[permafrostData, on = .(ID)]
    permafrostRas <- rasterToMatch
    permafrostRas[permafrostData$pixelIndex] <- permafrostData$Permafrost

    ## "rasterize" amount of thermokarst
    thermokarstData <- genericExtract(x = permafrost, y = pixIDsPoints, field = "Degree_Of_")
    thermokarstData[, thermokrst := factor(Degree_Of_, levels = c("None", "Low", "Medium", "High"),
                                           labels = c(0, 1, 2, 3))]
    thermokarstData[, thermokrst := as.integer(as.character(thermokrst))]
    ## points in between polys could touch more than one, take min value
    ## for a pessimistic estimate of thermokarst
    if (any(duplicated(thermokarstData$ID))) {
      suppressWarnings({
        thermokarstData <- thermokarstData[, list(thermokrst = max(thermokrst, na.rm = TRUE)),
                                         by = ID]
        thermokarstData[thermokrst == Inf, thermokrst := NA]
      })
    }
    thermokarstData <- pixIDs[thermokarstData, on = .(ID)]
    thermokarstRas <- rasterToMatch
    thermokarstRas[thermokarstData$pixelIndex] <- thermokarstData$thermokrst

    ## /!\ there are "NAs"/"NULLs" in thermokarst areas (Degree_Of) that have permafrost.
    ## These are probably true zeros
    ## note: this has to be done after rasterizing -- doing it before lead to mismatches
    trueZeros <- which(is.na(values(thermokarstRas)) & !is.na(values(permafrostRas)))
    values(thermokarstRas)[trueZeros] <- 0

    levelsthermokarstRas <- data.frame(VALUE = c(0, 1, 2, 3), LEGEND = c("None", "Low", "Medium", "High"))
    levels(thermokarstRas) <- levelsthermokarstRas

    if ("permafrost" %in% vars) {
      outs[["permafrost"]] <- permafrostRas
    }

    if ("thermokarst" %in% vars) {
      outs[["thermokarst"]] <- thermokarstRas
    }

    ## combine thermokarst and permafrost layers
    ## high thermokarst in low permafrost areas != than in high permafrost areas (especially for veg)
    ## perhaps multiplying is misleading too. is 25% PCC * 2 thermokarst (medium) the same as 50% permafrost * 1 (low) thermokarst?
    ## actually we're lucky with Gibson's data because the permafrost classes are not multiples of each other - however this approach is not ideal when they are...
    if ("thermXperm" %in% vars) {
      outs[["thermXperm"]] <- thermokarstRas*permafrostRas
    }
  }
  outs
}


#' Default future projections (forecasts) of climate variables
#'
#' Downloads a set of default climate variables' forecasts as rasters -- used
#' to model maxB, SEP and thermokarst (amongst others).
#'
#' @param vars character. Variables to download.
#' @param climateGCM character. The Global Circulation Model, GCM, (or ensemble) to use.
#' @param climateSSP character. The Socio-economic Shared Pathway, SSP (i.e. emisions scenario) to use.
#' @param periods character. The year periods to use.
#' @param studyArea a spatial polygon  (`sf`, `sp` or `SpatVector`) used to crop
#'   and mask data layers. Passed to `reproducible::Cache` and `reproducible::prepInputs`
#' @param rasterToMatch a raster layer (`RasterLayer` or `SpatRaster`). Only used to `rasterize`
#'   the default permafrost layer and fire perimeters.
#' @param userTags passed to `reproducible::Cache` and `reproducible::prepInputs`
#' @param destinationPath passed to `reproducible::Cache` and `reproducible::prepInputs`
#'
#' @details By default this function downloads four bioclimatic variables (mean annual temperature
#'   MAT, precipitation of warmest quarter, PPT_sm, precipitation of coldest quarter, PPT_wt, annual
#'   climate moisture index, CMI) obtained from AdaptWest for 20 year periods between 2001 and 2100
#'   (see https://adaptwest.databasin.org/pages/adaptwest-climatena/) - at 1 Km resolution].
#'
#' @return named list of `SpatRasters`, with names following `vars`.
#'
#' @importFrom reproducible Cache prepInputs postProcess
#' @export
defaultClimateDataProj <- function(vars = c("MAT", "PPT_wt", "PPT_sm", "CMI"),
                                   climateGCM = "ensemble_8GCMs", climateSSP = "585",
                                   periods = c("2001_2020", "2021_2040", "2041_2060", "2061_2080", "2081_2100"),
                                   studyArea, userTags, destinationPath, rasterToMatch) {
  ## checks
  vars <- match.arg(vars, several.ok = TRUE)

  if (length(climateGCM) > 1) {
    stop("Provide a single climateGCM")
  }
  climateGCM <- match.arg(climateGCM)

  if (length(climateSSP) > 1) {
    stop("Provide a single climateGCM")
  }
  climateSSP <- match.arg(climateSSP)

  bioClimURL <- "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/"

  if (grepl("ensemble", climateGCM)) {
    bioClimURL <- paste0(bioClimURL, "ensembles/", climateGCM, "_ssp", climateSSP)

  } else {
    bioClimURL <- paste0(bioClimURL, "selectedaogcms/", climateGCM, "_ssp", climateSSP)
  }

  bioClimURLs <- paste0(paste(bioClimURL, periods, sep = "_"),"_bioclim.zip")

  ## Mean annual temperature (mean of all of the monthly mean temperatures)
  if ("MAT" %in% vars) {
    MATfiles <- sub("bioclim.zip", "MAT.tif", basename(bioClimURLs))
    userTags <- lapply(MATfiles, function(x) c(userTags, x))
    MATProj <- Map(f = prepInputs,
                   archive = bioClimURLs,
                   targetFile = MATfiles,
                   userTags = userTags,
                   url = bioClimURLs,
                   MoreArgs =  list(
                     destinationPath = destinationPath,
                     alsoExtract = NA,
                     overwrite = TRUE,
                     # rasterToMatch = rasterToMatch,
                     studyArea = studyArea,   ## keep coarse res.
                     useSAcrs = FALSE,
                     omitArgs = c("userTags"))) |>
      Cache(userTags = c(userTags, "MATProj"),
            omitArgs = c("userTags"))
    names(MATProj) <- paste0("MAT_", sub("_.*", "", periods))
  }

  ## Precipitation of warmest quarter (total precipitation over warmest quarter of the year)
  ## assuming warmest quarter to be June-August
  if ("PPT_sm" %in% vars) {
    PPT_smfiles <- sub("bioclim.zip", "PPT_sm.tif", basename(bioClimURLs))
    userTags <- lapply(PPT_smfiles, function(x) c(userTags, x))
    PPT_smProj <- Map(f = prepInputs,
                      archive = bioClimURLs,
                      targetFile = PPT_smfiles,
                      userTags = userTags,
                      url = bioClimURLs,
                      MoreArgs =  list(
                        destinationPath = destinationPath,
                        alsoExtract = NA,
                        overwrite = TRUE,
                        # rasterToMatch = rasterToMatch,
                        studyArea = studyArea,   ## keep coarse res.
                        useSAcrs = FALSE,
                        omitArgs = c("userTags"))) |>
      Cache(userTags = c(userTags, "PPT_smProj"),
            omitArgs = c("userTags"))
    names(PPT_smProj) <- paste0("PPT_sm_", sub("_.*", "", periods))
  }

  ## Precipitation of coldest quarter (total precipitation over the coldest quarter of the year)
  ## assuming coldest quarter to be December-February
  if ("PPT_wt" %in% vars) {
    PPT_wtfiles <- sub("bioclim.zip", "PPT_wt.tif", basename(bioClimURLs))
    userTags <- lapply(PPT_wtfiles, function(x) c(userTags, x))
    PPT_wtProj <- Map(f = prepInputs,
                      archive = bioClimURLs,
                      targetFile = PPT_wtfiles,
                      userTags = userTags,
                      url = bioClimURLs,
                      MoreArgs =  list(
                        destinationPath = destinationPath,
                        alsoExtract = NA,
                        overwrite = TRUE,
                        # rasterToMatch = rasterToMatch,
                        studyArea = studyArea,   ## keep coarse res.
                        useSAcrs = FALSE,
                        omitArgs = c("userTags"))) |>
      Cache(userTags = c(userTags, "PPT_wtProj"),
            omitArgs = c("userTags"))
    names(PPT_wtProj) <- paste0("PPT_wt_", sub("_.*", "", periods))
  }

  ## Annual climate moisture index (CMI = P PET, where P is annual precipitation and PET is annual potential evapotranspiration)
  ## this is Hogg's CMI
  if ("CMI" %in% vars) {
    CMIfiles <- sub("bioclim.zip", "CMI.tif", basename(bioClimURLs))
    userTags <- lapply(CMIfiles, function(x) c(userTags, x))
    CMIProj <- Map(f = prepInputs,
                   archive = bioClimURLs,
                   targetFile = CMIfiles,
                   userTags = userTags,
                   url = bioClimURLs,
                   MoreArgs =  list(
                     destinationPath = destinationPath,
                     alsoExtract = NA,
                     overwrite = TRUE,
                     # rasterToMatch = rasterToMatch,
                     studyArea = studyArea,   ## keep coarse res.
                     useSAcrs = FALSE,
                     omitArgs = c("userTags"))) |>
      Cache(userTags = c(userTags, "CMIProj"),
            omitArgs = c("userTags"))
    names(CMIProj) <- paste0("CMI_", sub("_.*", "", periods))
  }

  outs <- c(as.list(MATProj), as.list(PPT_wtProj), as.list(PPT_smProj), as.list(CMIProj))
  outs
}
