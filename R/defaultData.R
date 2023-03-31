utils::globalVariables(c(
  "TSLF", "latitude", "longitude", "Degree_Of_", "Permafrost", "thermokarst",
  "dists"
))

## -----------------------------------------------
## FUNCTIONS TO SOURCE AND PREPARE DEFAULT DATA
## -----------------------------------------------

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
  if (any(c("permafrost", "thermokarst", "thermXperm") %in% vars)) {
    permafrost <- getPermafrostDataGibson(destinationPath, studyArea, userTags)
    outs[["permafrostPoly"]] <- permafrost
  }
  outs
}


#' Download permafrost/thermokarst layer
#'
#' Downloads and crops/masks permafrost peatland complex cover and amount of thermokarst
#' spatial data from Gibson et al. (2021)
#'
#' @param studyArea a spatial polygon  (`sf`, `sp` or `SpatVector`) used to crop
#'   and mask data layers. Passed to `reproducible::Cache` and `reproducible::prepInputs`
#' @param destinationPath passed to `reproducible::Cache` and `reproducible::prepInputs`
#' @param userTags passed to `reproducible::Cache` and `reproducible::prepInputs`
#'
#' @return Carolyn Gisbon's output polygon layer of % PPC and thermokarst state
#'
#' @importFrom reproducible prepInputs normPath
#' @export
getPermafrostDataGibson <- function(destinationPath, studyArea, cacheTags) {

  # outs <- list()
  ## the URL supplied cannot be used to retrieve the data:
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
                      userTags = c(cacheTags, "permafrost"),
                      omitArgs = c("userTags"))

  # ## "rasterize" amount of permafrost
  # ## make a points object to extract values in native CRS
  # pixIDs <- which(!is.na(as.vector(rasterToMatch[])))
  # pixIDsCoords <- as.data.frame(xyFromCell(rasterToMatch, cell = pixIDs))
  # pixIDsPoints <- vect(pixIDsCoords, geom = c("x", "y"), crs = crs(rasterToMatch))
  # pixIDs <- data.table(ID = 1:length(pixIDs), pixelIndex = pixIDs)   ## these will be the points IDs
  #
  # permafrostData <- genericExtract(x = permafrost, y = pixIDsPoints, field = "Permafrost")
  # permafrostData[is.na(Permafrost), Permafrost := 0] ## In NA's in RTM mean 0% permafrost (as rasterize(..., background = 0))
  #
  # ## points in between polys could touch more than one, take min value
  # ## for a conservative estimate of permafrost
  # if (any(duplicated(permafrostData$ID))) {
  #   suppressWarnings({
  #     permafrostData <- permafrostData[, list(Permafrost = min(Permafrost, na.rm = TRUE)),
  #                                      by = ID]
  #   })
  # }
  # permafrostData <- pixIDs[permafrostData, on = .(ID)]
  # permafrostRas <- rasterToMatch
  # permafrostRas[permafrostData$pixelIndex] <- permafrostData$Permafrost
  # permafrostRas <- mask(permafrostRas, rasterToMatch)
  #
  # ## "rasterize" amount of thermokarst
  # thermokarstData <- genericExtract(x = permafrost, y = pixIDsPoints, field = "Degree_Of_")
  # thermokarstData[, thermokarst := factor(Degree_Of_, levels = c("None", "Low", "Medium", "High"),
  #                                        labels = c(0, 1, 2, 3))]
  # thermokarstData[, thermokarst := as.integer(as.character(thermokarst))]
  #
  # ## /!\ there are "NAs"/"NULLs" in thermokarst areas (Degree_Of) that have permafrost.
  # ## These are probably true zeros
  # ## note: this has to be done after rasterizing -- doing it before lead to mismatches
  # thermokarstData <- thermokarstData[permafrostData, on = .(ID)]
  # thermokarstData[is.na(thermokarst) & !is.na(Permafrost), thermokarst := 0L]
  #
  # ## points in between polys could touch more than one, take min value
  # ## for a pessimistic estimate of thermokarst
  # if (any(duplicated(thermokarstData$ID))) {
  #   suppressWarnings({
  #     thermokarstData <- thermokarstData[, list(thermokarst = max(thermokarst, na.rm = TRUE)),
  #                                      by = ID]
  #   })
  # }
  # thermokarstData <- pixIDs[thermokarstData, on = .(ID)]
  # thermokarstRas <- rasterToMatch
  # thermokarstRas[thermokarstData$pixelIndex] <- thermokarstData$thermokarst
  # thermokarstRas <- mask(thermokarstRas, rasterToMatch)
  #
  # levelsthermokarstRas <- data.frame(VALUE = c(0, 1, 2, 3), LEGEND = c("None", "Low", "Medium", "High"))
  # levels(thermokarstRas) <- levelsthermokarstRas
  #
  # if ("permafrost" %in% vars) {
  #   outs[["permafrost"]] <- permafrostRas
  # }
  #
  # if ("thermokarst" %in% vars) {
  #   outs[["thermokarst"]] <- thermokarstRas
  # }
  #
  # ## combine thermokarst and permafrost layers
  # ## high thermokarst in low permafrost areas != than in high permafrost areas (especially for veg)
  # ## perhaps multiplying is misleading too. is 25% PCC * 2 thermokarst (medium) the same as 50% permafrost * 1 (low) thermokarst?
  # ## actually we're lucky with Gibson's data because the permafrost classes are not multiples of each other - however this approach is not ideal when they are...
  # if ("thermXperm" %in% vars) {
  #   outs[["thermXperm"]] <- thermokarstRas*permafrostRas
  # }

  # outs[["permafrostPoly"]] <- permafrost
  permafrost
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



#' Make a mask of areas that suitable/unsuitable for permafrost
#'
#' Uses a land-cover map and a wetlands maps to create a mask
#'  of areas suitable for permafrost presence. By default, it
#'  assumes suitable and unsuitable land-cover classification
#'  follows [Hermosilla et al. (2022)](https://www.sciencedirect.com/science/article/pii/S0034425721005009?via%3Dihub#bb0465).
#'  See details below.
#'
#' @param rstLCC a land-cover raster. Should already have been cropped to
#'   `studyArea`.
#' @param wetlands a wetlands raster (with `1L` and `NA`s only). Should already
#'   have been cropped to `studyArea`.
#' @param suitableCls land-cover classes (in `rstLCC`) suitable
#'  for permafrost presence
#' @param unsuitableCls land-cover classes (in `rstLCC`) unsuitable
#'  for permafrost presence
#' @param studyArea a polygon of the study area to mask the
#'  output raster to.
#' @param cacheTags passed to `reproducible::Cache`
#'
#' @details Suitable areas are classified as `1L` and unsuitable
#'  areas as `0L` in the output raster.
#'  Unsuitable land-cover classes are:
#'    * non-forested ecozone, which was not mapped for LC (0)
#'    * water (20)
#'    * snow/ice (31)
#'    * rock/rubble (32)
#'    * exposed/barren land (33)
#'    * wetland (80)
#'    * wetland-treed (81)
#'  Wetland areas in `wetlands` (coded as 1L, i.e. wetland presence)
#'  are also coded as unsuitable.
#'
#' @return a raster layer with `1`s and `0`s for pixels with suitable
#'  and unsuitable land cover for permafrost, respectively.
#'
#' @importFrom terra classify subst mask
#' @importFrom raster raster
#' @importFrom reproducible Cache postProcess
#'
#' @export
makeSuitForPerm <- function(rstLCC, wetlands, suitableCls = c(40, 50, 100, 210, 220, 230),
                            unsuitableCls = c(0, 20, 31, 32, 33, 80, 81),
                            studyArea, cacheTags) {
  unsuit <- cbind(unsuitableCls, NA_integer_)
  suit <- cbind(suitableCls, 1L)
  m <- rbind(unsuit, suit)

  rstLCC <- Cache(classify,
                  x = rstLCC,
                  rcl = m,
                  userTags = c(cacheTags, "rstLCC", "suitablePermafrost"),
                  omitArgs = c("userTags"))
  rstLCC <- Cache(mask,
                  x = rstLCC,
                  mask = wetlands,
                  inverse = TRUE,
                  userTags = c(cacheTags, "rstLCC", "suitablePermafrost"),
                  omitArgs = c("userTags"))

  if (getOption("LandR.assertions", TRUE)) {
    test <- unique(rstLCC[!is.na(as.vector(values(wetlands)))])
    if (any(!is.na(test))) {
      stop("Something went wrong. All wetland areas should have 'NA' in suitable areas for permafrost")
    }
  }

  ## covert NAs to zeros then mask again
  rstLCC <- subst(rstLCC, NA, 0L)
  rstLCC <- postProcess(x = rstLCC,
                        studyArea = studyArea,
                        useSAcrs = FALSE,
                        userTags = c(cacheTags, "rstLCC", "suitablePermafrost"),
                        omitArgs = c("userTags"))

  return(rstLCC)
}



#' Create permafrost P/A layer
#'
#' Creates a permafrost presence/absence raster layer based
#'  on information about land-cover suitability for permafrost
#'  (a raster layer at high resolution) and % of permafrost
#'  present (a polygon layer that groups several cells of the
#'  land-cover layer, i.e. at a larger spatial scale). The function
#'  has been designed to be looped over polygons, and so only processes
#'  only polygon at a time.
#'
#' @param gridPoly a `SpatVector`, `PackedSpatVector` or `character`.
#'  The full polygon layer containing information about % permafrost
#'  (as a `SpatVector` or `PackedSpatVector`) or the file name to that
#'  layer.
#' @param ras a `SpatRaster`, `PackedSpatRaster` or `character`.
#'  The land-cover suitability for permafrost raster layer (as a
#'  `SpatRaster` or `PackedSpatRaster`) or the file name to that layer.
#'  Suitable cells are coded as `rasClass`.
#' @param saveOut logical. Should the processed rasters for the focal
#'  polygon be saved (as the file name returned) or directly returned?
#'  If parallelising, choose to save, as `SpatRasters` cannot be serialized.
#' @param saveDir character. The directiory to save output rasters.
#' @param id character or numeric. The polygon ID in `gridPoly[IDcol]`
#'  to process (the focal polygon).
#' @param IDcol character. Column name in `gridPoly` containing polygon IDs.
#' @param rasClass Class
#'
#' @return a file name or the permafrost presence/absence raster for the focal
#'  polygon.
#'
#' @importFrom terra rast vect unwrap writeRaster as.polygons as.points as.lines
#' @importFrom terra mask crop cellFromXY crds disagg distance expanse nearby fillHoles
#' @importFrom data.table as.data.table
#' @importFrom crayon cyan
#' @export
assignPermafrost <- function(gridPoly, ras, saveOut = TRUE, saveDir = NULL,
                             id = NULL, IDcol = "OBJECTID", rasClass = 1L) {
  message(cyan("Preparing layers..."))

  if (is(gridPoly, "PackedSpatVector")) {
    gridPoly <- unwrap(gridPoly)
  }
  if (is(ras, "PackedSpatRaster")) {
    ras <- unwrap(ras)
  }

  if (is(gridPoly, "character")) {
    gridPoly <- vect(gridPoly)
  }

  if (is(ras, "character")) {
    ras <- rast(ras)
  }

  if (is.null(id)) {
    id <- gridPoly[[IDcol]][1,]
  }
  idd <- which(gridPoly[[IDcol]] == id)
  landscape <- gridPoly[idd]
  sub_ras <- crop(ras, landscape, mask = TRUE, touches = FALSE)
  ## certain polygons may be in ras boundaries and so small
  ## that they don't overlap any cell centroids. So try again with touches = TRUE
  if (isFALSE(any(!is.na(as.vector(sub_ras[]))))) {
    sub_ras <- crop(ras, landscape, mask = TRUE)
  }

  ## the solution may not have worked and there may simply be no
  ## raster cells intersected by this polygon.
  skipThisPoly <- isFALSE(any(!is.na(as.vector(sub_ras[]))))

  ## make storage raster
  sub_rasOut <- sub_ras
  sub_rasOut[] <- NA_integer_

  if (skipThisPoly) {
    warning(paste0("None of 'ras' touch polygon id ", id,
                   ".\n  You may want to consider extending 'ras'."))
  } else {
    permpercent <- as.numeric(landscape[["Permafrost"]])

    ## and in pixels
    suitablePixNo <- sum(as.vector(sub_ras[]) == rasClass, na.rm = TRUE)
    permpercentPix <- (permpercent/100)*ncell(sub_ras)  ## don't round here, otherwise values <0.5 become 0.

    ## use max(..., 1) to guarantee that values lower than 1, get one pixel.
    if (permpercentPix > 0)
      permpercentPix <- max(permpercentPix, 1)

    permpercentPix <- round(permpercentPix)

    ## if there's less permafrost than suitable areas
    ## find the largest patch and a point that is distant from its edge
    ## then assign permafrost starting from this point
    ## until the percentage is reached

    message(cyan("Assigning permafrost..."))
    ## if there are more suitable areas than permafrost start
    ## "filling in" from focal pixels using distance to edges
    ## as the probability of a pixel being selected (more distant = higher prob)
    if (suitablePixNo > permpercentPix & permpercentPix > 0) {
      ## make polygons, so that we can identify separate patches in raster
      sub_poly <- as.polygons(sub_ras) |> disagg()   ## 0s in sub_ras are ignored
      names(sub_poly) <- "patchType"
      ## subset to patches of interest
      sub_poly <- sub_poly[sub_poly$patchType == rasClass, ]
      sub_poly$ID <- 1:nrow(sub_poly)

      ## now make sub_ras with poly IDs
      sub_ras2 <- rasterize(sub_poly, sub_ras, field = "ID")
      # terra::plot(sub_ras2, col = viridis::viridis(length(sub_poly)))

      ## compute distances to edges.
      ## first make an "inverse" raster with NAs
      sub_rasDist <- sub_ras2
      sub_rasDist[] <- 1L
      sub_rasDist <- mask(sub_rasDist, sub_ras2, inverse = TRUE)   ## use sub_ras2, because 0s were NAed
      sub_rasDist <- distance(sub_rasDist)

      ## make "probabilities" by scaling 0-1
      spreadProb <- sub_rasDist/max(sub_rasDist[], na.rm = TRUE)
      # terra::plot(spreadProb, col = viridis::inferno(100))

      ## we may need to try several times until we get the number of pixels
      ## at each attempt increase spreadProb
      sub_rasOut <- assignPresences(assignProb = spreadProb,
                                    landscape = sub_ras2,
                                    pixToConvert = permpercentPix,
                                    probWeight = 0.5, numStartsDenom = 10)
      # terra::plot(sub_rasOut, col = viridis::inferno(100))
    }

    ## if there's more permafrost than suitable areas
    ## assign all suitable areas as permafrost, then increase
    ## with a buffer (starting with largest patch)
    if (suitablePixNo <= permpercentPix & permpercentPix > 0) {
      pixToConvert <- permpercentPix

      if (suitablePixNo > 0) {
        ## convert all available cells
        cellIDs <- which(as.vector(sub_ras[]) == rasClass)
        sub_rasOut[cellIDs] <- 1L

        ## if there are not enough points within the patches,
        ## try to fill neighbouring pixels (leave holes alone s we want to be
        ## able to start from a swiss cheese pattern)
        pixToConvert2 <- pixToConvert - length(cellIDs)

        while (pixToConvert2 > 0) {
          sub_poly <- as.polygons(sub_rasOut) |> disagg()

          sub_poly_filled <- fillHoles(sub_poly)
          sub_poly_filled <- rasterize(sub_poly_filled, sub_ras)

          ## buffer around patch
          sub_poly_filledBuffer <- buffer(sub_poly_filled, width = unique(res(sub_poly_filled)),
                                          background = NA)
          cellIDs <- which(as.vector(sub_poly_filledBuffer[]) == 1)

          ## remove cellIDs that have been converted already
          vals <- sub_rasOut[cellIDs][]
          cellIDs <- cellIDs[is.na(vals)]

          ## check that we don't have too many cells. if we do, keep
          ## pixels around larger patches
          if (length(cellIDs) > pixToConvert2) {
            sub_rasOut2 <- sub_rasOut
            sub_rasOut2[cellIDs] <- 1L

            sub_poly2 <- as.polygons(sub_rasOut2) |> disagg()
            sub_poly2$area <- expanse(sub_poly2)
            sub_rasOutArea <- rasterize(sub_poly2, sub_rasOut2, field = "area")

            DT <- data.table(cells = 1:ncell(sub_rasOutArea), area = as.vector(sub_rasOutArea[]))
            DT <- DT[complete.cases(DT)][cells %in% cellIDs]
            DT <- DT[order(area, decreasing = TRUE)]
            cellIDs <- DT[1:pixToConvert2, cells]
          }

          ## we may have exhausted areas to fill outside the holes
          ## so we fill holes (preferentially the smaller ones)
          if (length(cellIDs) < 1) {
            sub_poly <- as.polygons(sub_rasOut) |> disagg()
            sub_poly_holes <- fillHoles(sub_poly, inverse = TRUE) |> disagg()
            sub_poly_holes$area <- expanse(sub_poly_holes)

            sub_rasHolesArea <- rasterize(sub_poly_holes, sub_rasOut, field = "area")
            ## there may be patches inside holes (like islands) that need
            ## to be masked out, as they are ignored by fillHoles
            sub_rasHolesArea <- mask(sub_rasHolesArea, sub_rasOut, inverse = TRUE)

            DT <- data.table(cells = 1:ncell(sub_rasHolesArea), area = as.vector(sub_rasHolesArea[]))
            DT <- DT[complete.cases(DT)]
            setorder(DT, area, cells)

            cellIDs <- DT[1:pixToConvert2, cells]
          }

          sub_rasOut[cellIDs] <- 1L

          pixToConvert2 <- pixToConvert2 - length(cellIDs)
        }
      } else {
        ## there may be no available pixels, in which case permafrost can be assigned
        ## starting in a random polygon within unsuitable areas (hence equal prob below)
        spreadProb <- subst(sub_ras, c(NA, 0L), c(0L, 1L))

        sub_rasOut <- assignPresences(assignProb = spreadProb,
                                      landscape = sub_ras,
                                      pixToConvert = pixToConvert,
                                      probWeight = 1, numStartsDenom = 10)

      }
    }

    if (sum(!is.na(sub_rasOut[])) < permpercentPix) {
      warning(paste("Couldn't assign permafrost to enough pixels. id:", id))
    }
    message(cyan("Done!"))
  }

  if (saveOut) {
    tmpFile <- paste0("permafrost_polyID", id, ".tif")
    if (!is.null(saveDir)) {
      if (!dir.exists(saveDir)) dir.create(saveDir, showWarnings = FALSE)
      tmpFile <- file.path(saveDir, basename(tmpFile))
    }
    writeRaster(sub_rasOut, tmpFile, overwrite = TRUE)
    return(tmpFile)
  } else {
    return(sub_rasOut)
  }
}

#' Assign presences to patches based on a raster or presence
#'   probabilities
#'
#' @param assignProb a SpatRaster of presence probabilities
#' @param landscape a SpatRaster of the entire lanscape with NAs outside
#'   patches
#' @param pixToConvert numeric. Number of pixels to convert to presence
#'   across `landscape`. If `NULL`, `round(sum(!is.na(landscape[]))/2)`
#'   pixels will be converted to presences.
#' @param probWeight numeric. If `pixToConvert` cannot be reached
#'  `assignProb` will be weighted using `assignProb^probWeight`
#'  and the algorothm will try again.
#' @param numStartsDenom integer. Used to calculate the number of starting pixels
#'  to assign presences (at each try; as `pixToConvert/numStartsDenom`)
#'
#' @details This function attempts to iteratively assign
#'  presences starting in `pixToConvert/numStartsDenom` pixels
#'  sampled from areas with high probabilities in `assignProb` (weighted if
#'  `probWeight != 1`).
#'  For each starting pixel, it assigns a maximum of `numStartsDenom * probWeight`
#'  pixels. If not enough pixels are assigned presences (i.e. the number of
#'  presences does not reach `pixToConvert`), the function tries again
#'  after increasing `assignProb` by a factor of `1.5` (with values > 1 capped at 1).
#'  If too many pixels are assigned presences, pixels are sampled according to
#'  `assignProb ^ 10` (where `assignProb` corresponds to the original supplied values.)
#'
#' @return a SpatRaster
#' @importFrom SpaDES.tools spread
#' @importFrom raster raster
#' @importFrom terra rast
#' @importFrom data.table data.table
#'
#' @export
assignPresences <- function(assignProb, landscape, pixToConvert = NULL, probWeight = 1,
                            numStartsDenom = 10) {
  if (is.null(pixToConvert)) {
    pixToConvert <- round(sum(!is.na(landscape[]))/2)
  }

  if (probWeight < 0.5 || probWeight > 7)
    stop("probWeight must be between 0.5 and 7")

  ## save original probabilites for later
  assignProbOrig <- as.vector(assignProb[])

  ## if the mean is too high, then bring it down to 0.35 to avoid creating
  ## square patches
  meanSP <- mean(assignProbOrig, na.rm = TRUE)
  if (meanSP > 0.35)
    assignProb <- assignProb / meanSP * 0.35

  convertedPix <- 0

  while (convertedPix < pixToConvert) {
    ## exponentiate probabilities to provide more weight to pixels further from edges
    assignProbEx <- assignProb^probWeight
    # terra::plot(assignProbEx, col = viridis::inferno(100))

    ## try to spread in from many focal pixels
    numStarts <- ceiling(pixToConvert/numStartsDenom)
    startPoints <- sample(1:ncell(assignProbEx), size = numStarts, prob = assignProbEx[])

    # temp <- assignProbEx
    # temp[] <- NA_integer_
    # temp[startPoints] <- 1L
    # terra::plot(assignProbEx, col = viridis::inferno(100))
    # terra::plot(temp, add = TRUE, col = "blue")

    outRas <- SpaDES.tools::spread(landscape = raster(landscape),
                                   loci = startPoints,
                                   assignProb = raster(assignProbEx),
                                   maxSize = numStartsDenom * probWeight)
    outRas <- rast(outRas) |>
      mask(mask = landscape)
    convertedPix <- sum(!is.na(outRas[]))

    ## increase spread probabilities in case we need to try again
    assignProb <- min(assignProb * 1.5, 1)
  }

  ## if we spread too much remove pixels that are closest to edges
  if (convertedPix > pixToConvert) {
    ## "convert" to distances
    DT <- data.table(cells = 1:ncell(outRas), dists = assignProbOrig)
    DT <- DT[dists > 0 & !is.na(dists)][order(dists, decreasing = TRUE)]

    # pixToRm <- DT[(pixToConvert + 1):nrow(DT), cells] ## creates concave patches
    pixToKeep <- sample(DT$cells, pixToConvert, prob = DT$dists^10)  ## creates a little more noise

    outRas[] <- NA_integer_
    outRas[pixToKeep] <- 1L
  }
  convertedPix <- sum(!is.na(outRas[]))
  message(convertedPix, " were converted")
  return(outRas)
}
