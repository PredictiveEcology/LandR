utils::globalVariables(c(
  "Year"
))

#' Install \pkg{BioSIM} to retrieve climate and other projections using `BioSIM`
#'
#' @inheritParams utils::install.packages
#'
#' @export
installBioSIM <- function(lib) {
  install.packages("https://sourceforge.net/projects/repiceasource/files/latest",
                   lib, repos = NULL,  type = "source")
  install.packages("https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest",
                   lib, repos = NULL,  type = "source")
}

#' Extract point locations from DEM raster to pass to `BioSIM` functions
#'
#' @param x A digital elevation model (DEM) `RasterLayer`.
#'
#' @return `data.table` with columns `Name`, `Long`, `Lat`, `Elev`.
#' @export
#' @seealso [prepInputsCanDEM()]
BioSIM_extractPoints <- function(x) {
  nonNA <- which(!is.na(x[]))
  xy <- xyFromCell(x, cell = nonNA)
  spxy <- SpatialPoints(xy, proj4string = crs(x))
  sfxy <- sf::st_as_sf(spxy)
  sfxy <- sf::st_transform(sfxy, crs = 4326)

  dt <- data.table(Name = paste0("ID", 1:NROW(sfxy)), st_coordinates(sfxy), Elev = x[nonNA])
  setnames(dt, "X", "Long")
  setnames(dt, "Y", "Lat")
  dt
}

#' Get annual historic and projected wind maps from `BioSIM`
#'
#' @param dem `RasterLayer` of elevation data (m).
#' @param years numeric vector corresponding to the years to retrieve.
#' @param climModel climate model to use. one of `"GCM4"` or `"RCM4"`.
#' @param rcp RCP scenario to use. one of `"RCP45"` or `"RCP85"`.
#'
#' @return `RasterStack`
#' @export
BioSIM_getWindAnnual <- function(dem, years, climModel = "GCM4", rcp = "RCP45") {
  if (requireNamespace("BioSIM", quietly = TRUE)) {
    locations <- BioSIM_extractPoints(dem)

    # Do call to BioSIM using "ClimaticWind_Annual"
    windModel <- Cache(BioSIM::getModelList)[16] ## TODO: until this gets fixed in J4R, need to init java server here
    st <- system.time({
      ## TODO: need to split, apply, and recombine when nrow(locations) > 5000
      wind <- Cache(
        BioSIM::getModelOutput,
        fromYr = years[1],
        toYr = rev(years)[1],
        id = locations$Name,
        latDeg = locations$Lat,
        longDeg = locations$Long,
        elevM = locations$Elev,
        modelName = windModel,
        rcp = rcp,
        climModel = climModel
      )
    })

    message("Fetched ", NROW(wind), " locations in ", st[3], "s.")

    # Make RasterStack
    setDT(wind)
    windStk <- stack(raster(dem))
    for (yr in unique(wind$Year)) {
      yrChar <- paste0("X", yr)
      windYr <- wind[Year == yr]

      # Convert BioSIM data to Vector dataset
      sfWind <- SpatialPoints(wind[Year == yr, c("Longitude", "Latitude")],
                              proj4string = CRS("+init=epsg:4326")) |>
        sf::st_as_sf() |>
        sf::st_transform(crs = sf::st_crs(dem))
      cells <- cellFromXY(dem, sf::st_coordinates(sfWind))
      cols <- grep("^W[[:digit:]]", colnames(windYr), value = TRUE)

      # Convert to single main direction
      dirs <- (apply(windYr[, ..cols], 1, which.max) - 1) * 10

      # Convert to Raster
      windYrRas <- raster(dem)
      windYrRas[cells] <- dirs
      windStk[[yrChar]] <- windYrRas
    }
    windStk
  } else {
    stop("Package BioSIM not installed. Use `installBioSIM()` to install it.")
  }
}

#' Get monthly historic and projected wind maps from `BioSIM`
#'
#' @param dem `RasterLayer` of elevation data (m).
#' @param years numeric vector corresponding to the years to retrieve.
#' @param months numeric vector corresponding to the months to retrieve
#'               (e.g., `6:8` for June through August).
#' @param climModel climate model to use. one of `"GCM4"` or `"RCM4"`.
#' @param rcp RCP scenario to use. one of `"RCP45"` or `"RCP85"`.
#'
#' @return `RasterStack`
#'
#' @export
BioSIM_getWindMonthly <- function(dem, years, months, climModel = "GCM4", rcp = "RCP45") {
  if (requireNamespace("BioSIM", quietly = TRUE)) {
    locations <- BioSIM_extractPoints(dem)

    # Do call to BioSIM using "ClimaticWind_Annual"
    windModel <- Cache(BioSIM::getModelList)[16] ## TODO: until this gets fixed in J4R, need to init java server here
    st <- system.time({
      ## TODO: need to split, apply, and recombine when nrow(locations) > 5000
      wind <- Cache(
        BioSIM::getModelOutput,
        fromYr = years[1],
        toYr = rev(years)[1],
        id = locations$Name,
        latDeg = locations$Lat,
        longDeg = locations$Long,
        elevM = locations$Elev,
        modelName = windModel,
        rcp = rcp,
        climModel = climModel
      )
    })

    message("Fetched ", NROW(wind), " locations in ", st[3], "s.")

    # Make RasterStack
    browser() ## TODO: pull in Eliot's code from mpbClimateData

  } else {
    stop("Package BioSIM not installed. Use `installBioSIM()` to install it.")
  }
}

#' Get annual historic and projected MPB climate suitability maps from `BioSIM`
#'
#' Raster stacks for all 9 MPB climate indices. See `BioSIM::getModelHelp("MPB_SLR")`.
#'
#' @param dem `RasterLayer` of elevation data (m).
#' @param years numeric vector corresponding to the years to retrieve.
#' @param SLR character. Specifies which climate suitability index to extract.
#'            Currently, one of `"S"`, `"L"`, `"R"`, or `"G"`,
#'            corresponding to Safranyik-P3P4, Logan-2b, Régnière Cold Tolerance Survival, or
#'            their Geometric product (S\*L\*R), respectively.
#' @param climModel climate model to use. one of `"GCM4"` or `"RCM4"`.
#' @param rcp RCP scenario to use. one of `"RCP45"` or `"RCP85"`.
#'
#' @note Although the BioSIM MPB_SLR model provides several other indices
#'       (see `BioSIM::getModelHelp("MPB_SLR")`), only 4 are currently used here.
#'
#' @return `RasterStack`
#' @export
BioSIM_getMPBSLR <- function(dem, years, SLR = "R", climModel = "GCM4", rcp = "RCP45") {
  if (requireNamespace("BioSIM", quietly = TRUE)) {
    SLR2use <- switch(SLR,
                      S = "Safranyik_p_34",
                      L = "Logan_P_2b",
                      R = "CT_Survival",
                      G = "Geo_prod_pL2b_pS34_pC",
                      stop("SLR must be one of 'S', 'L', 'R', 'G'."))

    locations <- BioSIM_extractPoints(dem)

    # Do call to BioSIM using "MPB_SLR"
    mpbSLRmodel <- Cache(BioSIM::getModelList)[46] ## TODO: until this gets fixed in J4R, need to init java server here
    st <- system.time({
      ## TODO: need to split, apply, and recombine when nrow(locations) > 5000
      slr <- lapply(years, function(yr) { ## TODO: use future_lapply?
        Cache(
          BioSIM::getModelOutput,
          fromYr = yr - 1,
          toYr = yr,
          id = locations$Name,
          latDeg = locations$Lat,
          longDeg = locations$Long,
          elevM = locations$Elev,
          modelName = mpbSLRmodel,
          rep = 1, ## TODO: how many?
          rcp = rcp,
          climModel = climModel
        )
        setDT(slr)
      }) |>
        rbindlist()
    })

    message("Fetched ", NROW(slr), " locations in ", st[3], "s.")

    # Make RasterStack
    slrStk <- stack(raster(dem))
    for (yr in unique(slr$Year)) {
      yrChar <- paste0("X", yr)
      slrYr <- slr[Year == yr]

      colID <- which(colnames(slrYr) == SLR2use)

      # Convert BioSIM data to Vector dataset
      sfSLR <- SpatialPoints(slr[Year == yr, c("Longitude", "Latitude")],
                             proj4string = CRS("+init=epsg:4326")) |>
        sf::st_as_sf() |>
        sf::st_transform(crs = sf::st_crs(dem))
      cells <- cellFromXY(dem, sf::st_coordinates(sfSLR))

      # Convert to Raster
      slrYrRas <- raster(dem)
      slrYrRas[cells] <- slrYr[[colID]]
      slrStk[[yrChar]] <- slrYrRas
    }
    slrStk
  } else {
    stop("Package BioSIM not installed. Use `installBioSIM()` to install it.")
  }
}
