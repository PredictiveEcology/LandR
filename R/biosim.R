utils::globalVariables(c(
  "Year"
))

#' Install \pkg{BioSIM} to retrieve climate and other projections using \code{BioSIM}
#'
#' @inheritParams utils::install.packages
#'
#' @export
#' @importFrom utils install.packages
installBioSIM <- function(lib) {
  install.packages("https://sourceforge.net/projects/repiceasource/files/latest",
                   lib, repos = NULL,  type = "source")
  install.packages("https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest",
                   lib, repos = NULL,  type = "source")
}

#' Extract point locations from DEM raster to pass to \code{BioSIM} functions
#'
#' @param x A digital elevation model (DEM) \code{RasterLayer}.
#'
#' @return \code{data.table} with columns \code{Name}, \code{Long}, \code{Lat}, \code{Elev}.
#' @export
#' @importFrom data.table data.table setnames
#' @importFrom raster crs xyFromCell
#' @importFrom sf st_as_sf st_coordinates
#' @importFrom sp SpatialPoints
#' @seealso prepInputsCanDEM
BioSIM_extractPoints <- function(x) {
  nonNA <- which(!is.na(x[]))
  xy <- xyFromCell(x, cell = nonNA)
  spxy <- SpatialPoints(xy)
  crs(spxy) <- crs(x)
  sfxy <- sf::st_as_sf(spxy)
  sfxy <- sf::st_transform(sfxy, crs = 4326)

  dt <- data.table(Name = paste0("ID", 1:NROW(sfxy)), st_coordinates(sfxy), Elev = x[nonNA])
  setnames(dt, "X", "Long")
  setnames(dt, "Y", "Lat")
  dt
}

#' Get annual historic and projected wind maps from \code{BioSIM}
#'
#' @param dem \code{RasterLayer} of elevation data (m).
#' @param years numeric vector corresponding to the years to retrieve.
#' @param climModel climate model to use. one of \code{"GCM4"} or \code{"RCM4"}.
#' @param rcp RCP scenario to use. one of \code{"RCP45"} or \code{"RCP85"}.
#'
#' @return \code{RasterStack}
#' @export
#' @importFrom data.table setDT
#' @importFrom raster cellFromXY raster stack
#' @importFrom reproducible Cache
#' @importFrom sf st_as_sf st_coordinates st_crs st_transform
#' @importFrom sp CRS SpatialPoints
BioSIM_getWind <- function(dem, years, climModel = "GCM4", rcp = "RCP45") {
  if (requireNamespace("BioSIM", quietly = TRUE)) {
    locations <- BioSIM_extractPoints(dem)

    # Do call to BioSIM using "ClimaticWind_Annual"
    windModel <- Cache(BioSIM::getModelList)[16] ## TODO: until this gets fixed in J4R, need to init java server here
    stWind <- system.time({
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

    message("Fetched ", NROW(wind), " locations in ", stWind[3], "s.")

    # Make RasterStack
    setDT(wind)
    windStk <- stack(raster(dem))
    for (yr in unique(wind$Year)) {
      yrChar <- paste0("X", yr)
      windYr <- wind[Year == yr]

      # Convert BioSIM data to Vector dataset
      spWind <- SpatialPoints(wind[Year == yr, c("Longitude", "Latitude")],
                              proj4string = CRS("+init=epsg:4326")) %>%
        sf::st_as_sf(.) %>%
        sf::st_transform(., crs = sf::st_crs(dem))
      cells <- cellFromXY(dem, sf::st_coordinates(spWind))
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
