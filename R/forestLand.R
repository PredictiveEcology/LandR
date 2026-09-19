## NFIS publishes the FAO forest layer for these years only
.faoForestYears <- c(2019L, 2022L)

## Land-cover classes that are treed. 81 (wetland-treed) does not occur in the SCANFI
## product, which has no wetland classes at all, so the same vector serves both sources.
.treedLCCClasses <- c(81L, 210L, 220L, 230L)

#' Is a pixel forest land?
#'
#' Forest land is a land *use*: ground that grows trees, including ground that has none
#' right now because it burned or was harvested. The land-cover maps answer a different
#' question -- what covered the ground in one year -- so a recent burn is "shrubs" there.
#' Telling the temporarily open forest apart from land that is permanently open (a bog, a
#' rock barren) needs a second source, and there are two to choose from.
#'
#' `"fao"` uses the NFIS FAO forest layer: code 1 is treed in the layer's own year, code 2
#' is not treed then but with a fire or harvest recorded since 1984. It is one small layer.
#' Its limit is its fixed year: a stand that was open in 1995 and had grown back by then is
#' code 1, which the old rule (code 2 only) discarded, and openings whose disturbance
#' predates 1984 are code 0.
#'
#' `"lccYears"` calls a pixel forest land if it is treed in *any* of several land-cover
#' years. That sees the recovered stands and the pre-1984 openings, at the cost of reading
#' one land-cover layer per year.
#'
#' `"both"` takes the union, and is the default.
#'
#' @param lccList A list of land-cover `SpatRaster`s, all aligned. May be empty.
#' @param faoRas An FAO forest `SpatRaster` aligned with `lccList`, or `NULL`.
#' @param treedClasses Land-cover classes that count as treed.
#' @param faoForestCodes FAO codes that count as forest land. Both 1 and 2 are forest land;
#'   code 2 alone means "still an opening in the FAO layer's year".
#'
#' @return A `SpatRaster` of 1 (forest land) and 0, with no `NA`s.
#'
#' @export
forestLandMask <- function(lccList = list(), faoRas = NULL,
                           treedClasses = .treedLCCClasses,
                           faoForestCodes = c(1L, 2L)) {
  masks <- lapply(lccList, function(r) .naTo0(terra::`%in%`(r, treedClasses)))
  if (!is.null(faoRas)) {
    masks <- c(masks, list(.naTo0(terra::`%in%`(faoRas, faoForestCodes))))
  }
  if (!length(masks)) {
    stop("forestLandMask() needs at least one of `lccList` or `faoRas`.")
  }
  ## `|` returns a logical layer while a single `%in%` returns 1/0, so flatten either to 1/0
  .naTo0(Reduce(function(a, b) a | b, masks))
}

#' Label the forest land that has no trees this year
#'
#' Every non-treed class is eligible by default, exposed barren land and rock included: a
#' pixel is only relabelled where the evidence says it is forest land, and a rock or barren
#' pixel that the record shows as treed in other years is exposed by something temporary,
#' such as a severe fire. Water and snow/ice are in that set too, so a pixel that was treed
#' in an early year and is a reservoir later will be called forest land; pass
#' `convertibleClasses` to narrow the set if that matters for a given study area.
#'
#' @param lcc A land-cover `SpatRaster` for the year being prepared.
#' @param forestLand A mask from [forestLandMask()], aligned with `lcc`.
#' @param treedClasses The classes that are already treed, and so never relabelled.
#' @param convertibleClasses The classes that may become `disturbedCode`, or `NULL`
#'   (default) for every class that is not in `treedClasses`. `NA` pixels are never
#'   relabelled.
#' @param disturbedCode The class given to forest land that is not treed in this year.
#' @param filename Optional file to write to.
#'
#' @return A `SpatRaster` like `lcc`, with `disturbedCode` where it applies.
#'
#' @keywords internal
.applyForestLand <- function(lcc, forestLand, treedClasses = .treedLCCClasses,
                             convertibleClasses = NULL, disturbedCode = 240L,
                             filename = NULL) {
  isConvertible <- if (is.null(convertibleClasses)) {
    (!.naTo0(terra::`%in%`(lcc, treedClasses))) & !is.na(lcc)
  } else {
    .naTo0(terra::`%in%`(lcc, convertibleClasses))
  }
  if (is.null(filename)) {
    filename <- tempfile(fileext = ".tif")
  }
  terra::ifel(
    forestLand & isConvertible,
    disturbedCode,
    lcc,
    overwrite = TRUE,
    filename = filename,
    wopt = list(datatype = "INT1U", gdal = c("COMPRESS=ZSTD", "TILED=YES"))
  )
}

## terra's `%in%` returns NA for NA cells, and an NA in `ifel`'s condition blanks the
## pixel, so masks are flattened to 0/1 before they are combined or used.
.naTo0 <- function(x) terra::ifel(is.na(x), 0, x)

#' Obtain the NFIS FAO forest layer
#'
#' @param year One of `r .faoForestYears`.
#' @param to Passed to `prepInputs`, to align the layer with the land cover.
#' @param destinationPath Passed to `prepInputs`.
#' @param method Resampling method passed to `prepInputs`.
#' @param ... Other arguments passed to `prepInputs`, e.g. `maskTo = NA` to align without masking.
#'
#' @return a `SpatRaster` of FAO forest codes: 0 non-forest, 1 forest, 2 forest land whose
#'   trees were removed by fire or harvest since 1984.
#'
#' @export
prepInputs_FAO_forest <- function(year = 2022, to = NULL, destinationPath = NULL,
                                  method = "near", ...) {
  if (!(year %in% .faoForestYears)) {
    stop("The FAO forest layer is published for ", paste(.faoForestYears, collapse = " and "),
         " only; got ", year, ".")
  }
  prepInputs(
    url = paste0("https://opendata.nfis.org/downloads/forest_change/CA_FAO_forest_", year, ".zip"),
    method = method,
    destinationPath = destinationPath,
    to = to,
    ...
  )
}

#' The default land-cover years used to decide forest land
#'
#' Roughly one year per decade plus the most recent, which is enough because a stand that
#' recovers stays treed for decades. Every returned year is one the source publishes.
#'
#' @param availableYears The years the source has.
#'
#' @return An integer vector of years.
#'
#' @keywords internal
.defaultForestLandYears <- function(availableYears) {
  availableYears <- sort(unique(as.integer(availableYears)))
  wanted <- c(1985L, 1995L, 2005L, 2015L, max(availableYears))
  unique(vapply(wanted, function(y) {
    availableYears[which.min(abs(availableYears - y))]
  }, integer(1)))
}

#' Build the forest-land mask for one of the land-cover sources
#'
#' @param lcc The land-cover `SpatRaster` for the year being prepared; everything is
#'   aligned to it.
#' @param forestLandFrom One of `"both"`, `"fao"`, `"lccYears"`. See [forestLandMask()].
#' @param forestLandYears The land-cover years scanned when `forestLandFrom` is
#'   `"both"` or `"lccYears"`.
#' @param faoYear The year of the FAO layer, when one is used.
#' @param lccFor A function of one argument (a year) returning that year's land cover,
#'   aligned with `lcc` but not masked by it (`maskTo = NA`). Each source passes its own.
#' @param lccSource A name for the land-cover source `lccFor` reads, e.g. `"SCANFI V2"`. It is
#'   part of the cache key of each year's layer, so it must differ between sources.
#' @param treedClasses Passed to [forestLandMask()].
#' @param destinationPath Passed to `prepInputs`.
#' @param resampleMethod Passed to `prepInputs`.
#'
#' @details
#' Each forest-land input is `Cache()`d with `useCache = "always"`, so a caller that prepares
#' several years of one study area, as `fireSense` does, prepares each input once; `"always"`
#' holds even when Cache is otherwise off, e.g. under `spades.useCache = "eventsOnly"`, and when
#' the call is nested in a `Cache()` that is. The inputs are aligned with `lcc` but not masked by
#' it, so they depend on its geometry and not its values; the key is that geometry (crs, extent,
#' dimensions), never `lcc` itself. `lcc`'s own `NA`s are excluded by [.applyForestLand()].
#' There is no switch; delete the entries (`reproducible::clearCache()`) to recompute.
#'
#' @return A `SpatRaster` mask, as [forestLandMask()].
#'
#' @keywords internal
.forestLandFor <- function(lcc, forestLandFrom = c("both", "fao", "lccYears"),
                           forestLandYears, faoYear, lccFor, lccSource,
                           treedClasses = .treedLCCClasses,
                           destinationPath = NULL, resampleMethod = "near") {
  forestLandFrom <- match.arg(forestLandFrom)
  ## the key: `lcc`'s geometry, never its values, which differ between years
  geometry <- list(crs = terra::crs(lcc), ext = as.vector(terra::ext(lcc)), dim = dim(lcc)[1:2])

  lccList <- list()
  if (forestLandFrom %in% c("both", "lccYears")) {
    lccList <- lapply(forestLandYears, function(y) {
      message("  ... forest land: land cover for ", y)
      Cache(lccFor(y), useCache = "always", .functionName = "forestLand_landCover",
            .cacheExtra = list(geometry, lccSource, resampleMethod))
    })
  }

  faoRas <- NULL
  if (forestLandFrom %in% c("both", "fao")) {
    message("  ... forest land: FAO forest ", faoYear)
    faoRas <- Cache(
      prepInputs_FAO_forest(year = faoYear, to = lcc, maskTo = NA,
                            destinationPath = destinationPath, method = resampleMethod),
      useCache = "always", omitArgs = "to", .functionName = "forestLand_FAO",
      .cacheExtra = geometry
    )
  }

  forestLandMask(lccList = lccList, faoRas = faoRas, treedClasses = treedClasses)
}
