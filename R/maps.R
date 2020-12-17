utils::globalVariables(c(
  ".", "..pgdAndScAndLeading", ":=", "B", "HQ", "leading", "LQ", "mixed", "N",
  "pixelGroup", "pure", "speciesCode", "speciesGroupB", "speciesProportion", "SPP",
  "totalB", "totalcover", "Type"
))

#' Define flammability map
#'
#' @param LandCoverClassifiedMap A \code{Raster} that represents land cover
#' (e.g., Land Cover Classified map from 2005 or 2010 from the Canadian Forest Service).
#'
#' @param nonFlammClasses numeric vector defining which classes in \code{LandCoverClassifiedMap}.
#'
#' @param mask A raster to use as a mask (see \code{\link[raster]{mask}}).
#'
#' @param filename2 See \code{\link[reproducible]{postProcess}}. Default \code{NULL}.
#'
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom quickPlot setColors<-
#' @importFrom raster maxValue minValue ratify reclassify writeRaster
defineFlammable <- function(LandCoverClassifiedMap = NULL,
                            nonFlammClasses = c(0L, 25L, 30L, 33L,  36L, 37L, 38L, 39L),
                            mask = NULL, filename2 = NULL) {
  if (!is.null(mask))
    if (!is(mask, "Raster")) stop("mask must be a raster layer")
  if (!is(LandCoverClassifiedMap, "RasterLayer")) {
    stop("Need a classified land cover map. Currently only accepts 'LCC2005'")
  }
  if (!is.integer(LandCoverClassifiedMap[]))
    stop("LandCoverClassifiedMap must be an integer")
  if (is.null(nonFlammClasses))
    stop("Need nonFlammClasses, which are the classes that cannot burn in",
         "the LandCoverClassifiedMap")

  oldClass <- minValue(LandCoverClassifiedMap):maxValue(LandCoverClassifiedMap)
  newClass <- ifelse(oldClass %in% nonFlammClasses, 0L, 1L) ## NOTE: 0 codes for NON-flammable
  flammableTable <- cbind(oldClass, newClass)
  rstFlammable <- ratify(reclassify(LandCoverClassifiedMap, flammableTable))
  if (!is.null(filename2))
    rstFlammable <- writeRaster(rstFlammable, filename = filename2, overwrite = TRUE)

  setColors(rstFlammable, n = 2) <- colorRampPalette(c("blue", "red"))(2)
  if (!is.null(mask)) rstFlammable[is.na(mask[])] <- NA_integer_

  rstFlammable[] <- as.integer(rstFlammable[])
  rstFlammable
}

#' Simple \code{prepInputs} for LCC2005 or LCC2010 data
#'
#' A wrapper around \code{prepInputs} for the Canadian Land Cover Classification product(s).
#'
#' @inheritParams reproducible::cropInputs
#' @inheritParams reproducible::postProcess
#' @inheritParams reproducible::prepInputs
#'
#' @param year Numeric, either 2005 or 2010
#'
#' @export
#' @importFrom reproducible asPath prepInputs
prepInputsLCC <- function(year = 2005,
                          destinationPath = asPath("."),
                          studyArea = NULL,
                          rasterToMatch = NULL,
                          filename2 = NULL, ...) {

  dots <- list(...)
  if (is.null(dots$url)) {
    if (identical(as.integer(year), 2005L)) {
      #url <- "https://drive.google.com/file/d/1g9jr0VrQxqxGjZ4ckF6ZkSMP-zuYzHQC/view?usp=sharing"
      url <- paste0("ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/",
                    "LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip")
      filename <- asPath("LCC2005_V1_4a.tif")
      archive <- asPath("LandCoverOfCanada2005_V1_4.zip")
    } else {
      if (identical(as.integer(year), 2010L)) {
        url <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/",
                      "Land-cover_Couverture-du-sol/canada-landcover_canada-couverture-du-sol/",
                      "CanadaLandcover2010.zip")
        filename <- asPath("CAN_LC_2010_CAL.tif")
        archive <- asPath("CanadaLandcover2010.zip")
      } else {
        stop("Other LCC covers don't exist yet.")
      }
    }
  }

  Cache(prepInputs, targetFile = filename,
        archive = archive,
        url = url,
        destinationPath = asPath(destinationPath),
        studyArea = studyArea,
        rasterToMatch = rasterToMatch,
        method = "ngb",
        datatype = "INT2U",
        filename2 = filename2, ...)
}

#' Make a vegetation type map from a stack of species abundances
#'
#' @description
#' \code{makeVegTypeMap} is a wrapper around \code{vegTypeMapGenerator}
#' that works from a species stack of percent cover. These do not have
#' to sum to 100%
#'
#' @param speciesStack A \code{RasterStack} of species abundances.
#'                     This must be one \code{RasterLayer} per species.
#' @param vegLeadingProportion See \code{vegTypeMapGenerator}.
#' @param mixed Deprecated. See \code{mixedType} argument to \code{vegTypeMapGenerator}.
#' @param ... Other arguments passed to \code{vegTypeMapGenerator}, i.e.,
#'   \code{vegLeadingProportion}, \code{mixedType}, \code{sppEquiv},
#'   \code{sppEquivCol}, \code{colors}, \code{pixelGroupColName}, and \code{doAssertion}
#'
#' @return A factor raster
#'
#' @export
#' @importFrom quickPlot numLayers
#' @importFrom raster levels maxValue raster
#' @importFrom SpaDES.tools inRange
#' @rdname LandR-deprecated
makeVegTypeMap <- function(speciesStack, vegLeadingProportion, mixed, ...) {
  .Deprecated("vegTypeMapGenerator")

  if (isTRUE(mixed)) mixed <- 2

  vegTypeMapGenerator(x = speciesStack, vegLeadingProportion = vegLeadingProportion,
                      mixedType = as.numeric(mixed), ...)
}

#' Generate vegetation type map
#'
#' @param x Either a \code{cohortData} object or a \code{speciesCover} \code{RasterStack}
#'
#' @template pixelGroupMap
#'
#' @param vegLeadingProportion Numeric between 0-1, determining the relative biomass
#'               threshold a species needs to pass to be considered "leading".
#'
#' @param mixedType An integer defining whether mixed stands are of any kind of species
#'                  admixture (1), or only when deciduous mixed with conifer (2).
#'                  Defaults to 2.
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param colors A named vector of colour codes. The names MUST match the names of species
#'               in \code{cohortData$speciesCode}, plus an optional "Mixed" colour.
#'
#' @param pixelGroupColName Name of the column in \code{pixelGroup} to use.
#'
#' @template doAssertion
#'
#' @param ... Additional arguments.
#'
#' @author Eliot McIntire, Ceres Barros, Alex Chubaty
#' @export
#' @importFrom data.table copy data.table setkey setorderv
#' @importFrom pemisc factorValues2
#' @importFrom raster getValues projection projection<- setValues
#' @importFrom SpaDES.tools rasterizeReduced
#' @importFrom utils data
#' @rdname vegTypeMapGenerator
vegTypeMapGenerator <- function(x, ...) {
  UseMethod("vegTypeMapGenerator")
}

#' @export
#' @rdname vegTypeMapGenerator
vegTypeMapGenerator.RasterStack <- function(x, ..., doAssertion = getOption("LandR.doAssertion", TRUE)) {
  pixelTable <- suppressMessages(makePixelTable(x, printSummary = FALSE, doAssertion = doAssertion))
  cohortTable <- suppressMessages(.createCohortData(pixelTable, rescale = FALSE, doAssertion = doAssertion))
  cohortTable <- cohortTable[cover > 0]
  pixelGroupMap <- raster(x)
  pixelGroupMap[pixelTable[["pixelIndex"]]] <- pixelTable[["initialEcoregionCode"]]
  vegTypeMap <- vegTypeMapGenerator(x = cohortTable,
                                    pixelGroupMap = pixelGroupMap,
                                    pixelGroupColName = "initialEcoregionCode",
                                    doAssertion = doAssertion,
                                    ...)

  if (FALSE) { # This is the old version -- Eliot & Alex July 11, 2019
    sumVegPct <- sum(speciesStack) ## TODO: how is the sum >100 ?

    if (isTRUE(mixed)) {
      ## create "mixed" layer, which is given a value slightly higher than any other layer,
      ## if it is deemed a mixed pixel.
      ## All layers must be below vegLeadingProportion to be called Mixed.
      ## This check turns stack to binary: 1 if < vegLeadingProportion; 0 if more than.
      ## Then, sum should be numLayers of all are below vegLeadingProportion
      whMixed <- which(sum(speciesStack < (100 * vegLeadingProportion))[] == numLayers(speciesStack))
      MixedRas <- speciesStack[[1]]
      MixedRas[!is.na(speciesStack[[1]][])] <- 0
      MixedRas[whMixed] <- max(maxValue(speciesStack)) * 1.01

      speciesStack$Mixed <- MixedRas
    }

    a <- speciesStack[]
    nas <- is.na(a[, 1])
    maxes <- apply(a[!nas, ], 1, function(x) {
      whMax <- which(x == max(x, na.rm = TRUE))
      if (length(whMax) > 1) {
        whMax <- sample(whMax, size = 1)
      }
      return(whMax)
    })

    vegTypeMap <- raster(speciesStack[[1]])

    vegTypeMap[!nas] <- maxes

    layerNames <- names(speciesStack)
    names(layerNames) <- layerNames
    levels(vegTypeMap) <- data.frame(ID = seq(layerNames), Species = names(layerNames),
                                     stringsAsFactors = TRUE)
    vegTypeMap

  }
  return(vegTypeMap)
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom SpaDES.tools inRange
#' @rdname vegTypeMapGenerator
vegTypeMapGenerator.data.table <- function(x, pixelGroupMap, vegLeadingProportion = 0.8,
                                           mixedType = 2, sppEquiv = NULL, sppEquivCol, colors,
                                           pixelGroupColName = "pixelGroup",
                                           doAssertion = getOption("LandR.assertions", TRUE), ...) {
  if (!inRange(vegLeadingProportion, 0, 1))
    stop("vegLeadingProportion must be a proportion")

  nrowCohortData <- NROW(x)

  leadingBasedOn <- if ("B" %in% colnames(x)) {
    message("Using B to derive leading type")
    "B"
  } else if ("cover" %in% colnames(x)) {
    message("Using cover to derive leading type, as there is no B column")
    "cover"
  } else {
    stop("x must have either B or cover to determine leading species/type")
  }
  assert_that(nrowCohortData > 0)

  if (isTRUE(doAssertion))
    message("LandR::vegTypeMapGenerator: NROW(x) == ", nrowCohortData)

  if (mixedType == 2) {
    if (is.null(sppEquiv)) {
      sppEquiv <- get(data("sppEquivalencies_CA", package = "LandR", envir = environment()),
                      inherits = FALSE)

      # Find the sppEquivCol that best matches what you have in x
      sppEquivCol <- names(sort(sapply(sppEquiv, function(xx) sum(xx %in% unique(x$species))),
                                decreasing = TRUE)[1])
      message(paste0("Using mixedType == 2, but no sppEquiv provided. ",
                     "Attempting to use data('sppEquivalencies_CA', 'LandR') ",
                     "and sppEquivCol == '", sppEquivCol, "'"))
    }
  }

  ## use new vs old algorithm based on size of x. new one (2) is faster in most cases.
  ## enable assertions to view timings for each algorithm before deciding which to use.
  algo <- ifelse(nrowCohortData > 3.5e6, 1, 2)

  pgdAndSc <- c(pixelGroupColName, "speciesCode")
  pgdAndScAndLeading <- c(pgdAndSc, leadingBasedOn)
  totalOfLeadingBasedOn <- paste0("total", leadingBasedOn)
  speciesOfLeadingBasedOn <- paste0("speciesGroup", leadingBasedOn)
  if (algo == 1 || isTRUE(doAssertion)) {
    # slower -- older, but simpler Eliot June 5, 2019
    # 1. Find length of each pixelGroup -- don't include pixelGroups in "by" that have only 1 cohort: N = 1
    cohortData1 <- copy(x)
    systimePre1 <- Sys.time()
    pixelGroupData1 <- cohortData1[, list(N = .N), by = pixelGroupColName]

    # Calculate speciesProportion from cover or B
    pixelGroupData1 <- cohortData1[, ..pgdAndScAndLeading][pixelGroupData1, on = pixelGroupColName]
    set(pixelGroupData1, NULL, totalOfLeadingBasedOn, pixelGroupData1[[leadingBasedOn]])

    if (identical(leadingBasedOn, "cover")) {
      pixelGroupData1[N != 1, (totalOfLeadingBasedOn) := sum(cover, na.rm = TRUE), by = pixelGroupColName]
      pixelGroupData1 <- pixelGroupData1[, list(sum(cover, na.rm = TRUE), totalcover[1]), by = pgdAndSc]
    } else {
      pixelGroupData1[N != 1, (totalOfLeadingBasedOn) := sum(B, na.rm = TRUE), by = pixelGroupColName]
      pixelGroupData1 <- pixelGroupData1[, list(sum(B, na.rm = TRUE), totalB[1]), by = pgdAndSc]
    }
    setnames(pixelGroupData1, old = c("V1", "V2"), new = c(speciesOfLeadingBasedOn, totalOfLeadingBasedOn))

    set(pixelGroupData1, NULL, "speciesProportion", pixelGroupData1[[speciesOfLeadingBasedOn]] /
          pixelGroupData1[[totalOfLeadingBasedOn]])
    systimePost1 <- Sys.time()

    setorderv(pixelGroupData1, pixelGroupColName)
  }

  if (algo == 2 || isTRUE(doAssertion)) {
    # Replacement algorithm to calculate speciesProportion
    #  Logic is similar to above --
    #  1. sort by pixelGroup
    #  2. calculate N, use this to repeat itself (instead of a join above)
    #  3. calculate speciesProportion, noting to calculate with by only if N > 1, otherwise
    #     it is a simpler non-by calculation
    cohortData2 <- copy(x)
    systimePre2 <- Sys.time()
    setkeyv(cohortData2, pgdAndSc)
    # setorderv(x, pixelGroupColName)
    pixelGroupData2 <- cohortData2[, list(N = .N), by = pixelGroupColName]
    cohortData2 <- cohortData2[, ..pgdAndScAndLeading]

    N <- rep.int(pixelGroupData2$N, pixelGroupData2$N)
    wh1 <- N == 1
    set(cohortData2, which(wh1), totalOfLeadingBasedOn, cohortData2[[leadingBasedOn]][wh1])
    if (identical(leadingBasedOn, "cover")) {
      totalBNot1 <- cohortData2[!wh1, list(N = .N, totalcover = sum(cover, na.rm = TRUE)), by = pixelGroupColName]
    } else {
      totalBNot1 <- cohortData2[!wh1, list(N = .N, totalB = sum(B, na.rm = TRUE)), by = pixelGroupColName]
    }
    totalBNot1 <- rep.int(totalBNot1[[totalOfLeadingBasedOn]], totalBNot1[["N"]])
    set(cohortData2, which(!wh1), totalOfLeadingBasedOn, totalBNot1)

    b <- cohortData2[, list(N = .N), by = pgdAndSc]
    b <- rep.int(b[["N"]], b[["N"]])
    GT1 <- (b > 1)
    if (any(GT1)) {
      pixelGroupData2List <- list()
      cohortData2[GT1, speciesProportion := sum(B, na.rm = TRUE) / totalB[1], by = pgdAndSc]
      cohortData2[!GT1, speciesProportion := B / totalB]
      #pixelGroupData2List[[2]] <- cohortData2[!GT1]
      #pixelGroupData2 <- rbindlist(pixelGroupData2List)
    } else {
      # cols <- c(pixelGroupColName, "speciesCode", "speciesProportion")
      set(cohortData2, NULL, "speciesProportion", cohortData2[[leadingBasedOn]] /
            cohortData2[[totalOfLeadingBasedOn]])
      # pixelGroupData2[[NROW(pixelGroupData2) + 1]] <- cohortData2[!GT1, ..cols]
    }
    pixelGroupData2 <- cohortData2
    systimePost2 <- Sys.time()
  }

  if (isTRUE(doAssertion)) {
    ## slower -- older, but simpler Eliot June 5, 2019
    ## TODO: these algorithm tests should be deleted after a while. See date on prev line.
    if (!exists("oldAlgoVTM", envir = .pkgEnv)) .pkgEnv$oldAlgoVTM <- 0
    if (!exists("newAlgoVTM", envir = .pkgEnv)) .pkgEnv$newAlgoVTM <- 0
    .pkgEnv$oldAlgoVTM <- .pkgEnv$oldAlgoVTM + (systimePost1 - systimePre1)
    .pkgEnv$newAlgoVTM <- .pkgEnv$newAlgoVTM + (systimePost2 - systimePre2)
    message("LandR::vegTypeMapGenerator: new algo ", .pkgEnv$newAlgoVTM)
    message("LandR::vegTypeMapGenerator: old algo ", .pkgEnv$oldAlgoVTM)
    setorderv(pixelGroupData2, pgdAndSc)
    whNA <- unique(unlist(sapply(pixelGroupData2, function(xx) which(is.na(xx)))))
    pixelGroupData1 <- pixelGroupData1[!pixelGroupData2[whNA], on = pgdAndSc]
    setkeyv(pixelGroupData1, pgdAndSc)
    setkeyv(pixelGroupData2, pgdAndSc)
    aa <- pixelGroupData1[pixelGroupData2, on = pgdAndSc]
    if (!isTRUE(all.equal(aa[["speciesProportion"]], aa[["i.speciesProportion"]])))
      stop("Old algorithm in vegMapGenerator is different than new map")
  }

  if (algo == 1) {
    x <- cohortData1
    pixelGroupData <- pixelGroupData1
    rm(pixelGroupData1)
  } else if (algo == 2) {
    x <- cohortData2
    pixelGroupData <- pixelGroupData2
    rm(pixelGroupData2)
  }

  ########################################################
  #### Determine "mixed"
  ########################################################
  if (FALSE) {
    ## old algorithm; keep this code as reference -- it's simpler to follow
    b1 <- Sys.time()
    pixelGroupData4 <- x[, list(totalB = sum(B, na.rm = TRUE),
                                speciesCode, B), by = pixelGroup]
    pixelGroupData4 <- pixelGroupData4[, .(speciesGroupB = sum(B, na.rm = TRUE),
                                           totalB = totalB[1]),
                                       by = pgdAndSc]
    set(pixelGroupData4, NULL, "speciesProportion", pixelGroupData4$speciesGroupB /
          pixelGroupData4$totalB)
    pixelGroupData4[, speciesProportion := speciesGroupB / totalB]
    b2 <- Sys.time()
    mussage(b2 - b1)
    all.equal(pixelGroupData4[, .(pixelGroup, speciesCode, totalB)], pixelGroupData[
      , .(pixelGroup, speciesCode, totalB)])
  }

  if (mixedType == 1) {
    ## create "mixed" class #    -- Eliot May 28, 2019 -- faster than previous below
    ## 1. anything with >= vegLeadingProportion is "pure"
    ## 2. sort on pixelGroup and speciesProportion, reverse so that 1st row of each pixelGroup is the largest
    ## 3. Keep only first row in each pixelGroup
    ## 4. change column names and convert pure to mixed ==> mixed <- !pure
    pixelGroupData3 <- pixelGroupData[, list(pure = speciesProportion >= vegLeadingProportion,
                                             speciesCode, pixelGroup, speciesProportion)]
    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName]
    pixelGroupData3[pure == FALSE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    pixelGroupData3[, pure := !pure]
    setnames(pixelGroupData3, "pure", "mixed")

    ## Old algorithm for above, this is ~43 times slower
    # a2 <- Sys.time()
    # pixelGroupData2 <- pixelGroupData[, list(mixed = all(speciesProportion < vegLeadingProportion),
    #                                          leading = speciesCode[which.max(speciesProportion)]),
    #                                   by = pixelGroupColName]
    # pixelGroupData2[mixed == TRUE, leading := "Mixed"]
    # b2 <- Sys.time()
  } else if (mixedType == 2) {
    if (!sppEquivCol %in% colnames(sppEquiv))
      stop(sppEquivCol, " is not in sppEquiv. Please pass an existing sppEquivCol")

    sppEq <- data.table(sppEquiv[[sppEquivCol]], sppEquiv[["Type"]])

    names(sppEq) <- c("speciesCode", "Type")
    setkey(pixelGroupData, speciesCode)

    # don't need all columns now
    colsToDelete <- c("rasterToMatch", leadingBasedOn, totalOfLeadingBasedOn)
    colsToDelete <- colsToDelete[colsToDelete %in% colnames(pixelGroupData)]
    set(pixelGroupData, NULL, colsToDelete, NULL)
    pixelGroupData3 <- merge(pixelGroupData, sppEq[!duplicated(sppEq)], all.x = TRUE)
    setkeyv(pixelGroupData, pgdAndSc)

    setkeyv(pixelGroupData3, pgdAndSc)
    mixedType2Condition <- quote(Type == "Deciduous" &
                                   speciesProportion < vegLeadingProportion &
                                   speciesProportion > 1 - vegLeadingProportion)
    pixelGroupData3[, mixed := FALSE]

    if (algo == 2 || isTRUE(doAssertion)) {
      b <- pixelGroupData3[, list(N = .N), by = pixelGroupColName]
      b <- rep.int(b[["N"]], b[["N"]])
      GT1 <- b > 1


      pgd3GT1 <- pixelGroupData3[GT1]
      pgd3NGT1 <- pixelGroupData3[!GT1]

      pgd3GT1[eval(mixedType2Condition), mixed := TRUE, by = pixelGroupColName]
      pgd3GT1[, mixed := any(mixed), by = pixelGroupColName]
      pixelGroupData3 <- rbindlist(list(pgd3NGT1, pgd3GT1))
    } else {
      pixelGroupData3[eval(mixedType2Condition), mixed := TRUE, by = pixelGroupColName]
      pixelGroupData3[, mixed := any(mixed), by = pixelGroupColName]
    }

    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    set(pixelGroupData3, NULL, "Type", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName] ## sp. w/ highest prop. per pixelGroup
    pixelGroupData3[mixed == TRUE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    set(pixelGroupData3, NULL, "leading", factor(pixelGroupData3[["leading"]]))
  } else {
    stop("invalid mixedType! Must be one of '1' or '2'.")
  }

  if ("pixelIndex" %in% colnames(pixelGroupData3)) { #(!any(duplicated(pixelGroupData3[[pixelGroupColName]]))) {
    if (!is.factor(pixelGroupData3[["leading"]])) {
      pixelGroupData3[["leading"]] <- factor(pixelGroupData3[["leading"]])
    }

    vegTypeMap <- raster(pixelGroupMap)

    vegTypeMap[pixelGroupData3[["pixelIndex"]]] <- pixelGroupData3[["leading"]]
    levels(vegTypeMap) <- data.frame(ID = seq_along(levels(pixelGroupData3[["leading"]])),
                                     species = levels(pixelGroupData3[["leading"]]),
                                     stringsAsFactors = TRUE)
  } else {
    if (is.factor(pixelGroupData3[["initialEcoregionCode"]])) {
      f <- pixelGroupData3[["initialEcoregionCode"]]
      pixelGroupData3[["initialEcoregionCode"]] <- as.numeric(levels(f))[f]
    }

    vegTypeMap <- rasterizeReduced(pixelGroupData3, pixelGroupMap, "leading", pixelGroupColName)
  }

  if (missing(colors)) {
    colors <- if (!"Mixed" %in% sppEquiv[[sppEquivCol]])
      sppColors(sppEquiv, sppEquivCol, newVals = "Mixed", palette = "Accent")
    else
      sppColors(sppEquiv, sppEquivCol, palette = "Accent")
  }

  levels(vegTypeMap) <- cbind(levels(vegTypeMap)[[1]],
                              colors = colors[match(levels(vegTypeMap)[[1]][[2]], names(colors))],
                              stringsAsFactors = FALSE)
  setColors(vegTypeMap, n = length(colors)) <- levels(vegTypeMap)[[1]][, "colors"]

  if (isTRUE(doAssertion)) {
    if (sum(!is.na(vegTypeMap[])) < 100)
      ids <- sample(which(!is.na(vegTypeMap[])), sum(!is.na(vegTypeMap[])), replace = FALSE)
    else
      ids <- sample(which(!is.na(vegTypeMap[])), 100, replace = FALSE)

    pgs <- pixelGroupMap[][ids]
    dups <- duplicated(pgs)
    if (any(dups)) {
      pgs <- pgs[!dups]
      ids <- ids[!dups]
    }
    leadingTest <- factorValues2(vegTypeMap, vegTypeMap[ids], att = 2)
    names(pgs) <- leadingTest
    pgTest <- pixelGroupData[get(pixelGroupColName) %in% pgs]
    whNA <- unique(unlist(sapply(pgTest, function(xx) which(is.na(xx)))))
    whNAPG <- pgTest[whNA][[pixelGroupColName]]
    pgs <- pgs[!pgs %in% whNAPG]
    pgTest <- pgTest[!whNA]
    pgTest[, Type := equivalentName(speciesCode, sppEquiv, "Type")]
    if (mixedType == 1) {
      pgTest2 <- pgTest[, list(mixed = all(speciesProportion < vegLeadingProportion),
                               leading = speciesCode[which.max(speciesProportion)]),
                        by = pixelGroupColName]
      out <- pgTest2[mixed == TRUE, leading := "Mixed"]
    } else if (mixedType == 2) {
      pgTest2 <- pgTest[, list(mixed = eval(mixedType2Condition),
                               leading = speciesCode[which.max(speciesProportion)]),
                        by = pixelGroupColName]
      pgTest2[, mixed := any(mixed), by = pixelGroupColName]
      pgTest2[mixed == TRUE, leading := "Mixed"]
      pgTest2 <- pgTest2[!duplicated(pgTest2)]
      out <- pgTest2
    } ## TODO: make sure mixedType = 0 works!! currently not implemented

    length(unique(out[[pixelGroupColName]]))
    length(pgs %in% unique(pgTest2[[pixelGroupColName]]))
    pgs2 <- pgs[pgs %in% unique(pgTest2[[pixelGroupColName]])]
    if (!isTRUE(all(setkeyv(pgTest2, pixelGroupColName)[["leading"]] == names(pgs)[order(pgs)])))
      stop("The vegTypeMap is incorrect. Please debug LandR::vegTypeMapGenerator")
  }

  return(vegTypeMap)
}

#' Load kNN species layers from online data repository
#'
#' TODO: description needed
#'
#' @param dPath path to the data directory
#'
#' @template rasterToMatch
#'
#' @template studyArea
#'
#' @template sppEquiv
#'
#' @param knnNamesCol character string indicating the column in \code{sppEquiv}
#'                    containing kNN species names.
#'                    Default \code{"KNN"} for when \code{sppEquivalencies_CA} is used.
#'
#' @template sppEquivCol
#'
#' @param thresh the minimum number of pixels where the species must have
#'               \code{biomass > 0} to be considered present in the study area.
#'               Defaults to 1.
#'
#' @param url the source url for the data, passed to \code{\link[reproducible]{prepInputs}}
#'
#' @param ... Additional arguments passed to \code{\link[reproducible]{Cache}}
#'            and \code{\link{equivalentName}}. Also valid: \code{outputPath}, and \code{studyAreaName}.
#'
#' @return A raster stack of percent cover layers by species.
#'
#' @export
#' @importFrom httr config with_config
#' @importFrom magrittr %>%
#' @importFrom raster ncell raster
#' @importFrom reproducible Cache .prefix preProcess basename2
#' @importFrom tools file_path_sans_ext
#' @importFrom utils capture.output untar
#' @importFrom RCurl getURL
#' @importFrom XML getHTMLLinks
loadkNNSpeciesLayers <- function(dPath, rasterToMatch, studyArea, sppEquiv,
                                 knnNamesCol = "KNN", sppEquivCol, thresh = 1, url, ...) {
  dots <- list(...)
  oPath <- if (!is.null(dots$outputPath)) dots$outputPath else dPath

  sppEquiv <- sppEquiv[, lapply(.SD, as.character)]
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  ## remove empty names
  sppNameVector <- sppNameVector[sppNameVector != ""]

  sppMerge <- unique(sppEquiv[[sppEquivCol]][duplicated(sppEquiv[[sppEquivCol]])])
  sppMerge <- sppMerge[nzchar(sppMerge)]
  if ("cachePath" %in% names(dots)) {
    cachePath <- dots$cachePath
  } else {
    cachePath <- getOption("reproducible.cachePath")
  }

  ## get all files in url folder
  fileURLs <- getURL(url, dirlistonly = TRUE,
                     .opts = list(followlocation = TRUE,
                                  ssl.verifypeer = 0L)) ## TODO: re-enable verify
  fileNames <- getHTMLLinks(fileURLs)
  ## get all kNN species - names only
  allSpp <- grep("2001_kNN_Species_.*\\.tif$", fileNames, value = TRUE) %>%
    sub("_v1.tif", "", .) %>%
    sub(".*Species_", "", .)

  if (getRversion() < "4.0.0") {
    if (length(allSpp) == 0)
      stop("Incomplete file list retrieved from server.")
  } else {
    stopifnot("Incomplete file list retrieved from server." = length(allSpp) > 1)
  }

  ## get all species layers from .tar
  if (length(sppNameVector) == 1) ## avoids a warning in next if
    if (sppNameVector == "all")
      sppNameVector <- allSpp

  ## Make sure spp names are compatible with kNN names
  kNNnames <- as.character(equivalentName(sppNameVector, sppEquiv, column = knnNamesCol,
                                          multi = TRUE))
  sppNameVector <- as.character(equivalentName(sppNameVector, sppEquiv, column = sppEquivCol,
                                               multi = TRUE))

  ## if there are NA's, that means some species can't be found in kNN database
  if (any(is.na(kNNnames))) {
    warning(paste0("Can't find ", sppNameVector[is.na(kNNnames)], " in `sppEquiv$",
                   knnNamesCol, ".\n",
                   "Will use remaining matching species, but check if this is correct."))
    ## select only available species
    sppNameVector <- sppNameVector[!is.na(kNNnames)]
    kNNnames <- kNNnames[!is.na(kNNnames)]
  }

  emptySppNames <- kNNnames == ""
  if (any(emptySppNames)) {
    ## select only available species
    kNNnames <- kNNnames[!emptySppNames]
    sppNameVector <- sppNameVector[!emptySppNames]
  }

  ## same as above
  if (any(!kNNnames %in% allSpp)) {
    warning(paste0("Can't find ", sppNameVector[!kNNnames %in% allSpp], " in kNN database.\n",
                   "Will use remaining matching species, but check if this is correct."))
    sppNameVector <- sppNameVector[kNNnames %in% allSpp]
    kNNnames <- kNNnames[kNNnames %in% allSpp]
  }

  ## define suffix to append to file names
  suffix <- if (basename(cachePath) == "cache") {
    paste0(as.character(ncell(rasterToMatch)), "px")
  } else {
    basename(cachePath)
  }
  suffix <- paste0("_", suffix)

  ## select which targetFiles to extract
  targetFiles <- grep(paste(kNNnames, ".*\\.tif$", sep = "", collapse = "|"), fileNames, value = TRUE)
  postProcessedFilenames <- .suffix(targetFiles, suffix = suffix)
  postProcessedFilenamesWithStudyAreaName <- .suffix(postProcessedFilenames, paste0("_", dots$studyAreaName))

  message("Running prepInputs for ", paste(kNNnames, collapse = ", "))
  if (length(kNNnames) > 15) {
    message("This looks like a lot of species;",
            " did you mean to pass only a subset of this to sppEquiv?\n",
            " You can use the list above to choose species, then select only those rows",
            " in sppEquiv before passing here.")
  }
  with_config(config = config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
    speciesLayers <- Cache(Map,
                           targetFile = asPath(targetFiles),
                           filename2 = postProcessedFilenamesWithStudyAreaName,
                           url = paste0(url, targetFiles),
                           MoreArgs = list(destinationPath = asPath(dPath),
                                           fun = "raster::raster",
                                           studyArea = studyArea,
                                           rasterToMatch = rasterToMatch,
                                           method = "bilinear",
                                           datatype = "INT2U",
                                           overwrite = TRUE,
                                           userTags = dots$userTags
                           ),
                           prepInputs, quick = TRUE) # don't need to digest all the "targetFile"
  })
  names(speciesLayers) <- unique(kNNnames) ## TODO: see #10

  # remove "no data" first
  noData <- sapply(speciesLayers, function(xx) is.na(maxValue(xx)))
  if (any(noData)) {
    message(names(noData)[noData], " has no data in this study area; omitting it")
    speciesLayers <- speciesLayers[!noData]
  }

  # remove "little data" first
  layersWdata <- sapply(speciesLayers, function(xx) if (maxValue(xx) < thresh) FALSE else TRUE)
  if (sum(!layersWdata) > 0) {
    sppKeep <- names(speciesLayers)[layersWdata]
    if (length(sppKeep)) {
      message("removing ", sum(!layersWdata), " species because they had <", thresh,
              " % cover in the study area\n",
              "  These species are retained (and could be further culled manually, if desired):\n",
              paste(sppKeep, collapse = " "))
    } else {
      message("no pixels for ", paste(names(layersWdata), collapse = " "),
              " were found with >=", thresh, " % cover in the study area.",
              "\n  No species layers were retained. Try lowering the threshold",
              " to retain species with low % cover")
    }
  }
  speciesLayers <- speciesLayers[layersWdata]
  if (!is.null(sppMerge)) {
    if (length(sppMerge) == 0) {
      lapply(1:length(speciesLayers), function(i, rasters = speciesLayers,
                                               filenames = postProcessedFilenamesWithStudyAreaName) {
        writeRaster(rasters[[i]], file.path(oPath, paste0(filenames[i], '.tif')), overwrite = TRUE)
      })
    } else {
      speciesLayers <- mergeSppRaster(sppMerge = sppMerge, speciesLayers = speciesLayers,
                                      sppEquiv = sppEquiv, column = "KNN", suffix = suffix,
                                      dPath = oPath)
    }
  }
  ## Rename species layers - There will be 2 groups -- one
  nameChanges <- equivalentName(names(speciesLayers), sppEquiv, column = sppEquivCol)
  nameChangeNA <- is.na(nameChanges)
  names(speciesLayers)[!nameChangeNA] <- nameChanges[!nameChangeNA]

  nameChangesNonMerged <- equivalentName(names(speciesLayers)[nameChangeNA],
                                         sppEquiv, column = sppEquivCol)
  names(speciesLayers)[nameChangeNA] <- nameChangesNonMerged

  ## return stack and updated species names vector
  if (length(speciesLayers)) {
    stack(speciesLayers)
  }
}

#' Load kNN species layers from online data repository
#'
#' Downloads the 2011 kNN species cover layers from the Canadian Forestry Service,
#' National Inventory System, for validation.
#'
#' @param dPath path to the data directory
#'
#' @template rasterToMatch
#'
#' @template studyArea
#'
#' @template sppEquiv
#'
#' @param knnNamesCol character string indicating the column in \code{sppEquiv}
#'                    containing kNN species names.
#'                    Default \code{"KNN"} for when \code{sppEquivalencies_CA} is used.
#'
#' @template sppEquivCol
#'
#' @param thresh the minimum number of pixels where the species must have
#'               \code{biomass > 0} to be considered present in the study area.
#'               Defaults to 1.
#'
#' @param url the source url for the data, passed to \code{\link[reproducible]{prepInputs}}
#'
#' @param ... Additional arguments passed to \code{\link[reproducible]{Cache}}
#'            and \code{\link{equivalentName}}.
#'
#' @return A raster stack of percent cover layers by species.
#'
#' @export
#' @importFrom httr config with_config
#' @importFrom magrittr %>%
#' @importFrom raster ncell raster
#' @importFrom RCurl getURL
#' @importFrom reproducible basename2 Cache .prefix preProcess
#' @importFrom tools file_path_sans_ext
#' @importFrom utils capture.output untar
#' @importFrom XML getHTMLLinks
loadkNNSpeciesLayersValidation <- function(dPath, rasterToMatch, studyArea, sppEquiv,
                                           knnNamesCol = "KNN", sppEquivCol, thresh = 1, url, ...) {
  dots <- list(...)
  oPath <- if (!is.null(dots$outputPath)) dots$outputPath else dPath

  sppEquiv <- sppEquiv[, lapply(.SD, as.character)]
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  ## remove empty names
  sppNameVector <- sppNameVector[sppNameVector != ""]

  sppMerge <- unique(sppEquiv[[sppEquivCol]][duplicated(sppEquiv[[sppEquivCol]])])
  sppMerge <- sppMerge[nzchar(sppMerge)]
  if ("cachePath" %in% names(dots)) {
    cachePath <- dots$cachePath
  } else {
    cachePath <- getOption("reproducible.cachePath")
  }

  ## get all online file names
  fileURLs <- getURL(url, dirlistonly = TRUE,
                     .opts = list(followlocation = TRUE,
                                  ssl.verifypeer = 0L)) ## TODO: re-enable verify
  fileNames <- getHTMLLinks(fileURLs)
  fileNames <- grep("Species_.*.tif$", fileNames, value = TRUE)

  ## get all kNN species - names only
  allSpp <- fileNames %>%
    sub("_v1.tif", "", .) %>%
    sub(".*Species_", "", .)

  ## get all species layers from .tar
  if (length(sppNameVector) == 1) ## avoids a warning in next if
    if (sppNameVector == "all")
      sppNameVector <- allSpp

  ## Make sure spp names are compatible with kNN names
  kNNnames <- as.character(equivalentName(sppNameVector, sppEquiv,
                                          column = knnNamesCol, multi = TRUE))
  sppNameVector <- as.character(equivalentName(sppNameVector, sppEquiv,
                                               column = sppEquivCol, multi = TRUE))

  ## if there are NA's, that means some species can't be found in kNN data base
  if (any(is.na(kNNnames))) {
    warning(paste0("Can't find ", sppNameVector[is.na(kNNnames)], " in `sppEquiv$",
                   knnNamesCol, ".\n",
                   "Will use remaining matching species, but check if this is correct."))
    ## select only available species
    sppNameVector <- sppNameVector[!is.na(kNNnames)]
    kNNnames <- kNNnames[!is.na(kNNnames)]
  }

  emptySppNames <- kNNnames == ""
  if (any(emptySppNames)) {
    ## select only available species
    kNNnames <- kNNnames[!emptySppNames]
    sppNameVector <- sppNameVector[!emptySppNames]
  }

  ## same as above
  if (any(!kNNnames %in% allSpp)) {
    warning(paste0("Can't find ", kNNnames[!kNNnames %in% allSpp], " in kNN database.\n",
                   "Will use remaining matching species, but check if this is correct."))
    sppNameVector <- sppNameVector[kNNnames %in% allSpp]
    kNNnames <- kNNnames[kNNnames %in% allSpp]
  }

  ## define suffix to append to file names
  suffix <- if (basename(cachePath) == "cache") {
    paste0(as.character(ncell(rasterToMatch)), "px")
  } else {
    basename(cachePath)
  }
  suffix <- paste0("_", suffix)

  ## select which archives/targetFiles to extract -- because there was "multi" above, need unique here
  targetFiles <- fileNames[allSpp %in% kNNnames]
  postProcessedFilenames <- .suffix(targetFiles, suffix = suffix)
  postProcessedFilenamesWithStudyAreaName <- .suffix(postProcessedFilenames, paste0("_", dots$studyAreaName))
  message("Running prepInputs for ", paste(kNNnames, collapse = ", "))
  if (length(kNNnames) > 15) {
    message("This looks like a lot of species; did you mean to pass only a subset of this to sppEquiv?",
            "\n  You can use the list above to choose species, then select only those rows ",
            "\n  in sppEquiv before passing here")
  }

  speciesLayers <- list()

  with_config(config = config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
    speciesLayers <- Cache(Map,
                           targetFile = asPath(targetFiles),
                           filename2 = postProcessedFilenamesWithStudyAreaName,
                           url = paste0(url, targetFiles),
                           MoreArgs = list(destinationPath = asPath(dPath),
                                           fun = "raster::raster",
                                           studyArea = studyArea,
                                           rasterToMatch = rasterToMatch,
                                           method = "bilinear",
                                           datatype = "INT2U",
                                           overwrite = TRUE,
                                           userTags = dots$userTags
                           ),
                           prepInputs, quick = TRUE) # don't need to digest all the "targetFile"
  })

  names(speciesLayers) <- unique(kNNnames) ## TODO: see #10

  # remove "no data" first
  noData <- sapply(speciesLayers, function(xx) is.na(maxValue(xx)))
  if (any(noData)) {
    message(names(noData)[noData], " has no data in this study area; omitting it")
    speciesLayers <- speciesLayers[!noData]
  }

  layersWdata <- sapply(speciesLayers, function(xx) if (maxValue(xx) < thresh) FALSE else TRUE)
  if (sum(!layersWdata) > 0) {
    sppKeep <- names(speciesLayers)[layersWdata]
    if (length(sppKeep)) {
      message("removing ", sum(!layersWdata), " species because they had <",thresh,
              " % cover in the study area.\n",
              " These species are retained (and could be further culled manually, if desired):\n  ",
              paste(sppKeep, collapse = " "))
    } else {
      message("no pixels for ", paste(names(layersWdata), collapse = " "),
              " were found with >=", thresh, " % cover in the study area.\n",
              "No species layers were retained. Try lowering the threshold",
              " to retain species with low % cover")
    }
  }
  speciesLayers <- speciesLayers[layersWdata]
  if (!is.null(sppMerge)) {
    if (length(sppMerge) == 0) {
      lapply(1:length(speciesLayers), function(i, rasters = speciesLayers,
                                               filenames = postProcessedFilenamesWithStudyAreaName) {
        writeRaster(rasters[[i]], file.path(oPath, paste0(filenames[i], '.tif')), overwrite = TRUE)
      })
    } else {
      speciesLayers <- mergeSppRaster(sppMerge = sppMerge, speciesLayers = speciesLayers,
                                      sppEquiv = sppEquiv, column = "KNN", suffix = suffix,
                                      dPath = oPath)
    }
  }
  ## Rename species layers - There will be 2 groups -- one
  nameChanges <- equivalentName(names(speciesLayers), sppEquiv, column = sppEquivCol)
  nameChangeNA <- is.na(nameChanges)
  names(speciesLayers)[!nameChangeNA] <- nameChanges[!nameChangeNA]

  nameChangesNonMerged <- equivalentName(names(speciesLayers)[nameChangeNA],
                                         sppEquiv, column = sppEquivCol)
  names(speciesLayers)[nameChangeNA] <- nameChangesNonMerged

  ## return stack and updated species names vector
  if (length(speciesLayers)) {
    stack(speciesLayers)
  }
}

#' Function to sum rasters of species layers
#'
#' @param speciesLayers stack of species layers rasters
#' @param layersToSum names/indices of layers to be summed - optional
#' @param filenameToSave file path to save output raster
#' @param newLayerName name of the output raster layer
#'
#' @export
#' @importFrom raster calc stack writeRaster
sumRastersBySpecies <- function(speciesLayers, layersToSum, filenameToSave, newLayerName) {
  out <- raster::calc(raster::stack(speciesLayers[layersToSum]), sum)
  names(out) <- newLayerName
  writeRaster(out, filename = filenameToSave, datatype = "INT2U", overwrite = TRUE)
  out # Work around for Cache
}

#' Overlay layers within raster stacks
#'
#' Overlays rasters of different data resolution by filling in gaps in the highest
#' resolution raster with data available in lowest resolution one.
#' If only high or low resolution data are available, it will use it without
#' attempting to overlay.
#'
#' @param highQualityStack      high quality list/stack of rasters
#'                              (will be used preferentially)
#' @param lowQualityStack       low quality list/stack of rasters
#'                              (will be used to fill \code{NA}s in \code{highQualityStack})
#' @param outputFilenameSuffix  file suffix to save raster if there was overlaying.
#'                              Defaults to \code{"overlay"}.
#' @param destinationPath       directory for saved rasters
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom quickPlot layerNames
#' @importFrom raster ncell res stack
overlayStacks <- function(highQualityStack, lowQualityStack, outputFilenameSuffix = "overlay",
                          destinationPath) {
  ## check if there are any layers/values in the lowQualityStack
  ## if not return the HQ one
  if (class(lowQualityStack) != "RasterStack" &
      all(is.na(getValues(lowQualityStack)))) {
    highQualityStack
  } else {
    ## check if HQ resolution > LQ resolutions
    hqLarger <- ncell(lowQualityStack) * prod(res(lowQualityStack)) <
      ncell(highQualityStack) * prod(res(highQualityStack))

    ## make table of species layers in HQ and LQ
    dt1 <- data.table(SPP = layerNames(highQualityStack), HQ = layerNames(highQualityStack))
    dt2 <- data.table(SPP = layerNames(lowQualityStack), LQ = layerNames(lowQualityStack))
    setkey(dt1, SPP); setkey(dt2, SPP)
    dtj <- merge(dt1, dt2, all = TRUE)
    dtj[, c("HQ", "LQ") := list(!is.na(HQ), !is.na(LQ))]

    ## check which layers have species info in HQ and LQ
    #dtj[, HQ := any(!is.na(highQualityStack[[SPP]][])), by = 1:nrow(dtj)] #nolint
    #dtj[, LQ := any(!is.na(lowQualityStack[[SPP]][])), by = 1:nrow(dtj)] #nolint

    stackRas <- list()
    for (x in seq(nrow(dtj))) {
      stackRas[[x]] <- dtj[x, .overlay(SPP, HQ, LQ, hqLarger = hqLarger,
                                       highQualityStack = highQualityStack,
                                       lowQualityStack = lowQualityStack,
                                       outputFilenameSuffix = outputFilenameSuffix,
                                       destinationPath = destinationPath)]
    }
    names(stackRas) <- dtj$SPP

    stack(stackRas)
  }
}

#' Overlaying function
#'
#' Used internally in \code{overlayStacks}. Function to be applied to each row
#' of a \code{data.table} containing information of whether the species layer
#' exists in the HQ and LQ data.
#' Only overlays if data exists in both layers, otherwise returns the layer with data.
#'
#' @inheritParams overlayStacks
#' @param SPP \code{data.table} column of species layer name
#' @param HQ \code{data.table} column of whether \code{SPP} is present in HQ layers
#' @param LQ \code{data.table} column of whether \code{SPP} is present in LQ layers
#'
#' @importFrom gdalUtils gdalwarp
#' @importFrom raster compareRaster crs extent filename res projectExtent raster
#' @importFrom raster writeRaster xmax xmin ymax ymin
#' @keywords internal
.overlay <- function(SPP, HQ, LQ, hqLarger, highQualityStack, lowQualityStack, #nolint
                     outputFilenameSuffix = "overlay", destinationPath) {
  ## if HQ & LQ have data, pool
  if (HQ & LQ) {
    ## check equality of raster attributes and correct if necessary
    if (!all(
      isTRUE(all.equal(extent(lowQualityStack), extent(highQualityStack))),
      isTRUE(all.equal(crs(lowQualityStack), crs(highQualityStack))),
      isTRUE(all.equal(res(lowQualityStack), res(highQualityStack))))) {
      message("  ", SPP, " extents, or resolution, or projection did not match; ",
              "using gdalwarp to make them overlap")
      LQRastName <- basename(tempfile(fileext = ".tif"))
      if (!nzchar(filename(lowQualityStack[[SPP]]))) {
        LQCurName <- basename(tempfile(fileext = ".tif"))
        lowQualityStack[[SPP]][] <- as.integer(lowQualityStack[[SPP]][])
        lowQualityStack[[SPP]] <- writeRaster(lowQualityStack[[SPP]],
                                              filename = LQCurName,
                                              datatype = "INT2U")
      }

      LQRastInHQcrs <- projectExtent(lowQualityStack, crs = crs(highQualityStack))
      # project LQ raster into HQ dimensions
      gdalwarp(overwrite = TRUE,
               dstalpha = TRUE,
               s_srs = as.character(crs(lowQualityStack[[SPP]])),
               t_srs = as.character(crs(highQualityStack[[SPP]])),
               multi = TRUE, of = "GTiff",
               tr = res(highQualityStack),
               te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                      xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
               filename(lowQualityStack[[SPP]]), ot = "Byte",
               LQRastName)

      LQRast <- raster(LQRastName)
      LQRast[] <- LQRast[]
      unlink(LQRastName)

      try(unlink(LQCurName), silent = TRUE)

      if (hqLarger) {
        tmpHQName <- basename(tempfile(fileext = ".tif"))

        gdalwarp(overwrite = TRUE,
                 dstalpha = TRUE,
                 s_srs = as.character(crs(highQualityStack[[SPP]])),
                 t_srs = as.character(crs(highQualityStack[[SPP]])),
                 multi = TRUE, of = "GTiff",
                 tr = res(highQualityStack),
                 te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                        xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
                 filename(highQualityStack[[SPP]]), ot = "Byte", tmpHQName)
        HQRast <- raster(tmpHQName)
        HQRast[] <- HQRast[]
        HQRast[HQRast[] == 255] <- NA_integer_
        unlink(tmpHQName)
      } else {
        HQRast <- highQualityStack[[SPP]]
      }
    } else {
      LQRast <- lowQualityStack[[SPP]]
      HQRast <- highQualityStack[[SPP]]
    }

    message("  Writing new, overlaid ", SPP, " raster to disk.")
    if (!compareRaster(LQRast, HQRast))
      stop("Stacks not identical, something is wrong with overlayStacks function.")

    NAs <- is.na(HQRast[])

    ## complete missing HQ data with LQ data
    HQRast[NAs] <- LQRast[][NAs]
    HQRast <- writeRaster(HQRast, datatype = "INT1U",
                          filename = file.path(destinationPath,
                                               paste0(SPP, "_", outputFilenameSuffix, ".tif")),
                          overwrite = TRUE)
    names(HQRast) <- SPP
    return(HQRast)
  } else {
    ## if only HQ/LQ exist return one of them
    ## if none have data return one of the empty to keep all layers
    if (HQ) {
      HQRast <- highQualityStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    } else if (LQ) {
      LQRast <- lowQualityStack[[SPP]]
      names(LQRast) <- SPP
      return(LQRast)
    } else {
      HQRast <- highQualityStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    }
  }
}

#' Merge species percent-cover rasters
#'
#' Used internally in \code{overlayStacks}.
#'
#' @param sppMerge TODO
#' @param speciesLayers stack of species layers rasters
#' @template sppEquiv
#' @param column TODO
#' @param dPath destination path TODO
#' @param suffix TODO
#' @param ... Additional arguments TODO
#'
#' @importFrom raster calc stack writeRaster
#' @importFrom reproducible .suffix prepInputs
#' @keywords internal
mergeSppRaster <- function(sppMerge, speciesLayers, sppEquiv, column, suffix, dPath, ...) {
  stopifnot(!is.null(sppMerge), length(sppMerge) > 0)

  ## make sure species names and list names are in the right formats
  names(sppMerge) <- sppMerge
  sppMerges <- lapply(sppMerge, FUN = function(x) {
    unique(equivalentName(x, sppEquiv,  column = column, multi = TRUE))
  })

  ## keep species present in the data
  sppMerges <- lapply(sppMerges, FUN = function(x) x[x %in% names(speciesLayers)])

  for (i in seq(length(sppMerges))) {
    sumSpecies <- sppMerges[[i]]
    if (length(sumSpecies) > 1) {
      newLayerName <- names(sppMerges)[i]

      fname <- .suffix(file.path(dPath, paste0("kNN", newLayerName, ".tif")), suffix)
      a <- calc(stack(speciesLayers[sumSpecies]), sum, na.rm = TRUE)
      names(a) <- newLayerName
      a <- writeRaster(a, filename = fname, overwrite = TRUE, ...)
      ## replace spp rasters by the summed one
      speciesLayers[sumSpecies] <- NULL
      speciesLayers[[newLayerName]] <- a
      message("  Merging ", paste(sumSpecies, collapse = ", "), "; becoming: ", newLayerName)
    }
  }

  return(speciesLayers)
}

#' Rasterize polygons using \code{fasterize}
#'
#' @param sp a shapefile to rasterize
#' @param raster the template raster to use
#' @param fieldName the field to use (will be ignored if the shapefile has no fields)
#'
#' @return TODO: is it a \code{RasterLayer}?
#'
#' @export
#' @importFrom fasterize fasterize
#' @importFrom sf st_as_sf
fasterizeFromSp <- function(sp, raster, fieldName) {
  ## check if projections are the same
  if (!identical(crs(sp), crs(raster)))
    stop("fasterize will probably be wrong, as shp and raster projections do not match")

  tempSf <- sf::st_as_sf(sp)

  if (all(names(tempSf) == "geometry")) {
    ## if there are no fields, ignore fieldname
    fasterize::fasterize(tempSf, raster)
  } else
    fasterize::fasterize(tempSf, raster, field = fieldName)
}
