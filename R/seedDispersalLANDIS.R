utils::globalVariables(c(
  "..colnamesST", "RcvCommunity", "speciesCode2", "..speciesCodeCols"
))

#' Simulate a LANDIS-II dispersal process on a landscape.
#'
#' Simulate seed dispersal using user defined function. This is a "receiving pixel" focused
#' dispersal approach.
#' It is the "potentially receiving" cell that looks around itself for potential seed sources.
#' If it finds a single seed source, that passes the probability function described by the
#' `dispersalFn`.
#' If this passes a comparison to a uniform random draw, then the receiving cell is deemed to have
#' a "successful" dispersal for that species.
#' This function can therefore only be used for a relatively specific situation
#' where there is a yes/no returned for each potential receiving cell, i.e., not abundance.
#' This function is also not cumulative, i.e,. there is no higher abundance of seeds received if
#' a receiving cell has lots of seed sources around it vs. a single seed source.
#' The difference will come with a higher probability of successfully receiving a "seed".
#'
#' `dispersalFn` (temporarily unused as code is converted to Rcpp -- the
#' default `dispersalFn` is hard coded within the `spiralSeedDispersal`
#' function that uses C++) must be an expression that returns a probability
#' distribution. Because it is a dispersal kernel, it must be a probability
#' distribution. The expression that can take an argument named "dis" (without
#' quotes) as this will be calculated internally and represents the distance
#' from the initial (receiving) pixel and all active pixels within that cluster
#' of active pixels. `SpaDES` includes the [Ward()] kernel as
#' defined in the LANDIS-II documentation.
#'
#' @param dtSrc data.table
#'
#' @param dtRcv data.table
#'
#' @template pixelGroupMap
#'
#' @template speciesTable
#'
#' @param dispersalFn  An expression that can take a "dis" argument. See details.
#'   Default is "Ward" (temporarily unused, as it is hard coded inside Rcpp function)
#'
#' @param plot.it  Deprecated. If TRUE, then plot the raster at every interaction,
#'                 so one can watch the `LANDISDisp` event grow.
#' @param b  LANDIS Ward seed dispersal calibration coefficient (set to 0.01 in LANDIS)
#'
#' @param k  LANDIS Ward seed dispersal the probability that seed will disperse within
#' the effective distance (e.g., 0.95)
#'
#' @template successionTimestep
#'
#' @template verbose
#'
#' @param ...   Additional parameters. Currently none.
#'
#' @return A numeric vector of raster pixel indices, in the same resolution and extent as
#' `seedSrc` raster.
#'
#' @author Eliot McIntire
#' @export
#' @name LANDISDisp
#' @rdname LANDISDisp
#'
#' @examples
#' if (require("googledrive")) {
#'   seed <- sample(1e6, 1)
#'   set.seed(seed)
#'   library(data.table)
#'
#'   # keep this here for interactive testing with a larger raster
#'   rasterTemplate <- LandR:::rasterRead(terra::ext(0, 2500, 0, 2500), res = 100)
#'
#'   # make a pixelGroupMap
#'   pgs <- 4 # make even just because of approach below requires even
#'   pixelGroupMap <- SpaDES.tools::randomPolygons(rasterTemplate, numTypes = pgs)
#'   pixelGroupMap[1:100] <- NA # emulate a mask at the start
#'
#'   # Make a receive pixels table -- need pixelGroup and species
#'   nSpecies <- 3
#'   maxNSpeciesPerPixel <- min(5, nSpecies)
#'   rcvSpByPG <- lapply(seq_len(pgs / 2), function(pg) {
#'     data.table(speciesCode = sample(nSpecies, size = sample(maxNSpeciesPerPixel, 1)))
#'   })
#'   seedReceive <- rbindlist(rcvSpByPG, idcol = "pixelGroup")
#'
#'   # Make a source pixels table -- need pixelGroup and species
#'   srcSpByPG <- lapply(seq_len(pgs / 2), function(pg) {
#'     data.table(speciesCode = sample(nSpecies, size = sample(maxNSpeciesPerPixel, 1)))
#'   })
#'   seedSource <- rbindlist(srcSpByPG, idcol = "pixelGroup")
#'   # make source pixels not same pixelGroups as receive
#'   seedSource[, pixelGroup := pixelGroup + pgs / 2]
#'
#'   # Get a species table -- if using in Canada, can use this
#'   speciesTable <- getSpeciesTable(dPath = tempdir())
#'   speciesTable <- speciesTable[Area == "BSW"]
#'   speciesTable[, speciesCode := as.factor(LandisCode)]
#'   speciesTable[, seeddistance_eff := SeedEffDist]
#'   speciesTable[, seeddistance_max := SeedMaxDist]
#'
#'   speciesTable <- speciesTable
#'   speciesTable <- data.table(speciesTable)[, speciesCode := seq_along(LandisCode)]
#'   seedReceiveFull <- speciesTable[seedReceive, on = "speciesCode"]
#'   output <- LANDISDisp(
#'     dtRcv = seedReceiveFull, plot.it = interactive(),
#'     dtSrc = seedSource,
#'     speciesTable = speciesTable,
#'     pixelGroupMap,
#'     verbose = TRUE,
#'     successionTimestep = 10
#'   )
#'   # Summarize
#'   output[, .N, by = speciesCode]
#'
#'   ## Plot the maps
#'   if (interactive()) {
#'     library(quickPlot)
#'     clearPlot()
#'     spMap <- list()
#'     spMap$pixelGroupMap <- pixelGroupMap
#'     for (sppp in unique(output$speciesCode)) {
#'       spppChar <- paste0("Sp_", sppp)
#'       spMap[[spppChar]] <- LandR:::rasterRead(pixelGroupMap)
#'       ss <- unique(seedSource[speciesCode == sppp], on = c("pixelGroup", "speciesCode"))
#'       spMap[[spppChar]][pixelGroupMap[] %in% ss$pixelGroup] <- 1
#'
#'       receivable <- LandR:::rasterRead(pixelGroupMap)
#'       srf <- unique(seedReceiveFull[speciesCode == sppp], on = c("pixelGroup", "speciesCode"))
#'       receivable[pixelGroupMap[] %in% srf$pixelGroup] <- 1
#'
#'       forest <- which(!is.na(pixelGroupMap[]))
#'       src <- which(!is.na(spMap[[spppChar]][]))
#'       recvable <- which(!is.na(receivable[]))
#'       rcvd <- output[speciesCode == sppp]$pixelIndex
#'
#'       spMap[[spppChar]][forest] <- 0
#'       spMap[[spppChar]][recvable] <- 2
#'       spMap[[spppChar]][src] <- 1
#'       spMap[[spppChar]][rcvd] <- 3
#'       spMap[[spppChar]][intersect(src, rcvd)] <- 4
#'
#'       levels(spMap[[spppChar]]) <- data.frame(ID = 0:4,
#'                                               type = c("OtherForest", "Source", "Didn't receive",
#'                                                        "Received", "Src&Rcvd"))
#'     }
#'     Plot(spMap, cols = "Set2")
#'
#'     # A summary
#'     rr <- apply(rast(spMap)[[-1]][] + 1, 2, tabulate)
#'     rownames(rr) <- raster::levels(spMap[[2]])[[1]][,"type"][1:NROW(rr)]
#'     # next line only works if there are some places that are both source and potential to receive
#'     # rr <- rbind(rr, propSrcRcved = round(rr[5,]/ (rr[5,]+rr[2,]), 2))
#'  }
#' }
#'
LANDISDisp <- function(dtSrc, dtRcv, pixelGroupMap, speciesTable,
                       dispersalFn = Ward, b = 0.01, k = 0.95, plot.it = FALSE,
                       successionTimestep,
                       verbose = getOption("LandR.verbose", TRUE),
                       ...) {

  if ((NROW(dtSrc) > 0) & (NROW(dtRcv) > 0)) {
    ####### Assertions #############
    if (!( (is.numeric(dtSrc$speciesCode) && is.numeric(dtRcv$speciesCode) && is.numeric(speciesTable$speciesCode)) ||
           (is.factor(dtSrc$speciesCode) && is.factor(dtRcv$speciesCode) && is.factor(speciesTable$speciesCode))))
      stop("In LANDISDisp, dtSrc and dtRcv and speciesTable must each have columns for speciesCode which ",
           "must be both integer or both factor; they are not. Please correct this.")

    dtSrc <- data.table::copy(dtSrc)
    dtRcv <- data.table::copy(dtRcv)
    speciesTable <- data.table::copy(speciesTable)

    origClassWasNumeric <- is.numeric(speciesTable[["speciesCode"]])
    if (is(dtSrc$speciesCode, "numeric")) {
      #if (length(unique(dtSrc$speciesCode)) != max(dtSrc$speciesCode)) {
      # This is "numerics" that are no contiguous from 1
      set(dtSrc, NULL, c("speciesCode"), factor(dtSrc[["speciesCode"]]))
      origLevels <- levels(dtSrc[["speciesCode"]])
      speciesTable <- speciesTable[as.numeric(origLevels), ]
      set(speciesTable, NULL, c("speciesCode"),
          factor(speciesTable[["speciesCode"]], levels = origLevels))
      set(dtRcv, NULL, c("speciesCode"),
          factor(dtRcv[["speciesCode"]], levels = origLevels))
      #}
    }

    if (is.factor(dtSrc$speciesCode)) {
      if (!identical(levels(dtSrc$speciesCode), levels(dtRcv$speciesCode)) &&
          identical(levels(dtSrc$speciesCode), levels(speciesTable$speciesCode)))
        stop("In LANDISDisp, dtSrc$speciesCode and dtRcv$speciesCode and speciesTable$speciesCode ",
             "are all factors (good), ",
             "but they have different levels (bad). They must have the same factor levels.")
      origLevels <- levels(dtSrc$speciesCode)
      dtSrc[, speciesCode2 := as.integer(speciesCode)]
      dtRcv[, speciesCode2 := as.integer(speciesCode)]
      speciesTable[, speciesCode2 := as.integer(speciesCode)]
      set(dtSrc, NULL, "speciesCode", NULL)
      set(dtRcv, NULL, "speciesCode", NULL)
      set(speciesTable, NULL, "speciesCode", NULL)
      setnames(dtSrc, "speciesCode2", "speciesCode")
      setnames(dtRcv, "speciesCode2", "speciesCode")
      setnames(speciesTable, "speciesCode2", "speciesCode")
      if (!"species" %in% colnames(speciesTable))
        set(speciesTable, NULL, "species", paste0("Spp_", speciesTable[["speciesCode"]]))
      setorderv(speciesTable, "speciesCode")
      setorderv(dtSrc, "speciesCode")
      setorderv(dtRcv, "speciesCode")
    }

    if (is.character(dtSrc$speciesCode)) {
      stop("LANDISDisp expects that dtSrc$speciesCode and dtRcv$speciesCode are either factor or integer")
    }
    ####### End Assertions #############

    # Create srcPixelMatrix -- which is the matrix representation
    #    of the dtSrc x speciesCode with NAs everywhere there is no
    #    species present
    pgv <- pixelGroupMap[]

    rasVectorTemplate <- rep(NA_integer_, ncell(pixelGroupMap))
    rasTemplate <- rasterRead(pixelGroupMap)
    srcSpeciesCodes <- sort(unique(dtSrc$speciesCode))
    names(srcSpeciesCodes) <- as.character(srcSpeciesCodes)
    cellsCanSrc <- which(pgv %in% dtSrc$pixelGroup)
    dtSrcLong <- data.table(pixelGroup = pgv[cellsCanSrc], pixelIndex = cellsCanSrc)
    dtSrcLong <- dtSrc[, c("pixelGroup", "speciesCode")][dtSrcLong, on = "pixelGroup", allow.cartesian = TRUE]
    set(dtSrcLong, NULL, "pixelGroup", NULL)
    setkeyv(dtSrcLong, "speciesCode") # sort ascending

    srcSpeciesByIndex <- split(dtSrcLong, by = "speciesCode")
    speciesSrcRasterVecList <- lapply(srcSpeciesCodes, function(sc) {
      rasVectorTemplate[srcSpeciesByIndex[[as.character(sc)]][["pixelIndex"]]] <- sc
      rasVectorTemplate
    })
    maxSpCode <- max(as.integer(srcSpeciesCodes))
    speciesSrcRasterVecList <- lapply(seq_len(maxSpCode), function(ind) {
      if (as.character(ind) %in% names(speciesSrcRasterVecList)) {
        speciesSrcRasterVecList[[as.character(ind)]]
      }
    })
    srcPixelMatrix <- do.call(cbind, speciesSrcRasterVecList)

    #### cellSize -- faster than
    cellSize <- unique(res(pixelGroupMap))
    if (length(cellSize) > 1) {
      ## check for equal cell sizes that "aren't" due to floating point error
      res <-
        vapply(cellSize, function(x) {
          isTRUE(all.equal(x, cellSize[1]))
        }, logical(1))
      if (all(res)) {
        cellSize <- cellSize[1]
      } else {
        stop("pixelGroupMap resolution must be same in x and y dimension")
      }
    }

    #  Remove any species in dtRcv that are not available in dtSrc
    dtRcvNew <- dtRcv[unique(dtSrc[, "speciesCode"], by = "speciesCode"), on = "speciesCode",
                      nomatch = NULL]
    cellsCanRcv <- which(pgv %in% dtRcvNew$pixelGroup)
    rcvSpeciesCodes <- sort(unique(dtRcvNew$speciesCode))
    dtRcvLong <- data.table(pixelGroup = pgv[cellsCanRcv], pixelIndex = cellsCanRcv)
    dtRcvSmall <- dtRcvNew[, c("pixelGroup", "speciesCode")]
    dtSrcUniqueSP <- unique(dtSrc[, "speciesCode"], by = "speciesCode")
    dtRcvSmall1 <- dtRcvSmall[dtSrcUniqueSP, on = "speciesCode", nomatch = NULL]
    dtRcvLong <- dtRcvLong[dtRcvSmall, on = "pixelGroup", allow.cartesian = TRUE, nomatch = NULL]
    setorderv(dtRcvLong, c("pixelIndex", "speciesCode"))
    if (NROW(dtRcvLong)) {
      # There can be a case where a pixelGroup exists on map, with a species that is in Rcv but not in Src
      if (anyNA(dtRcvLong[["pixelIndex"]]))
        dtRcvLong <- na.omit(dtRcvLong)

      if (verbose >= 3) {
        message("numRcvPixels: ", length(unique(dtRcvLong$pixelIndex)),
                "; numSrcPixels: ", max(apply(srcPixelMatrix, 2, function(x) sum(!is.na(x)))))
      }
      dtRcvLong <- spiralSeedDispersalR(speciesTable, pixelGroupMap, dtRcvLong,
                                        srcPixelMatrix, cellSize, k, b, successionTimestep,
                                        verbose, dispersalFn = dispersalFn)
      if (exists("origLevels", inherits = FALSE)) {
        dtRcvLong[, speciesCode := factor(origLevels[speciesCode], levels = origLevels)]
        if (origClassWasNumeric) {
          set(dtRcvLong, NULL, "speciesCode", as.integer(as.character(dtRcvLong[["speciesCode"]])))
        }
      }
    }
  } else {
    dtRcvLong <- data.table(pixelIndex = integer(0), speciesCode = integer(0),
                            DistOfSuccess = numeric(0), ReasonForStop = character(0),
                            species = character(0))
  }
  return(dtRcvLong[])
}

speciesCodeFromCommunity <- function(num) {
  indices <- lapply(strsplit(intToBin2(num), split = ""), function(x) {
    rev(as.logical(as.numeric(x)))
  })

  speciesCode <- lapply(indices, function(x) (seq_len(length(x)) - 1)[x])
}

speciesComm <- function(num, sc) {
  speciesCode <- speciesCodeFromCommunity(num)
  a <- data.table(
    RcvCommunity = as.integer(rep(num, sapply(speciesCode, length))),
    speciesCode = unlist(speciesCode),
    key = "speciesCode"
  )[!is.na(speciesCode)]
  sc[a]
}

#' Ward Dispersal Kernel -- vectorized, optimized for speed
#'
#' A probability distribution used in LANDIS-II.
#'
#' @export
#' @docType methods
#'
#' @author Eliot McIntire
#' @param dist A vector of distances to evaluate kernel against
#' @param cellSize A numeric, length 1, of the cell resolution (e.g., res(raster))
#' @param effDist A vector of effective distance (parameter in kernel),
#'   with same length as `dist`
#' @param maxDist A vector of maximum distance (parameter in kernel),
#'   with same length as `dist`
#' @param k A parameter in the kernel
#' @param b A parameter in the kernel
#' @param algo Either 1 or 2. 2 is faster and is default. 1 is "simpler code" as it
#'   uses `ifelse`
#'
#' @name WardKernel
#' @rdname WardKernel
Ward <- function(dist, cellSize, effDist, maxDist, k, b, algo = 2) {
  if (length(maxDist) == 1) {
    if (length(dist) != 1)
      maxDist <- rep(maxDist, length(dist))
    else if (length(effDist) != 1)
      maxDist <- rep(maxDist, length(effDist))
  }
  if (length(dist) == 1) dist <- rep(dist, length(maxDist))
  if (length(effDist) == 1) effDist <- rep(effDist, length(maxDist))

  if (algo == 2) {
    # This is a 4 part ifelse statement
    # This is 3x faster, but less clear algorithm that uses 4 explicit logicals
    #   to create 4 indices, then evaluate each of the 4 expressions that go
    #   with the logicals
    out <- numeric(length(maxDist))

    # here are the expressions
    #  where a1, a2, b1, b2 are the 4 parts
    expra1 <- expression(
      exp((dista1wh - cellSize) * log(1 - k) / effDista1wh) -
        exp(dista1wh * log(1 - k) / effDista1wh))
    expra2 <- expression(
      (1 - k) * exp((dista2wh - cellSize - effDista2wh) * log(b) / maxDista2wh) -
        (1 - k) * exp((dista2wh - effDista2wh) * log(b) / maxDista2wh))
    exprb1 <- expression(
      exp((distb1wh - cellSize) * log(1 - k) / effDistb1wh) - (1 - k) *
        exp((distb1wh - effDistb1wh) * log(b) / maxDistb1wh))
    exprb2 <- expression(
      (1 - k) * exp((distb2wh - cellSize - effDistb2wh) * log(b) / maxDistb2wh) -
        (1 - k) * exp((distb2wh - effDistb2wh) * log(b) / maxDistb2wh))

    # Determine the indices for each one
    # Here are the 3 logicals, toplevel, then 2 nested
    a <- cellSize <= effDist
    awh <- which(a)
    bwh <- which(!a)
    a1 <- dist[awh] <= effDist[awh]
    b1 <- dist[bwh] <= cellSize

    a1wh <- awh[a1]
    a2wh <- awh[!a1]

    b1wh <- bwh[b1]
    b2wh <- bwh[!b1]

    # Build objects that are the dist, effDist, maxDist for each of the 4 subsets
    dista1wh <- dist[a1wh]
    effDista1wh <- effDist[a1wh]

    dista2wh <- dist[a2wh]
    effDista2wh <- effDist[a2wh]
    maxDista2wh <- maxDist[a2wh]

    distb1wh <- dist[b1wh]
    effDistb1wh <- effDist[b1wh]
    maxDistb1wh <- maxDist[b1wh]

    distb2wh <- dist[b2wh]
    effDistb2wh <- effDist[b2wh]
    maxDistb2wh <- maxDist[b2wh]

    # Create the content
    out[a1wh] <- eval(expra1)
    out[a2wh] <- eval(expra2)
    out[b1wh] <- eval(exprb1)
    out[b2wh] <- eval(exprb2)

    out
  } else {
    ifelse(cellSize <= effDist, {
      ifelse(
        dist <= effDist,
        exp((dist - cellSize) * log(1 - k) / effDist) -
          exp(dist * log(1 - k) / effDist),
        (1 - k) * exp((dist - cellSize - effDist) * log(b) / maxDist) -
          (1 - k) * exp((dist - effDist) * log(b) / maxDist)
      )
    }, {
      ifelse(
        dist <= cellSize,
        exp((dist - cellSize) * log(1 - k) / effDist) - (1 - k) *
          exp((dist - effDist) * log(b) / maxDist),
        (1 - k) * exp((dist - cellSize - effDist) * log(b) / maxDist) -
          (1 - k) * exp((dist - effDist) * log(b) / maxDist)
      )
    }
    )
  }
}

intToBin2 <- function(x) {
  y <- as.integer(x)
  class(y) <- "binmode"
  y <- as.character(y)
  dim(y) <- dim(x)
  y
}

spiralSeedDispersalR <- function(speciesTable, pixelGroupMap, dtRcvLong,
                                 srcPixelMatrix, cellSize, k, b,
                                 successionTimestep, verbose, dispersalFn) {
  # # quick test -- just first cell
  # if (missing(useMask))
  #   useMask <- isTRUE(is.na(srcPixelMatrix[1]) | is.na(all(tail(srcPixelMatrix, 1))))

  rcvFull <- dtRcvLong[, c("pixelIndex", "speciesCode")]
  rcvFull <- rcvFull[speciesTable[, c("seeddistance_max", "speciesCode")],
                     on = "speciesCode", nomatch = NULL]

  # # If there are entire rows that are NAs (quick test first -- just srcPixelMatrix --
  # #    then correct slower test inside this protected block)
  # #    then this section crops the pixelGroupMap and
  # #    adjusts the srcPixelMatrix accordingly
  # if (useMask) { # was an imprecise test, but quick test to get here; now reassess correct test
  #   useMask <- FALSE
  #   srcPixelMatrixNAs <- rowSums(srcPixelMatrix, na.rm = TRUE)
  #   whNAs <- which(srcPixelMatrixNAs > 0)
  #   numInitialNAs <- whNAs[1] - 1
  #   numFinalNAs <- tail(whNAs, 1)
  #   numColsPGM <- ncol(pixelGroupMap)
  #   newExtent <- extent(pixelGroupMap)
  #   rowsToRm <- floor(numInitialNAs / numColsPGM)
  #   ncells <- ncell(pixelGroupMap)
  #   rowsToRmEnd <- floor( (ncells - numFinalNAs) / numColsPGM)
  #   if (rowsToRm > 0) {
  #     cellsToRmInit <- rowsToRm * numColsPGM
  #     cellsToRmEnd <- rowsToRmEnd * numColsPGM
  #     cellsToRmInitInds <- seq(cellsToRmInit)
  #     cellsToRmEndInds <- (ncells - cellsToRmEnd + 1):ncells
  #     cellsToOmit <- c(cellsToRmInitInds, cellsToRmEndInds)
  #     whKeep <- sort(union(setdiff(seq(ncells), cellsToOmit),
  #                     unique(rcvFull$pixelIndex)))
  #     rc22 <- rowFromCell(pixelGroupMap, whKeep)
  #     rangeRC22 <- range(rc22)
  #     numRowsNew <- diff(rangeRC22) + 1
  #     if (numRowsNew < nrow(pixelGroupMap)) {
  #       newExtent@ymax <- newExtent@ymin + numRowsNew * res(pixelGroupMap)[1]
  #       srcPixelMatrix <- srcPixelMatrix[whKeep,]
  #       if (min(rangeRC22) > 1) {
  #         useMask <- TRUE
  #         pixelIndexAdj <- (min(rangeRC22) - 1) * ncol(pixelGroupMap)
  #         set(rcvFull, NULL, "pixelIndex", rcvFull$pixelIndex - pixelIndexAdj)
  #       }
  #       pixelGroupMap <- crop(pixelGroupMap, newExtent)
  #
  #     }
  #   }
  # }
  # print(paste0("number cells in pixelGroupMap: ", ncell(pixelGroupMap)))
  # print(paste0("number rows in srcPixelMatrix: ", NROW(srcPixelMatrix)))
  # print(paste0("number cols in srcPixelMatrix: ", NCOL(srcPixelMatrix)))

  speciesTable <- copy(speciesTable)
  set(speciesTable, NULL, "seeddistance_maxMinCellSize",
      pmax(cellSize, speciesTable[["seeddistance_max"]]))
  maxDis <- max(speciesTable[, "seeddistance_maxMinCellSize"])

  # Main matrix of distances that occur in a spiral
  preExistingSpiral <- paste0("spirals_max", round(maxDis, 6), "_cell", round(cellSize, 6))
  if (!exists(preExistingSpiral, envir = .pkgEnv)) {
    .pkgEnv[[preExistingSpiral]] <- spiralDistances(pixelGroupMap, maxDis, cellSize)
  }
  spiral <- .pkgEnv[[preExistingSpiral]]

  speciesTableSmall <- speciesTable[, c("speciesCode", "seeddistance_eff", "seeddistance_max")]
  uniqueDists <- unique(spiral[, "dists", drop = FALSE]) * cellSize
  numSp <- NROW(speciesTable)
  spSeq <- seq(numSp)
  distsBySpCode <- as.data.table(expand.grid(dists = uniqueDists, speciesCode = speciesTable[["speciesCode"]]))
  set(distsBySpCode, NULL, "seeddistance_max", speciesTableSmall[distsBySpCode[["speciesCode"]],
                                                                 "seeddistance_max"])
  set(distsBySpCode, NULL, "seeddistance_eff", speciesTableSmall[distsBySpCode[["speciesCode"]],
                                                                 "seeddistance_eff"])
  set(distsBySpCode, NULL, "wardProb",
      pmin(1, dispersalFn(dist = distsBySpCode$dists, cellSize = cellSize, effDist = distsBySpCode$seeddistance_eff,
                          maxDist = distsBySpCode$seeddistance_max, k = k, b = b)))
  set(distsBySpCode, NULL, c("seeddistance_max", "seeddistance_eff"), NULL)
  setorderv(distsBySpCode, c("dists", "speciesCode"))

  rcvFull <- dtRcvLong[, c("pixelIndex", "speciesCode")]
  rcvFull <- rcvFull[speciesTable[, c("seeddistance_max", "speciesCode")],
                     on = "speciesCode", nomatch = NULL]
  activeFullIndex <- seq.int(NROW(rcvFull))

  # activeIndexShrinking <- activeFullIndex
  rc1 <- rowColFromCell(pixelGroupMap, rcvFull[["pixelIndex"]])
  colnames(rc1) <- c("row", "col")   ## terra::rowColFromCell output has no colnames
  rowOrig <- rc1[, "row"]
  colOrig <- rc1[, "col"]

  speciesCode <- rcvFull[["speciesCode"]]

  activeSpecies <- speciesTable[, c("seeddistance_max", "seeddistance_eff", "speciesCode", "seeddistance_maxMinCellSize")]
  overallMaxDist <- max(speciesTable[["seeddistance_max"]])
  cantDoShortcutYet <- TRUE
  lastWardMaxProb <- 1

  nrowSrcPixelMatrix <- NROW(srcPixelMatrix)
  # dim(srcPixelMatrix) <- NULL # make a single vector -- a bit faster

  ## Assertions for inputs to spiral
  if (!is.numeric(k)) stop("not numeric k")
  if (!is.numeric(b)) stop("not numeric b")
  if (!is.numeric(successionTimestep)) stop("not numeric successTimestep")

  prevCurDist <- spiral[1, 3]
  spiralRow <- spiral[, "row"]
  spiralCol <- spiral[, "col"]
  prevSpiralRow <- spiralRow[1]
  prevSpiralCol <- spiralCol[1]
  newCurDist <- TRUE
  newActiveIndex <- TRUE
  tooLong <- FALSE

  curDists <- drop(spiral[, 3]) * cellSize
  uniqueDistCounter <- 0

  startTime <- Sys.time()
  cumSuccesses <- 0

  for (i in 1:NROW(spiral)) {
    curDist <- curDists[i]

    if (i > 1) {
      if (curDist > prevCurDist) {
        newCurDist <- TRUE
        uniqueDistCounter <- uniqueDistCounter + 1
        prevCurDist <- curDist
        spTooLong <- curDist > activeSpecies$seeddistance_maxMinCellSize
        if (any(spTooLong)) {
          activeSpecies <- activeSpecies[!spTooLong]
          tooLong <- curDist > ( pmax(cellSize, rcvFull[["seeddistance_max"]][activeFullIndex]) ) # * sqrt(2)) don't need this because spiral is sorted by distance
          if (any(tooLong, na.rm = TRUE)) {
            if (verbose >= 1) {
              tooLongFull <- curDist > rcvFull[["seeddistance_max"]]
              set(rcvFull, which(tooLongFull & is.na(rcvFull$ReasonForStop)), "ReasonForStop", "NoneRecdBeforeMaxDistReached")
            }
            activeFullIndex <- activeFullIndex[!tooLong]
            rowOrig <- rowOrig[!tooLong]
            colOrig <- colOrig[!tooLong]
            speciesCode <- speciesCode[!tooLong]
            newActiveIndex <- TRUE
          }
        }
      } else {
        newCurDist <- FALSE
      }
    }

    needNewCol <- TRUE
    if (newActiveIndex) {
      rowNow <- rowOrig
      colNow <-  colOrig
    } else {
      if (identical(prevSpiralCol, spiralCol[i])) {
        needNewCol <- FALSE
      }
    }

    row <- if (spiralRow[i] == 0) rowNow else rowNow + spiralRow[i] #spiral[i, "row"] # faster?
    if (needNewCol) {
      col <- if (spiralCol[i] == 0) colNow else colNow + spiralCol[i]
    }

    prevSpiralRow <- spiralRow[i]
    prevSpiralCol <- spiralCol[i]

    # I tried to get this faster; but the problem is omitting all the
    #   row/col combinations that are "off" raster. cellFromRowCol does this
    #   internally and is fast
    newPixelIndex <- as.integer(cellFromRowCol(row = row, col = col, object = pixelGroupMap))

    # lookup on src rasters
    if (newActiveIndex) {
      activeSpeciesCode <- speciesCode#[activeFullIndex]
    }
    srcPixelMatrixInd <- (activeSpeciesCode - 1) * nrowSrcPixelMatrix + newPixelIndex
    srcPixelValues <- srcPixelMatrix[srcPixelMatrixInd]
    # browser()
    hasSpNot <- is.na(srcPixelValues)
    # hasSp <- is.na(srcPixelValues)
    if (newCurDist) {
      rowsForDistsBySpCode <- as.integer(uniqueDistCounter * numSp + spSeq)
      # whUse <- which(distsBySpCode$dists == curDist) # old slower
      wardProbActual <- distsBySpCode$wardProb[rowsForDistsBySpCode]

      # The calculations for dispersal kernal are annual --
      #   so exponentiate for any other time step
      if (successionTimestep > 1) {
        wardProbActual = 1 - (1 - wardProbActual) ^ successionTimestep
      }
    }

    sumHasSp <- length(hasSpNot) - sum(hasSpNot)
    if (sumHasSp) {
      ran <- runifC(sumHasSp)

      whRanLTprevMaxProb <- which(ran <= lastWardMaxProb)

      if (length(whRanLTprevMaxProb)) {
        whHasSp <- which(!hasSpNot)
        if (i == 1) { # self pixels are 100%
          oo <- seq.int(length(ran))
        } else {
          wardRes <- wardProbActual[activeSpeciesCode[whHasSp][whRanLTprevMaxProb]]
          lastWardMaxProb <- min(1, max(wardRes))
          oo <- ran[whRanLTprevMaxProb] < wardRes
          oo <- whRanLTprevMaxProb[oo]
        }
        numSuccesses <- length(oo)
        if (verbose >= 2) {
          print(paste0(i, "; curDist: ",round(curDist,0),"; NumSuccesses: ",
                       numSuccesses, "; NumRows: ", NROW(na.omit(activeFullIndex)),
                       "; NumSp: ", length(unique(speciesCode[activeFullIndex]))))
        }
        notActiveSubIndex <- whHasSp[oo]
        if (length(notActiveSubIndex)) {
          notActiveFullIndex <- activeFullIndex[notActiveSubIndex]

          # This next block does 1 of 2 things: either resize the vectors
          #    (rowOrig, colOrig, speciesCode, activeFullIndex), or just set the
          #    values that are not active to NA. Apparently, setting NA is quite a
          #    bit faster, up to a point. So, only resize objects every once in a while
          cumSuccesses <- cumSuccesses + numSuccesses
          if (cumSuccesses > 2000) {
            cumSuccesses <- 0
            # if (i %% modulo == 0) {
            elapsedTime <- Sys.time() - startTime
            ss <- seq(activeFullIndex)
            ss <- ss[-notActiveSubIndex]
            activeFullIndexNAs <- is.na(activeFullIndex)
            ss <- ss[!activeFullIndexNAs[-notActiveSubIndex]]
            activeFullIndex <- activeFullIndex[ss]
            rowOrig <- rowOrig[ss]
            colOrig <- colOrig[ss]
            speciesCode <- speciesCode[ss]
          } else {
            # Don't need to do colOrig or speciesCode as these are automatically NAs
            #   downstream because rowOrig is NA -- cuts off 10% of computation time
            activeFullIndex[notActiveSubIndex] <- NA
            rowOrig[notActiveSubIndex] <- NA
          }

          newActiveIndex <- TRUE

          set(rcvFull, notActiveFullIndex, "Success", TRUE)
          if (verbose >= 1) {
            set(rcvFull, notActiveFullIndex, "DistOfSuccess", curDist)
            set(rcvFull, notActiveFullIndex, "ReasonForStop", "SuccessFullSeedRcvd")
          }
        } else {
          notActiveSubIndex <- integer()
          newActiveIndex <- FALSE
        }

      } else {
        newActiveIndex <- FALSE
      }
    }

    if (length(activeFullIndex) == 0) {
      break
    }
  }

  set(rcvFull, NULL, c("seeddistance_max"), NULL)

  if (!"Success" %in% colnames(rcvFull)) {
    if (verbose >= 1) {
      set(rcvFull, NULL, "ReasonForStop", "NoSuccesses")
      ReasonForStop <- rcvFull
    }
    rcvFull <- rcvFull[0]
  } else {
    # Remove the unsuccessful ones
    whSuccess <- which(!is.na(rcvFull[["Success"]]))
    if (verbose >= 1) {
      fails <- which(is.na(rcvFull[["DistOfSuccess"]]))
      if (length(fails)) {
        set(rcvFull, fails, "ReasonForStop", "RanOutOfDistance")
      }
      whFails <- sort(c(fails, whSuccess))
      ReasonForStop <- rcvFull[whFails]
    }
    rcvFull <- rcvFull[whSuccess]
    set(rcvFull, NULL, c("Success"), NULL)
  }

  speciesCodeCols <- intersect(c("species", "speciesCode"), colnames(speciesTable))
  rcvFull <- rcvFull[speciesTable[, ..speciesCodeCols], on = "speciesCode",
                     nomatch = NULL]

  if (verbose >= 1) {
    setattr(rcvFull, "ReasonForStop", ReasonForStop)
  }
  # if (useMask) {
  #   # pixelIndexAdj <- (min(rangeRC22) - 1) * ncol(pixelGroupMap)
  #   set(rcvFull, NULL, "pixelIndex", rcvFull$pixelIndex + pixelIndexAdj)
  # }
  rcvFull
}

spiralDistances <- function(pixelGroupMap, maxDis, cellSize) {
  spiral <- which(focalMat(pixelGroupMap, maxDis, type = "circle") > 0, arr.ind = TRUE) -
    ceiling(maxDis/cellSize) - 1
  spiral <- cbind(spiral, dists = sqrt( (0 - spiral[,1]) ^ 2 + (0 - spiral[, 2]) ^ 2))
  spiral <- spiral[order(spiral[, "dists"], apply(abs(spiral), 1, sum),
                         abs(spiral[, 1]), abs(spiral[, 2])),, drop = FALSE]
}
