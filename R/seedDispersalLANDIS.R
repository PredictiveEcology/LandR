utils::globalVariables(c(
  "..colnamesST", "RcvCommunity"
))

#' Simulate a LANDIS-II dispersal process on a landscape.
#'
#' Simulate seed dispersal using user defined function. This is a "receiving
#' pixel" focused dispersal approach. It is the "potentially receiving" cell
#' that looks around itself for potential seed sources. If it finds a single
#' seed source, that passes the probability function described by the
#' \code{dispersalFn}. If this passes a comparison to a uniform random draw,
#' then the receiving cell is deemed to have a "successful" dispersal for that
#' species.
#' This function can therefore only be used for a relatively specific situation
#' where there is a yes/no returned for each potential receiving cell, i.e., not
#' abundance. This function is also not cumulative, i.e,. there is no higher
#' abundance of seeds received if a receiving cell has lots of seed sources
#' around it vs. a single seed source. The difference will come with a higher
#' probability of successfully receiving a "seed".
#'
#' \code{dispersalFn} (temporarily unused as code is converted to Rcpp -- the
#' default \code{dispersalFn} is hard coded within the \code{spiralSeedDispersal}
#' function that uses C++) must be an expression that returns a probability
#' distribution. Because it is a dispersal kernal, it must be a probability
#' distribution. The expression that can take an argument named "dis" (without
#' quotes) as this will be calculated internally and represents the distance
#' from the initial (receiving) pixel and all active pixels within that cluster
#' of active pixels. \code{SpaDES} includes the \code{\link{Ward}} kernel as
#' defined in the LANDIS-II documentation.
#'
#' @param dtSrc data.table
#'
#' @param dtRcv data.table
#'
#' @param pixelGroupMap map
#'
#' @param speciesTable A data.table that should have at least 3 columns:
#'   \code{speciesCode} an integer representation of species, \code{seeddistance_max}
#'   a numeric with the maximum seed dispersal distance and \code{seeddistance_eff}
#'   the "effective" seed dispersal distance. These latter two are
#'   parameters passed to Ward dispersal kernel. This data.table can come from
#'   a \code{species} table from a LANDIS project.
#'
#' @param dispersalFn  An expression that can take a "dis" argument. See details.
#'   Default is "Ward" (temporarily unused, as it is hard coded inside Rcpp function)
#'
#' @param plot.it  If TRUE, then plot the raster at every interaction, so one can watch the
#' LANDISDisp event grow.
#' @param b  Landis Ward seed dispersal calibration coefficient (set to 0.01 in Landis)
#'
#' @param k  Landis Ward seed dispersal the probability that seed will disperse within
#' the effective distance (eg., 0.95)
#'
#' @param successionTimestep integer. The time in timeunits between succession (i.e., dispersal) events.
#'
#' @param verbose Numeric. \code{0} is not verbose, with increasing numbers indicating
#'   increasing levels of verbosity (currently up to 2)
#'
#' @param ...   Additional parameters. Currently none
#'
#' @return A numeric vector of raster pixel indices, in the same resolution and extent as
#' \code{seedSrc} raster.
#'
#' @importFrom R.utils intToBin
#' @importFrom magrittr %>%
#' @importFrom raster xyFromCell
#' @importFrom stats na.omit
#' @export
#' @docType methods
#'
#' @author Eliot McIntire
#'
#' @name LANDISDisp
#' @rdname LANDISDisp
#'
#' @examples
#' seed <- sample(1e6, 1)
#' seed <- 532597
#' set.seed(seed)
#' library(data.table)
#' library(raster)
#' # keep this here for interactive testing with a larger raster
#' rasterTemplate <- raster(extent(0, 2500, 0, 2500), res = 100)
#'
#' # make a pixelGroupMap
#' pgs <- 4 # make even just because of approach below requires even
#' pixelGroupMap <- SpaDES.tools::randomPolygons(rasterTemplate, numTypes = pgs)
#'
#' # Make a receive pixels table -- need pixelGroup and species
#' nSpecies <- 3
#' maxNSpeciesPerPixel <- min(5, nSpecies)
#' rcvSpByPG <- lapply(seq_len(pgs/2), function(pg) {
#'   data.table(speciesCode = sample(nSpecies, size = sample(maxNSpeciesPerPixel, 1)))
#' })
#' seedReceive <- rbindlist(rcvSpByPG, idcol = "pixelGroup")
#'
#' # Make a source pixels table -- need pixelGroup and species
#' srcSpByPG <- lapply(seq_len(pgs/2), function(pg) {
#'   data.table(speciesCode = sample(nSpecies, size = sample(maxNSpeciesPerPixel, 1)))
#' })
#' seedSource <- rbindlist(srcSpByPG, idcol = "pixelGroup")
#' # make source pixels not same pixelGroups as receive
#' seedSource[, pixelGroup := pixelGroup + pgs/2]
#'
#' # Get a species table -- if using in Canada, can use this
#' speciesTable <- getSpeciesTable(dPath = tempdir())
#' speciesTable <- speciesTable[Area == "BSW"]
#' speciesTable[, speciesCode := as.factor(LandisCode)]
#' speciesTable[, seeddistance_eff := SeedEffDist]
#' speciesTable[, seeddistance_max := SeedMaxDist]
#'
#' speciesTable <- speciesTable
#' speciesTable <- data.table(speciesTable)[, speciesCode := seq_along(LandisCode)]
#' seedReceiveFull <- speciesTable[seedReceive, on = "speciesCode"]
#' output <- LANDISDisp(dtRcv = seedReceiveFull, plot.it = interactive(),
#'                      dtSrc = seedSource,
#'                      speciesTable = speciesTable,
#'                      pixelGroupMap,
#'                      verbose = TRUE,
#'                      successionTimestep = 10)
#' # Summarize
#' output[, .N, by = speciesCode]
#'
LANDISDisp <- function(dtSrc, dtRcv, pixelGroupMap, speciesTable,
                       dispersalFn = WardFast, b = 0.01, k = 0.95, plot.it = FALSE,
                       successionTimestep,
                       verbose = getOption("LandR.verbose", TRUE),
                       ...) {
  if (TRUE) { # This is rewrite and MASSIVE simplification for spiralSeedDispersal
    # Setup Rcv components cellCoords and rcvSpeciesByIndex

    pgv <- pixelGroupMap[]

    # speciesSrcRasterVecList
    rasTemplate <- rep(NA_integer_, ncell(pixelGroupMap))
    srcSpeciesCodes <- sort(unique(dtSrc$speciesCode))
    names(srcSpeciesCodes) <- as.character(srcSpeciesCodes)
    cellsCanSrc <- which(pgv %in% dtSrc$pixelGroup)
    dt <- data.table(pixelGroup = pgv[cellsCanSrc], pixelIndex = cellsCanSrc)
    dt <- dtSrc[, c("pixelGroup", "speciesCode")][dt, on = "pixelGroup", allow.cartesian = TRUE] # $speciesCode
    srcSpeciesByIndex <- split(dt$pixelIndex, dt$speciesCode)
    speciesSrcRasterVecList <- lapply(srcSpeciesCodes, function(sc) {
      rasTemplate[srcSpeciesByIndex[[as.character(sc)]]] <- sc; rasTemplate
      })
    maxSpCode <- max(as.integer(names(speciesSrcRasterVecList)))
    speciesSrcRasterVecList <- lapply(seq_len(maxSpCode), function(ind) {
      if (as.character(ind) %in% names(speciesSrcRasterVecList))
        speciesSrcRasterVecList[[as.character(ind)]]
    })
    speciesVectorsList <- speciesSrcRasterVecList

    # # speciesVectorsList
    # # dtSrcShort <- dtSrc$pixelGroup
    # dtSrcNoDups <- unique(dtSrc, by = c("speciesCode"))
    # # setorderv(dtSrcNoDups, c("speciesCode", "pixelGroup"))
    # # rasTemplate <- raster(pixelGroupMap)
    # rasTemplate <- rep(NA_integer_, ncell(pixelGroupMap))
    # pixelGroupMapVec <- pixelGroupMap[]
    # spCodes <- sort(dtSrcNoDups$speciesCode)
    # names(spCodes) <- as.character(spCodes)
    # speciesSrcRasterVecList <- lapply(spCodes, function(sp) {
    #   ras <- rasTemplate
    #   hasSp <- dtSrc$speciesCode == sp
    #   # dtSrcSp <- dtSrc[hasSp]
    #   pgs <- dtSrc$pixelGroup[hasSp]
    #   pixels <- pixelGroupMapVec %in% pgs
    #   ras[pixels] <- sp
    #   ras
    # })

    # Raster metadata
    e <- pixelGroupMap@extent
    ymin <- e@ymin
    xmin <- e@xmin
    numCols <- pixelGroupMap@ncols
    numCells <- pixelGroupMap@ncols * pixelGroupMap@nrows #ncell(pixelGroupMap)
    numRows <- numCells / numCols
    cellSize <- xr <- (e@xmax - e@xmin)/numCols # cellSize <- res(pixelGroupMap) %>% unique()

    if (length(cellSize) > 1) {
      ## check for equal cell sizes that "aren't" due to floating point error
      res <-
        vapply(cellSize, function(x)
          isTRUE(all.equal(x, cellSize[1])), logical(1))
      if (all(res))
        cellSize <- cellSize[1]
      else
        stop("pixelGroupMap resolution must be same in x and y dimension")
    }

    # rcvSpeciesCodes <- sort(unique(unlist(rcvSpeciesByIndex)))
    srcSpeciesCodes <- seq_along(speciesSrcRasterVecList)
    srcSpeciesCodes <- srcSpeciesCodes[!unlist(lapply(speciesSrcRasterVecList, is.null))]

    # rcvSpeciesByIndex
    #  Remove any species in dtRcv that are not available in dtSrc
    dtRcv <- dtRcv[unique(dtSrc[, "speciesCode"], by = "speciesCode"), on = "speciesCode"]
    cellsCanRcv <- which(pgv %in% dtRcv$pixelGroup)
    rcvSpeciesCodes <- sort(unique(dtRcv$speciesCode))
    dt <- data.table(pixelGroup = pgv[cellsCanRcv], pixelIndex = cellsCanRcv)
    dtRcvSmall <- dtRcv[, c("pixelGroup", "speciesCode")]
    dtSrcUniqueSP <- unique(dtSrc[, "speciesCode"], by = "speciesCode")
    dtRcvSmall1 <- dtRcvSmall[dtSrcUniqueSP, on = "speciesCode"]
    dt <- dt[dtRcvSmall, on = "pixelGroup", allow.cartesian = TRUE] # $speciesCode
    setorderv(dt, c("pixelIndex", "speciesCode"))
    rcvSpeciesByIndex <- split(dt$speciesCode, dt$pixelIndex)

    # cellCoords
    cellCoords <- xyFromCell(pixelGroupMap, cellsCanRcv)

    # Removing cases ## 2 stages
    # 1st stage -- keep is for "keeping" only rcv pixels where at least 1 species is in the src pixels
    if (!all(rcvSpeciesCodes %in% srcSpeciesCodes)) {
      keep <- unlist(lapply(rcvSpeciesByIndex, function(rsbi) {
        any(rsbi %in% srcSpeciesCodes)
      }))
      rcvSpeciesByIndex <- rcvSpeciesByIndex[keep]

      cellCoords <- cellCoords[keep, , drop = FALSE]
    } else {
      keep <- rep(TRUE, length(rcvSpeciesByIndex))
    }

    # This 2nd stage removal is to remove individual cases where the more than one
    #   (but less than all -- was dealt with by keep) rcv species does not
    #   exist in any src pixel
    # rcvSpeciesByIndex1 <- lapply(rcvSpeciesByIndex, function(rsbi) {
    #   sort(rsbi[rsbi %in% srcSpeciesCodes])
    # })

    if (sum(keep) > 0) {
      maxDistColName <- grep("max", colnames(speciesTable), value = TRUE)
      effDistColName <- grep("eff", colnames(speciesTable), value = TRUE)
      colnamesST <- c("speciesCode", effDistColName, maxDistColName)
      speciesTableSmall <- speciesTable[, ..colnamesST]
      uniqueSTS <- unique(speciesTableSmall)
      speciesTableInner <- as.matrix(uniqueSTS)

      speciesTableInner2 <- lapply(seq_len(maxSpCode), function(ind) {
        hasRow <- (speciesTableInner[, "speciesCode"] %in% ind )
        if (any(hasRow)) {
          speciesTableInner[hasRow,]
        } else {
          as.matrix(t(rep(NA, NCOL(speciesTableInner))))
        }

      })
      speciesTableInner <- do.call(rbind, speciesTableInner2)
      speciesTableInner <- na.omit(speciesTableInner)

      # sam <- sample(1e6, 1)
      # print(sam)
      # set.seed(sam)
      ind <- seq(NROW(cellCoords))
      # ind <- sample(NROW(cellCoords), 15)
      out <- spiralSeedDispersal(cellCoords = cellCoords, #[ind,, drop = FALSE],
                                             rcvSpeciesByIndex = rcvSpeciesByIndex, #[ind],
                                             speciesTable = speciesTableInner,
                                             speciesVectorsList = speciesSrcRasterVecList,
                                             cellSize = cellSize, numCells = numCells, xmin = xmin,
                                             ymin = ymin, numCols = numCols, numRows = numRows, b = b, k = k,
                                             successionTimestep = successionTimestep,
                                             verbose = as.numeric(verbose))

      colNum <- seq(ncol(out))
      names(colNum) <- paste0("spCode", seq(colNum))
      seedsArrivedList <- lapply(colNum, function(col) {
        a <- out[, col]
        cells1 <- cellsCanRcv[keep][ind][a]
      })

      seedsArrived <- data.table(pixelIndex = do.call(c, seedsArrivedList),
                                 speciesCode = unlist(lapply(seq(seedsArrivedList),
                                                             function(x) rep(x, length(seedsArrivedList[[x]])))))
    } else {
      seedsArrived <- data.table(pixelIndex = integer(),
                                 speciesCode = integer())
    }

  } else {

    cellSize <- res(pixelGroupMap) %>% unique()
    if (length(cellSize) > 1) {
      ## check for equal cell sizes that "aren't" due to floating point error
      res <-
        vapply(cellSize, function(x)
          isTRUE(all.equal(x, cellSize[1])), logical(1))
      if (all(res))
        cellSize <- cellSize[1]
      else
        stop("pixelGroupMap resolution must be same in x and y dimension")
    }

    # NOTE new as.integer for speciesCode -- it is now a factor
    sc <- speciesTable[, list(
      speciesCode = as.integer(speciesCode),
      effDist = seeddistance_eff,
      maxDist = seeddistance_max
    )]
    dtSrc[, speciesCode := as.integer(speciesCode)]
    dtRcv[, speciesCode := as.integer(speciesCode)]

    setkey(sc, speciesCode)
    setkey(dtSrc, speciesCode)
    setkey(dtRcv, speciesCode)

    speciesSrcPool <-
      sc[dtSrc][, list(speciesSrcPool = sum(2 ^ speciesCode)), by = "pixelGroup"]
    setkeyv(speciesSrcPool, "pixelGroup")
    speciesSrcPool <- na.omit(speciesSrcPool)

    speciesRcvPool <-
      sc[dtRcv][, list(speciesRcvPool = sum(2 ^ speciesCode)), by = "pixelGroup"]
    setkeyv(speciesRcvPool, "pixelGroup")
    speciesRcvPool <- na.omit(speciesRcvPool)

    setkey(sc, speciesCode)
    spPool <- merge(speciesRcvPool, speciesSrcPool, all = TRUE)
    seedSourceMaps <- rasterizeReduced(
      spPool,
      fullRaster = pixelGroupMap,
      mapcode = "pixelGroup",
      newRasterCols = c("speciesSrcPool", "speciesRcvPool")
    )
    seedSourceMaps <-
      lapply(seedSourceMaps, function(x)
        setValues(x, as.integer(x[])))

    seedRcvOrig <- which(!is.na(seedSourceMaps$speciesRcvPool[]))
    seedSrcOrig <- which(seedSourceMaps$speciesSrcPool[] > 0)

    if (is.null(seedRcvOrig))  {
      # start it in the centre cell
      activeCell <- (nrow(seedSrcOrig) / 2L + 0.5) * ncol(seedSrcOrig)
    }
    ### should sanity check map extents
    if (plot.it) {
      wardSeedDispersalHab1 <- raster(seedSourceMaps$speciesSrcPool)
      wardSeedDispersalHab1[] <- NA
      Plot(wardSeedDispersalHab1, new = TRUE)
    }

    seedSourceMaps$speciesSrcPool <-
      NULL # don't need once xysAll are gotten

    potentialsOrig <- data.table(
      "fromInit" = seedRcvOrig,
      "RcvCommunity" = seedSourceMaps$speciesRcvPool[][seedRcvOrig],
      key = "fromInit"
    )

    if (verbose > 0)
      message("  Running seed dispersal")

    seedsArrived <- seedDispInnerFn(
      activeCell = seedRcvOrig, #theList$activeCell,# subSampList[[y]][[1]],
      potentials = potentialsOrig,#theList$potentials,# subSampList[[y]][[2]],
      n = cellSize,
      speciesRcvPool,
      sc,
      pixelGroupMap,
      cellSize,
      dtSrc,
      dispersalFn,
      k,
      b,
      # speciesComm,
      successionTimestep = successionTimestep,
      verbose = verbose
    )

    # COnvert speciesCode back to factor using original speciesTable object from top of this fn
    seedsArrived[, speciesCode := speciesTable$speciesCode[speciesCode]]
  }

  # setnames(seedsArrived, "fromInit", "pixelIndex")
  return(seedsArrived)
}

speciesCodeFromCommunity <- function(num) {
  indices <- lapply(strsplit(R.utils::intToBin(num), split = ""), function(x) {
    rev(as.logical(as.numeric(x)))
  })

  speciesCode <- lapply(indices, function(x) (seq_len(length(x)) - 1)[x])
}



speciesComm <- function(num, sc) {
  # indices <- lapply(strsplit(R.utils::intToBin(num), split = ""), function(x) {
  #   rev(as.logical(as.numeric(x)))
  # })
  #
  # speciesCode <- lapply(indices, function(x) (seq_len(length(x)) - 1)[x])
  speciesCode <- speciesCodeFromCommunity(num)
  data.table(RcvCommunity = as.integer(rep(num, sapply(speciesCode, length))),
             speciesCode = unlist(speciesCode),
             key = "speciesCode")[!is.na(speciesCode)] %>%
    sc[.]
}


#' @importFrom fpCompare %<=%
WardEqn <- function(dis, cellSize, effDist, maxDist, k, b) {
  if (cellSize %<=% effDist) {
    ifelse(
      dis %<=% effDist,
      exp((dis - cellSize) * log(1 - k) / effDist) -
        exp(dis * log(1 - k) / effDist),
      (1 - k) * exp((dis - cellSize - effDist) * log(b) / maxDist) -
        (1 - k) * exp((dis - effDist) * log(b) / maxDist)
    )
  } else {
    ifelse(
      dis %<=% cellSize,
      exp((dis - cellSize) * log(1 - k) / effDist) - (1 - k) *
        exp((dis - effDist) * log(b) / maxDist),
      (1 - k) * exp((dis - cellSize - effDist) * log(b) / maxDist) -
        (1 - k) * exp((dis - effDist) * log(b) / maxDist)
    )
  }
}


#' @inheritParams LANDISDisp
seedDispInnerFn <-
  function(activeCell, potentials, n, speciesRcvPool, sc,
           pixelGroupMap, cellSize, # xysAll,
           dtSrc, dispersalFn, k, b, # lociReturn,
           # speciesComm,
           successionTimestep,
           verbose = getOption("LandR.verbose", TRUE)) {

    spRcvCommCodes <-
      speciesComm(unique(speciesRcvPool$speciesRcvPool), sc = sc)
    setkey(spRcvCommCodes, RcvCommunity)
    setkey(potentials, RcvCommunity)

    numCells <- ncell(pixelGroupMap)
    numCols <- ncol(pixelGroupMap)

    # Build the speciesSrcRasterVecList -- the source object -- this will be a
    #   list where each element is a complete vector of all pixels, where NAs
    #   are "does not exist in that pixel" and otherwise it will be an integer
    #   of speciesCode
    dtSrcShort <- dtSrc$pixelGroup
    dtSrcNoDups <- unique(dtSrc, by = c("speciesCode"))
    speciesSrcRasterVecList <- by(dtSrcNoDups, INDICES = dtSrcNoDups$speciesCode, function(x)
      rasterizeReduced(x, pixelGroupMap, "speciesCode", "pixelGroup")[])
    speciesCodes <- as.character(dtSrcNoDups$speciesCode)
    names(speciesSrcRasterVecList) <- speciesCodes
    maxSpCode <- max(as.integer(names(speciesSrcRasterVecList)))
    speciesSrcRasterVecList <- lapply(seq_len(maxSpCode), function(ind) {
      if (as.character(ind) %in% names(speciesSrcRasterVecList))
        speciesSrcRasterVecList[[as.character(ind)]]
    })

    potentialReceivers <- potentials[, c("fromInit", "RcvCommunity")]
    speciesTableInner <- as.matrix(spRcvCommCodes)
    setorderv(potentialReceivers, c("fromInit", "RcvCommunity"))
    cells <- potentialReceivers$fromInit
    cellsXY <- xyFromCell(pixelGroupMap, cells)
    srcSpeciesCodes <- seq_along(speciesSrcRasterVecList)
    srcSpeciesCodes <- srcSpeciesCodes[!unlist(lapply(speciesSrcRasterVecList, is.null))]
    ymin <- pixelGroupMap@extent@ymin
    xmin <- pixelGroupMap@extent@xmin
    rcvSpeciesByIndex <- speciesCodeFromCommunity(potentialReceivers$RcvCommunity)
    rcvSpeciesCodes <- sort(unique(unlist(rcvSpeciesByIndex)))

    # Removing cases ## 2 stages
    # 1st stage -- keep is for "keeping" only rcv pixels where at least 1 species is in the src pixels
    if (!all(rcvSpeciesCodes %in% srcSpeciesCodes)) {
      keep <- unlist(lapply(rcvSpeciesByIndex, function(rsbi) {
        any(rsbi %in% srcSpeciesCodes)
      }))
      rcvSpeciesByIndex <- rcvSpeciesByIndex[keep]

      cellsXY <- cellsXY[keep, , drop = FALSE]
    } else {
      keep <- rep(TRUE, length(rcvSpeciesByIndex))
    }

    # This 2nd stage removal is to remove individual cases where the more than one
    #   (but less than all -- was dealt with by keep) rcv species does not
    #   exist in any src pixel
    rcvSpeciesByIndex <- lapply(rcvSpeciesByIndex, function(rsbi) {
      rsbi[rsbi %in% srcSpeciesCodes]
    })

    if (sum(keep) > 0) {
      speciesTableInner <- unique(speciesTableInner[, c("speciesCode", "effDist", "maxDist"), drop = FALSE], )

      speciesTableInner2 <- lapply(seq_len(maxSpCode), function(ind) {
        hasRow <- (speciesTableInner[, "speciesCode"] %in% ind )
        if (any(hasRow)) {
          speciesTableInner[hasRow,]
        } else {
          as.matrix(t(rep(NA, NCOL(speciesTableInner))))
        }

      })
      speciesTableInner <- do.call(rbind, speciesTableInner2)
      speciesTableInner <- na.omit(speciesTableInner)

      out <- spiralSeedDispersal(cellCoords = cellsXY,
                                 rcvSpeciesByIndex = rcvSpeciesByIndex,
                                 speciesTable = speciesTableInner,
                                 speciesVectorsList = speciesSrcRasterVecList,
                                 cellSize = cellSize, numCells = numCells, xmin = xmin,
                                 ymin = ymin, numCols = numCols, b = b, k = k,
                                 successionTimestep = successionTimestep,
                                 verbose = as.numeric(verbose))
      colNum <- seq(ncol(out))
      names(colNum) <- paste0("spCode", seq(colNum))
      seedsArrivedList <- lapply(colNum, function(col) {
        a <- out[, col]
        cells1 <- cells[keep][a]
      })

      seedsArrived <- data.table(pixelIndex = do.call(c, seedsArrivedList),
                                 speciesCode = unlist(lapply(seq(seedsArrivedList),
                                                             function(x) rep(x, length(seedsArrivedList[[x]])))))
    } else {
      seedsArrived <- data.table(pixelIndex = integer(),
                                 speciesCode = integer())
    }


    return(seedsArrived)

  }



#' Ward Seed Dispersal kernel
#'
#' A probability distribution used in LANDIS-II.
#'
#' @export
#' @docType methods
#'
#' @author Eliot McIntire
#'
#' @name Ward
#' @rdname Ward
Ward <- expression(if (cellSize <= effDist) {
  ifelse(
    dis <= effDist,
    exp((dis - cellSize) * log(1 - k) / effDist) -
      exp(dis * log(1 - k) / effDist),
    (1 - k) * exp((dis - cellSize - effDist) * log(b) / maxDist) -
      (1 - k) * exp((dis - effDist) * log(b) / maxDist)
  )
} else {
  ifelse(
    dis <= cellSize,
    exp((dis - cellSize) * log(1 - k) / effDist) - (1 - k) *
      exp((dis - effDist) * log(b) / maxDist),
    (1 - k) * exp((dis - cellSize - effDist) * log(b) / maxDist) -
      (1 - k) * exp((dis - effDist) * log(b) / maxDist)
  )
})

#' WardFast Seed Dispersal kernel
#'
#' A probability distribution used in LANDIS-II.
#'
#' @export
#' @docType methods
#'
#' @author Eliot McIntire
#'
#' @name WardFast
#' @rdname WardFast
WardFast <- expression(ifelse(cellSize <= effDist, {
  ifelse(
    dis <= effDist,
    exp((dis - cellSize) * log(1 - k) / effDist) -
      exp(dis * log(1 - k) / effDist),
    (1 - k) * exp((dis - cellSize - effDist) * log(b) / maxDist) -
      (1 - k) * exp((dis - effDist) * log(b) / maxDist)
  )
} , {
  ifelse(
    dis <= cellSize,
    exp((dis - cellSize) * log(1 - k) / effDist) - (1 - k) *
      exp((dis - effDist) * log(b) / maxDist),
    (1 - k) * exp((dis - cellSize - effDist) * log(b) / maxDist) -
      (1 - k) * exp((dis - effDist) * log(b) / maxDist)
  )
}))
