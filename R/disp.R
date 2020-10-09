adj2 <- function(pixelGroupMapVec, pixelGroupMap, potentialReceivers, numCols, numCells, effDist, maxDist, cellSize,
                 dispersalFn, k, b, successionTimestep, dtSrcShort,
                 speciesSrcRasterVecList, dists, spRcvCommCodesList, verbose) {
  setorderv(potentialReceivers, c("fromInit", "RcvCommunity"))
  cells <- potentialReceivers$fromInit
  cellsXY <- xyFromCell(pixelGroupMap, cells)
  srcSpeciesCodes <- seq_along(speciesSrcRasterVecList)
  srcSpeciesCodes <- srcSpeciesCodes[!unlist(lapply(speciesSrcRasterVecList, is.null))]
  ymin <- pixelGroupMap@extent@ymin
  xmin <- pixelGroupMap@extent@xmin
  rcvSpeciesByIndex <- speciesCodeFromCommunity(potentialReceivers$RcvCommunity)
  rcvSpeciesCodes <- sort(unique(unlist(rcvSpeciesByIndex)))
  if (!all(rcvSpeciesCodes %in% srcSpeciesCodes)) {
    keep <- unlist(lapply(rcvSpeciesByIndex, function(rsbi) {
      all(rsbi %in% srcSpeciesCodes)
    }))
    rcvSpeciesByIndex <- rcvSpeciesByIndex[keep]
    cellsXY <- cellsXY[keep, , drop = FALSE]
  } else {
    keep <- rep(TRUE, length(rcvSpeciesByIndex))
  }
  speciesTableInner <- dists
  speciesTableInner <- unique(speciesTableInner[, c("speciesCode", "effDist", "maxDist")], )
  maxSpCode <- max(speciesTableInner[, "speciesCode"])

  speciesTableInner2 <- lapply(seq_len(maxSpCode), function(ind) {
    hasRow <- (speciesTableInner[, "speciesCode"] %in% ind )
    if (any(hasRow)) {
      speciesTableInner[hasRow,]
    } else {
      as.matrix(t(rep(NA, NCOL(speciesTableInner))))
    }

  })
  speciesTableInner <- do.call(rbind, speciesTableInner2)

  out <- spiralSeedDispersal(cellCoords = cellsXY, rcvSpeciesByIndex = rcvSpeciesByIndex,
                             speciesTable = speciesTableInner,
                             speciesVectorsList = speciesSrcRasterVecList,
                             cellSize = cellSize, numCells = numCells, xmin = xmin,
                             ymin = ymin, numCols = numCols, b = b, k = k,
                             successionTimestep = successionTimestep,
                             verbose = verbose)
  colNum <- seq(ncol(out))
  names(colNum) <- paste0("spCode", seq(colNum))
  seedsArrivedList <- lapply(colNum, function(col) {
    a <- out[, col]
    cells1 <- cells[keep][a]
  })

  seedsArrived <- data.table(pixelIndex = do.call(c, seedsArrivedList),
                             speciesCode = unlist(lapply(seq(seedsArrivedList),
                                                         function(x) rep(x, length(seedsArrivedList[[x]])))))
  return(seedsArrived)
}
