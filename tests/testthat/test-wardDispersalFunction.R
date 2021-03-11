test_that("test Ward dispersal seeding algorithm", {

  verbose <- 0

  if (FALSE) {
    # devtools::load_all("~/GitHub/LandR")
    outSummary <- list()
  }
  library(LandR)
  library(data.table)
  library(raster)
  library(quickPlot)

  # keep this here for interactive testing with a larger raster
  doLarge <- if (interactive()) FALSE else FALSE
  if (doLarge) {
    set.seed(1234)
    print("Doing LARGE raster test -- should take more than 4 minutes")
    reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*18000,
                                   ymn = 50, ymx = 50 + 99*18000,
                                   res = c(250, 250), val = 2)
    pgs <- 10000
    proportionRcvCells <- 0.01
    if (FALSE) { # medium sized for interaactive use
      reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*300,
                                     ymn = 50, ymx = 50 + 99*300,
                                     res = c(250, 250), val = 2)
      proportionRcvCells <- 0.5
      pgs <- 30

    }
  } else {
    reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*25,
                                   ymn = 50, ymx = 50 + 99*25,
                                   res = c(100, 100), val = 2)
    pgs <- 30
    proportionRcvCells <- 0.5
  }

  reducedPixelGroupMap <- SpaDES.tools::randomPolygons(reducedPixelGroupMap, numTypes = pgs)
  Sum_of_species <- raster(reducedPixelGroupMap)
  td <- tempdir()
  speciesTable <- reproducible::Cache(getSpeciesTable, dPath = td)
  speciesTable <- speciesTable[Area == "BSW"]
  speciesTable[, speciesCode := as.factor(LandisCode)]
  speciesTable[, seeddistance_eff := SeedEffDist]
  speciesTable[, seeddistance_max := SeedMaxDist]
  speciesTable[seeddistance_max > 2000, seeddistance_max := 2000]


  for (jj in 1:2) {
    species <- data.table::copy(speciesTable)
    if (jj == 1) {
      rcvSpByPG <- lapply(seq_len(pgs * proportionRcvCells), function(pg) {
        data.table(speciesCode = sample(c(1, 3:11), size = sample(1:5, 1)))
      })
      srcSpByPG <- lapply(seq_len(pgs * (1 - proportionRcvCells)), function(pg) {
        data.table(speciesCode = sample(c(1, 3:11), size = sample(1:5, 1)))
      })
      species <- data.table(species)[, speciesCode := seq_along(LandisCode)]
    } else {
      speciesCodes <- factor(sample(paste0("spp_",LETTERS[c(1, 3:12)])))
      rcvSpByPG <- lapply(seq_len(pgs * proportionRcvCells), function(pg) {
        data.table(speciesCode = sample(speciesCodes[c(1:4, 6:11)], size = sample(1:5, 1)))
      })
      srcSpByPG <- lapply(seq_len(pgs * (1 - proportionRcvCells)), function(pg) {
        data.table(speciesCode = sample(speciesCodes[c(1:4, 8:11)], size = sample(1:5, 1)))
      })
      species <- data.table(species)[, speciesCode := speciesCodes[seq_along(LandisCode)] ]
    }


    successionTimestep = 10
    seedReceive <- rbindlist(rcvSpByPG, idcol = "pixelGroup")
    seedSource <- rbindlist(srcSpByPG, idcol = "pixelGroup")
    seedSource[, pixelGroup := pixelGroup + pgs/2]

    seedReceiveFull <- species[seedReceive, on = "speciesCode"]
    # objects <- list("species" = species)
    # mb <- profvis::profvis(replicate(10,
    #    interval = 0.2,
    # set.seed(seedOuter)
    # print(seedOuter)

    st1 <- system.time({
      output <- LANDISDisp(dtRcv = seedReceiveFull, plot.it = FALSE,
                           dtSrc = seedSource,
                           speciesTable = species,
                           reducedPixelGroupMap,
                           verbose = 1,
                           successionTimestep = successionTimestep)
    })

    if (interactive()) {
      print(output[, .N, by = speciesCode])
      print(st1)
    }

    pixelName <- grep("pixelIn", names(output), value = TRUE)
    outputSum <- output[, list(speciesCode = sum(as.integer(speciesCode))), by = pixelName]
    Sum_of_species[outputSum[[pixelName]]] <- outputSum$speciesCode

    # Plotting
    a <- reducedPixelGroupMap[] %in% seedReceive$pixelGroup
    sum(a)

    library(SpaDES.tools)
    if (is.factor(seedReceiveFull$speciesCode)) {
      sppNames <- levels(seedReceiveFull$speciesCode)
      sps <- sppNames
      names(sps) <- sppNames
    } else {
      sps <- sort(unique(as.integer(seedReceiveFull$speciesCode)))
      sppNames <- paste0("Species", as.character(sps))
      names(sps) <- sppNames
    }

    spsOut <- lapply(sps, function(sp) {
      ras <- raster(reducedPixelGroupMap)
      ras3 <- rasterizeReduced(seedReceiveFull[speciesCode == sp], reducedPixelGroupMap, "speciesCode", "pixelGroup")
      ras[!is.na(ras3[])] <- 1

      ras3 <- rasterizeReduced(seedSource[speciesCode == sp], reducedPixelGroupMap, "speciesCode", "pixelGroup")
      ras[!is.na(ras3[])] <- 2

      out <- output[speciesCode == sp]
      ras[out$pixelIndex] <- 3
      levels(ras) <- data.frame(ID = 1:3, c("Receive", "Source", "Successful"))
      ras
    })

    bigDispersers <- species$speciesCode[which(species$SeedMaxDist == 5000 & species$SeedEffDist == 1000)]
    lapply(spsOut, function(x) table(x[]))
    if (interactive()) {
      dev()
      clearPlot()
      Plot(reducedPixelGroupMap, Sum_of_species, new = TRUE)
      Plot(spsOut, legendRange = 1:3)
    }
    # i <<- i + 1;

    expect_true(all(unique(output$speciesCode) %in% unique(seedReceiveFull$speciesCode)))
    expect_true(all(is.na(Sum_of_species[reducedPixelGroupMap[] > 15]))) # nothing regenerates in the pgs that don't have receive available

    # Test whether each pixelGroup has only the species that could have arrived there
    output[, pixelGroup := reducedPixelGroupMap[pixelIndex]]
    joined <- seedReceiveFull[output, on = "pixelGroup", allow.cartesian = TRUE]
    joinedTest <- joined[, all(i.speciesCode %in% speciesCode) , by = "pixelGroup"]
    expect_true(all(joinedTest$V1))

  }
  if (!doLarge) {
    env <- new.env()
    env$cellSize = res(Sum_of_species)[1]
    env$b = 0.01
    env$k = 0.95

    testDists <- list()

    for (dis in 1:12) {
      message("Working on ", dis * env$cellSize, " m")
      ras2 <- raster(reducedPixelGroupMap)
      ras2[] <- 2
      ras2[SpaDES.tools::middlePixel(ras2)] <- 1
      ras2[SpaDES.tools::middlePixel(ras2) + dis] <- 3
      seedReceive <- data.table(pixelGroup = 3, speciesCode = species$speciesCode)
      seedSource <- data.table(pixelGroup = 1, speciesCode = species$speciesCode)
      output <- lapply(1:100, function(x)  {
        LANDISDisp(dtRcv = seedReceive, plot.it = FALSE,
                   dtSrc = seedSource,
                   speciesTable = species,
                   pixelGroupMap = ras2,
                   verbose = FALSE,
                   successionTimestep = successionTimestep)
      })
      output <- rbindlist(output)
      joined <- output[species, on = "speciesCode"]
      tooFar <- pmax(res(ras2)[1], joined$seeddistance_max) < dis * res(ras2)[1]
      expect_true(all(is.na(joined$pixelIndex[tooFar])))

      testDists[[dis]] <- sapply(unique(output$speciesCode), function(spCode) {
        env$effDist <- unique(joined[speciesCode == spCode]$seeddistance_eff)
        env$maxDist <- unique(joined[speciesCode == spCode]$seeddistance_max)
        env$dis <- dis * env$cellSize
        dispersalProb <- eval(Ward, envir = env)
        dispersalProb = 1 - (1 - dispersalProb)^successionTimestep
        probOfThatNumber <- dbinom(x = NROW(output[speciesCode == spCode]), size = 100, prob = dispersalProb)
      })
    }
    tests <- unlist(testDists)
    # Fairly conservative test -- the number of tests that fail at p < 0.01 should be about 5% ... really, it should be 1%
    expect_true(sum(tests < 0.01) / length(tests) <= 0.1)

    # Where rcv can receive a species, but it doesn't exist in Src
    seedReceive <- data.table(pixelGroup = 3, speciesCode = species$speciesCode[1])
    seedSource <- data.table(pixelGroup = 1, speciesCode = species$speciesCode[2])
    output <- LANDISDisp(dtRcv = seedReceive, plot.it = FALSE,
                         dtSrc = seedSource,
                         speciesTable = species,
                         pixelGroupMap = ras2,
                         verbose = FALSE,
                         successionTimestep = successionTimestep)
    expect_true(NROW(output) == 0)
  }
})

test_that("test large files", {
  if (interactive()) {
    whichTest <- 0 # 0 for full test (slow), 1 (manual interactive) or 2 (medium)
    dp <- "~/tmp"
  } else {
    whichTest <- 2
    dp <- tempdir()
  }
  library(reproducible)
  library(quickPlot)
  dtSrc <- prepInputs(url = 'https://drive.google.com/file/d/1MHA3LeBuPJXRPkPDp33M6iJmNpw7ePZI/view?usp=sharing',
                        targetFile = "dtSrc.rds",
                        fun = "readRDS",
                        destinationPath = dp, overwrite = TRUE)
  dtRcv <- prepInputs(url = 'https://drive.google.com/file/d/1MHA3LeBuPJXRPkPDp33M6iJmNpw7ePZI/view?usp=sharing',
                      targetFile = "dtRcv.rds",
                      fun = "readRDS",
                      destinationPath = dp)
  pixelGroupMap <- prepInputs(url = 'https://drive.google.com/file/d/1MHA3LeBuPJXRPkPDp33M6iJmNpw7ePZI/view?usp=sharing',
                              targetFile = "pixelGroupMap.rds",
                              fun = "readRDS",
                              destinationPath = dp)
  speciesTable <- prepInputs(url = 'https://drive.google.com/file/d/1MHA3LeBuPJXRPkPDp33M6iJmNpw7ePZI/view?usp=sharing',
                             targetFile = "speciesTable.rds",
                             fun = "readRDS",
                             destinationPath = dp)
  seed <- 1234
  set.seed(seed)
  dtSrc1 <- data.table::copy(dtSrc)
  dtRcv1 <- data.table::copy(dtRcv)
  sppKeep <- unique(dtRcv1$speciesCode)
  dtSrc1 <- dtSrc1[speciesCode %in% sppKeep]
  dtRcv1 <- dtRcv1[speciesCode %in% sppKeep]

  speciesTable1 <- data.table::copy(speciesTable)
  speciesTable1 <- speciesTable1[speciesCode %in% sppKeep]


  if (whichTest == 1) { # 1 is for manual, interactive testing

    both <- dtSrc1[dtRcv1, on = c("speciesCode", "pixelGroup"), nomatch = 0]

    pixelsWithSrcAndRcv <- which(pixelGroupMap[] %in% both$pixelGroup)

    pixelsWithSrcAndRcv <- pixelsWithSrcAndRcv[diff(pixelsWithSrcAndRcv) == 1 & c(FALSE, diff(diff(pixelsWithSrcAndRcv) == 1) == 0)]
    pix <- pixelsWithSrcAndRcv[4]
    pixGr <- pixelGroupMap[pix]
    pixGrs <- pixelGroupMap[pix + (-1:1)]

    dtSrc1 <- dtSrc1[pixelGroup %in% pixGr] # Abie_bal
    dtRcv2 <- dtRcv1[pixelGroup %in% (pixGrs)]

    # verify
    rcv <- which(pixelGroupMap[] %in% dtRcv2$pixelGroup)
    src <- which(pixelGroupMap[] %in% dtSrc1$pixelGroup)
    expect_true(src %in% rcv) # src is one of the rcv
    expect_true(sum(diff(rcv) == 1) > 1) # there are 3 adjacent cells
  } else if (whichTest == 2) {
    # subsetting -- but it doesn't seem to work for final test
    dtRcv2 <- dtRcv1[, .SD[sample(NROW(.SD), size = min(NROW(.SD), 300))], by = "speciesCode"]
  } else {
    dtRcv2 <- dtRcv1
  }
  st <- system.time(out <- LANDISDisp(dtSrc = dtSrc1,
                    dtRcv = dtRcv2,
                    pixelGroupMap = pixelGroupMap,
                    successionTimestep = 1,
                    speciesTable = speciesTable1))

  clearPlot()
  spMap <- list()
  spMap$pixelGroupMap <- pixelGroupMap
  for (sppp in unique(out$speciesCode)) {
    spMap[[sppp]] <- SpaDES.tools::rasterizeReduced(fullRaster = pixelGroupMap,
                                                    unique(dtSrc[speciesCode == sppp], on = c("pixelGroup", "speciesCode")),
                                                    newRasterCols = c("seeddistance_eff"))
    receivable <- SpaDES.tools::rasterizeReduced(fullRaster = pixelGroupMap,
                                                 unique(dtRcv[speciesCode == sppp], on = c("pixelGroup", "speciesCode")),
                                                 newRasterCols = c("seeddistance_eff"))

    forest <- which(!is.na(pixelGroupMap[]))
    src <- which(!is.na(spMap[[sppp]][]))
    recvable <- which(!is.na(receivable[]))
    rcvd <- out[speciesCode == sppp]$pixelIndex

    spMap[[sppp]][forest] <- 0
    spMap[[sppp]][recvable] <- 2
    spMap[[sppp]][src] <- 1
    spMap[[sppp]][rcvd] <- 3
    spMap[[sppp]][intersect(src, rcvd)] <- 4

    levels(spMap[[sppp]]) <- data.frame(ID = 0:4, type = c("OtherForest", "Source", "Didn't receive", "Received", "Src&Rcvd"))
  }
  if (interactive()) {
    clearPlot()
    sp <- spMap[-1]
    Plot(sp, cols = "Set2")
  }

  rr <- apply(raster::stack(spMap)[[-1]][] + 1, 2, function(x) {
    oo <- tabulate(x)
    if (length(oo) < 5)
      oo <- c(oo, rep(0, 5 - length(oo)))
    oo
    })
  rownames(rr) <- raster::levels(spMap[[2]])[[1]][,"type"]
  rr <- t(rr)
  rr <- as.data.frame(rr)
  rr <- cbind(rr, propSrcRcved = round(rr[,5]/ (rr[,5]+rr[,2]), 5))
  speciesTable[,c(1,5)]
  if (!(whichTest %in% 1:2)) {
    # This is a weak test -- that is often wrong with small samples -- seems to only
    #   work with full dataset
    corr <- cor(speciesTable[match(rownames(rr), species)]$shadetolerance, rr[, "propSrcRcved" ],
              method = "spearman")
    expect_true(corr > 0.8)
  }
  print(rr)
  print(st)


})
