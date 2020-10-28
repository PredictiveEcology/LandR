library(LandR)
library(data.table)
library(raster)
library(testthat)
if (TRUE) {
  library(SpaDES.core)
  library(SpaDES.tools)
  #library(LandR)
  library(fpCompare)
  source(file.path("~/GitHub", "Biomass_core", "R", "seedDispersalLANDIS.R"))
  # keep this here for interactive testing with a larger raster
  module <- list("Biomass_core")
  path <- list(modulePath="~/GitHub/",
               outputPath="~/GitHub/")
  mySim <- simInit(times = list(start = 0, end = 2),
                   #modules = module,
                   paths = path)
}
# keep this here for interactive testing with a larger raster

# keep this here for interactive testing with a larger raster
# test_that("test Ward dispersal seeding algorithm", {

  verbose <- 1

  if (FALSE) {
    devtools::load_all("~/GitHub/LandR")
    outSummary <- list()
  }
  library(data.table)
  library(raster)
  library(quickPlot)

  # keep this here for interactive testing with a larger raster
  doLarge <- if (interactive()) TRUE else FALSE
  if (doLarge) {
    set.seed(1234)
    print("Doing LARGE raster test -- should take more than 4 minutes")
    reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*18000,
                                   ymn = 50, ymx = 50 + 99*18000,
                                   res = c(250, 250), val = 2)
    pgs <- 10000
    proportionRcvCells <- 0.01
  } else {
    reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*25,
                                   ymn = 50, ymx = 50 + 99*25,
                                   res = c(100, 100), val = 2)
    pgs <- 30
    proportionRcvCells <- 0.5
  }

  # seeds[[ff]] <- seedOuter # 639394 # 368697
  # set.seed(seedOuter)
  # saveRDS(seeds, file = "seed.rds")
  # saveRDS(verbose, file = "verbose.rds")
  reducedPixelGroupMap <- SpaDES.tools::randomPolygons(reducedPixelGroupMap, numTypes = pgs)
  Sum_of_species <- raster(reducedPixelGroupMap)
  rcvSpByPG <- lapply(seq_len(pgs * proportionRcvCells), function(pg) {
    data.table(speciesCode = sample(1:11, size = sample(1:5, 1)))
  })
  srcSpByPG <- lapply(seq_len(pgs * (1 - proportionRcvCells)), function(pg) {
    data.table(speciesCode = sample(1:11, size = sample(1:5, 1)))
  })

  successionTimestep = 10
  seedReceive <- rbindlist(rcvSpByPG, idcol = "pixelGroup")
  seedSource <- rbindlist(srcSpByPG, idcol = "pixelGroup")
  seedSource[, pixelGroup := pixelGroup + pgs/2]

  td <- tempdir()
  speciesTable <- reproducible::Cache(getSpeciesTable, dPath = td)
  speciesTable <- speciesTable[Area == "BSW"]
  speciesTable[, speciesCode := as.factor(LandisCode)]
  speciesTable[, seeddistance_eff := SeedEffDist]
  speciesTable[, seeddistance_max := SeedMaxDist]
   speciesTable[seeddistance_max > 2000, seeddistance_max := 2000]

  species <- speciesTable
  species <- data.table(species)[, speciesCode := seq_along(LandisCode)]
  seedReceiveFull <- species[seedReceive, on = "speciesCode"]
  # objects <- list("species" = species)
  # mb <- profvis::profvis(replicate(10,
  #    interval = 0.2,
  # set.seed(seedOuter)
  # print(seedOuter)
  objects <- list("species" = species)
  #    interval = 0.2,
  st1 <- system.time(
    output <- LANDISDisp(dtRcv = seedReceiveFull, plot.it = FALSE,
                       dtSrc = seedSource,
                       pixelGroupMap = reducedPixelGroupMap, species = species,
                       useParallel = FALSE, sim = mySim,
                       verbose = verbose,
                       successionTimestep = successionTimestep)

  #st1 <- system.time(
    # output <- LANDISDisp(dtRcv = seedReceiveFull, plot.it = FALSE,
    #                      dtSrc = seedSource,
    #                      speciesTable = species,
    #                      reducedPixelGroupMap,
    #                      verbose = 1,
    #                      successionTimestep = successionTimestep))
  )
  if (interactive()) {
    print(output[, .N, by = speciesCode])
    print(st1)
  }

  pixelName <- grep("pixelIn", names(output), value = TRUE)
  outputSum <- output[, list(speciesCode = sum(speciesCode)), by = pixelName]
  Sum_of_species[outputSum[[pixelName]]] <- outputSum$speciesCode
  if (FALSE) {
    a <- reducedPixelGroupMap[] %in% seedReceive$pixelGroup
    sum(a)

    library(SpaDES.tools)
    sps <- 1:11
    names(sps) <- paste0("Species", as.character(sps))
    sps <- lapply(sps, function(sp) {
      ras <- rasterizeReduced(seedReceiveFull[speciesCode == sp], reducedPixelGroupMap, "speciesCode", "pixelGroup")
      ras[!is.na(ras[])] <- 1

      ras2 <- rasterizeReduced(seedSource[speciesCode == sp], reducedPixelGroupMap, "speciesCode", "pixelGroup")
      ras[!is.na(ras2[])] <- 2

      out <- output[speciesCode == sp]
      ras[out$pixelIndex] <- 3
      levels(ras) <- data.frame(ID = 1:3, c("Receive", "Source", "Successful"))
      ras
    })
    dev()
    clearPlot()
    Plot(reducedPixelGroupMap, Sum_of_species, new = TRUE)
    Plot(sps, legendRange = 1:3)
    i <<- i + 1;

  }

  expect_true(all(unique(output$speciesCode) %in% unique(seedReceiveFull$speciesCode)))
  expect_true(all(is.na(Sum_of_species[reducedPixelGroupMap[] > 15]))) # nothing regenerates in the pgs that don't have receive available

  # Test whether each pixelGroup has only the species that could have arrived there
  output[, pixelGroup := reducedPixelGroupMap[pixelIndex]]
  joined <- seedReceiveFull[output, on = "pixelGroup", allow.cartesian = TRUE]
  joinedTest <- joined[, all(i.speciesCode %in% speciesCode) , by = "pixelGroup"]
  expect_true(all(joinedTest$V1))

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
      seedReceive <- data.table(pixelGroup = 3, speciesCode = 1:11)
      seedSource <- data.table(pixelGroup = 1, speciesCode = 1:11)
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
      tooFar <- joined$seeddistance_max < dis * res(ras2)[1]
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
    expect_true(sum(tests < 0.01)/length(tests) <= 0.05)

    # Where rcv can receive a species, but it doesn't exist in Src
    seedReceive <- data.table(pixelGroup = 3, speciesCode = 1)
    seedSource <- data.table(pixelGroup = 1, speciesCode = 2)
    output <-
      LANDISDisp(dtRcv = seedReceive, plot.it = FALSE,
                 dtSrc = seedSource,
                 speciesTable = species,
                 pixelGroupMap = ras2,
                 verbose = FALSE,
                 successionTimestep = successionTimestep)
    expect_true(NROW(output) == 0)
  }


# })


#
#
#
# objects <- list("species" = species)
# #    interval = 0.2,
# output <- LANDISDisp(dtRcv = seedReceiveFull, plot.it = FALSE,
#                      dtSrc = seedSource,
#                      pixelGroupMap = reducedPixelGroupMap, species = species,
#                      useParallel = TRUE, sim = mySim,
#                      verbose = verbose,
#                      successionTimestep = successionTimestep)
# output[, .N, by = speciesCode]
#
# pixelName <- grep("pixelIn", names(output), value = TRUE)
# outputSum <- output[, list(speciesCode = sum(speciesCode)), by = pixelName]
# ras[outputSum[[pixelName]]] <- outputSum$speciesCode
# if (FALSE) {
#   dev()
#   clearPlot()
#   Plot(reducedPixelGroupMap, ras, new = TRUE, col = c("red", "blue"))
# }
#
# expect_true(all(unique(output$speciesCode) %in% unique(seedReceiveFull$speciesCode)))
# expect_true(all(is.na(ras[reducedPixelGroupMap[] > 15]))) # nothing regenerates in the pgs that don't have receive available
#
# # Test whether each pixelGroup has only the species that could have arrived there
# output[, pixelGroup := reducedPixelGroupMap[pixelIndex]]
# joined <- seedReceiveFull[output, on = "pixelGroup", allow.cartesian = TRUE]
# joinedTest <- joined[, all(i.speciesCode %in% speciesCode) , by = "pixelGroup"]
# expect_true(all(joinedTest$V1))
#
# if (!doLarge) {
#   env <- new.env()
#   env$cellSize = res(ras)[1]
#   env$b = 0.01
#   env$k = 0.95
#
#   testDists <- list()
#
#   for (dis in 1:12) {
#     message("Working on ", dis * env$cellSize, " m")
#     ras2 <- raster(reducedPixelGroupMap)
#     ras2[] <- 2
#     ras2[SpaDES.tools::middlePixel(ras2)] <- 1
#     ras2[SpaDES.tools::middlePixel(ras2) + dis] <- 3
#     seedReceive <- data.table(pixelGroup = 3, speciesCode = 1:11)
#     seedSource <- data.table(pixelGroup = 1, speciesCode = 1:11)
#     output <- lapply(1:100, function(x)  {
#       LANDISDisp(dtRcv = seedReceive, plot.it = FALSE,
#                  dtSrc = seedSource,
#                  # pixelGroupMap = reducedPixelGroupMap,
#                  species = species,
#                  useParallel = TRUE, sim = mySim,
#                  pixelGroupMap = ras2,
#                  verbose = verbose,
#                  successionTimestep = successionTimestep)
#     })
#     output <- rbindlist(output)
#     joined <- output[species, on = "speciesCode"]
#     tooFar <- joined$seeddistance_max < dis * res(ras2)[1]
#     expect_true(all(is.na(joined$pixelIndex[tooFar])))
#
#     testDists[[dis]] <- sapply(unique(output$speciesCode), function(spCode) {
#       env$effDist <- unique(joined[speciesCode == spCode]$seeddistance_eff)
#       env$maxDist <- unique(joined[speciesCode == spCode]$seeddistance_max)
#       env$dis <- dis * env$cellSize
#       dispersalProb <- eval(Ward, envir = env)
#       dispersalProb = 1 - (1 - dispersalProb)^successionTimestep
#       probOfThatNumber <- dbinom(x = NROW(output[speciesCode == spCode]), size = 100, prob = dispersalProb)
#     })
#   }
#   tests <- unlist(testDists)
#   # Fairly conservative test -- the number of tests that fail at p < 0.01 should be about 5% ... really, it should be 1%
#   expect_true(sum(tests < 0.01)/length(tests) <= 0.1)
#
#   # Where rcv can receive a species, but it doesn't exist in Src
#   seedReceive <- data.table(pixelGroup = 3, speciesCode = 1)
#   seedSource <- data.table(pixelGroup = 1, speciesCode = 2)
#   output <-
#     LANDISDisp(dtRcv = seedReceive, plot.it = FALSE,
#                dtSrc = seedSource,
#                # pixelGroupMap = reducedPixelGroupMap,
#                species = species,
#                useParallel = TRUE, sim = mySim,
#                pixelGroupMap = ras2,
#                verbose = verbose,
#                successionTimestep = successionTimestep)
#   expect_true(NROW(output) == 0)
# }
#
#
#
