test_that("test Ward dispersal seeding algorithm", {
  library(data.table)
  library(raster)
  reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*300,
                                 ymn = 50, ymx = 50 + 99*300,
                                 res = c(100, 100), val = 2)
  reducedPixelGroupMap <- raster(xmn = 50, xmx = 50 + 99*25,
                                 ymn = 50, ymx = 50 + 99*25,
                                 res = c(100, 100), val = 2)

  pgs <- 30
  reducedPixelGroupMap <- SpaDES.tools::randomPolygons(reducedPixelGroupMap, numTypes = pgs)
  rcvSpByPG <- lapply(seq_len(pgs/2), function(pg) {
    data.table(speciesCode = sample(1:11, size = sample(1:5, 1)))
  })
  srcSpByPG <- lapply(seq_len(pgs/2), function(pg) {
    data.table(speciesCode = sample(1:11, size = sample(1:5, 1)))
  })

  successionTimestep = 10
  seedReceive <- rbindlist(rcvSpByPG, idcol = "pixelGroup")
  seedSource <- rbindlist(srcSpByPG, idcol = "pixelGroup")
  seedSource[, pixelGroup := pixelGroup + pgs/2]

  speciesTable <- getSpeciesTable(dPath = ".")
  speciesTable <- speciesTable[Area == "BSW"]
  speciesTable[, speciesCode := as.factor(LandisCode)]
  speciesTable[, seeddistance_eff := SeedEffDist]
  speciesTable[, seeddistance_max := SeedMaxDist]

  species <- speciesTable
  species <- data.table(species)[, speciesCode := seq_along(LandisCode)]
  seedReceiveFull <- species[seedReceive, on = "speciesCode"]
  objects <- list("species" = species)
  mb <- profvis::profvis(
     interval = 0.2,
    output <- LANDISDisp(dtRcv = seedReceiveFull, plot.it = FALSE,
                         dtSrc = seedSource,
                         species = species,
                         reducedPixelGroupMap,
                         verbose = FALSE,
                         successionTimestep = successionTimestep)
   )
  output[, .N, by = speciesCode]

  ras <- raster(reducedPixelGroupMap)
  pixelName <- grep("pixelIn", names(output), value = TRUE)
  outputSum <- output[, list(speciesCode = sum(speciesCode)), by = pixelName]
  ras[outputSum[[pixelName]]] <- outputSum$speciesCode
  if (FALSE) {
    dev()
    clearPlot()
    Plot(reducedPixelGroupMap, ras, new= TRUE, col = c("red", "blue"))
  }

  expect_true(all(unique(output$speciesCode) %in% unique(seedReceiveFull$speciesCode)))
  expect_true(all(is.na(ras[reducedPixelGroupMap[] > 15]))) # nothing regenerates in the pgs that don't have receive available

  # Test whether each pixelGroup has only the species that could have arrived there
  output[, pixelGroup := reducedPixelGroupMap[pixelIndex]]
  joined <- seedReceiveFull[output, on = "pixelGroup", allow.cartesian = TRUE]
  joinedTest <- joined[, all(i.speciesCode %in% speciesCode) , by = "pixelGroup"]
  expect_true(all(joinedTest$V1))

  env <- new.env()
  env$cellSize = res(ras)[1]
  env$k =
  env$b = 0.01
  env$k = 0.95

  testDists <- list()

  for (dis in 1:12) {
    message("Working on ", dis * env$cellSize, " m")
    ras2 <- raster(reducedPixelGroupMap)
    ras2[] <- 2
    ras2[middlePixel(ras)] <- 1
    ras2[middlePixel(ras) + dis] <- 3
    seedReceive <- data.table(pixelGroup = 3, speciesCode = 1:11)
    seedSource <- data.table(pixelGroup = 1, speciesCode = 1:11)
    output <- lapply(1:100, function(x)  {
      LANDISDisp(dtRcv = seedReceive, plot.it = FALSE,
                 dtSrc = seedSource,
                 species = species,
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

})

