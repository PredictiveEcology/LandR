.defineLeading <- function(x, leadingPercentage = 0.8, totalCol) {
  colID <- which(x[-length(x)] > (leadingPercentage*x[[totalCol]]))
  if (length(colID) == 0) {
    # If we don't have a leading, we need to id conifer leading, or deciduous leading
    colID1 <- which.max(x[-length(x)])
    colID <- as.integer(paste0(length(x), colID1))
  }
  return(colID)
}

#' Leading species plots
#'
#' Plot effects on conifer-to-deciduous or deciduous-to-conifer conversions.
#'
#' @template summary_plots
#' @template Nreps
#' @param years TODO
#' @param treeSpecies TODO
#' @param defineLeading TODO
#' @param leadingPercentage TODO
#' @param treeType TODO
#' @template rasterToMatch
#'
#' @return list of filepaths corresponding to the images and/or objects written to disk
#'
#' @export
plotLeadingSpecies <- function(studyAreaName, climateScenario, Nreps, years, outputDir, treeSpecies,
                               defineLeading = .defineLeading, leadingPercentage = 0.8,
                               treeType = NULL, rasterToMatch) {
  if (requireNamespace("qs", quietly = TRUE)) {
    if (is.null(treeType)) {
      treeType <- data.frame(
        leading = as.integer(c(1:length(treeSpecies[["Species"]]),
                               paste0(length(treeSpecies[["Species"]]) + 1, 1:length(treeSpecies[["Species"]])))),
        landcover = c(treeSpecies[["Species"]], paste0("Mixed_", treeSpecies[["Species"]])),
        leadingType = c(tolower(treeSpecies[["Type"]]), rep("mixed", length(treeSpecies[["Species"]]))),
        stringsAsFactors = FALSE
      )
      treeType$newClass <- ifelse(treeType$leadingType == "deciduous", 1,
                                  ifelse(treeType$leadingType == "conifer", 0, 0.5))
    }

    # 1. for each rep within a scenario, calculate difference -->
    #    if conifer to decid = 1, if decid to conifer = -1, otherwise 0
    # 2. Create one single map of "proportion net conversion" sum of difference / Nreps
    allReps <- parallel::mclapply(1:Nreps, function(rep) {
      runName <- sprintf("%s_%s_run%02d", studyAreaName, climateScenario, rep)
      resultsDir <- file.path(outputDir, runName)

      bothYears <- lapply(years, function(year) {
        cohortData <- qs::qread(file = file.path(resultsDir, paste0("cohortData_", year, "_year", year, ".qs")))
        pixelGroupMap <- raster(file.path(resultsDir, paste0("pixelGroupMap_", year, "_year", year, ".tif")))

        cohortDataReduced <- cohortData[, list(sumBio = sum(B, na.rm = TRUE)),
                                        by = c("speciesCode", "pixelGroup")]

        biomassStack <- raster::stack(lapply(treeSpecies[["Species"]], function(tSp) {
          message(paste0("[", studyAreaName, "_", climateScenario, "]: creating biomass map for ",
                         tSp, " in year ", year, " [rep ", rep, "]"))
          r <- SpaDES.tools::rasterizeReduced(reduced = cohortDataReduced[speciesCode == tSp, ],
                                              fullRaster = pixelGroupMap,
                                              newRasterCols = "sumBio",
                                              mapcode = "pixelGroup")
          r[is.na(r[])] <- 0
          r[is.na(pixelGroupMap)] <- NA
          return(r)
        }))
        names(biomassStack) <- treeSpecies[["Species"]]

        biomassDT <- data.table(pixelID = 1:raster::ncell(biomassStack), raster::getValues(biomassStack))
        biomassDT[, totalBiomass := rowSums(.SD, na.rm = TRUE),
                  .SDcols = names(biomassDT)[names(biomassDT) != "pixelID"]]
        biomassDT <- biomassDT[totalBiomass != 0, ]
        biomassDT[, leading := apply(.SD, 1, defineLeading,
                                     leadingPercentage = leadingPercentage,
                                     totalCol = "totalBiomass"),
                  .SDcols = names(biomassDT)[names(biomassDT) != "pixelID"]]
        biomassDT <- merge(biomassDT, treeType[, c("leading","newClass")])
        allPixels <- data.table(pixelID = 1:raster::ncell(biomassStack))
        biomassDTfilled <- merge(allPixels, biomassDT, all.x = TRUE, by = "pixelID")
        leadingSpeciesRaster <- raster::setValues(raster(biomassStack), biomassDTfilled[["newClass"]])
        names(leadingSpeciesRaster) <- paste("biomassMap", studyAreaName, climateScenario, sep = "_")

        leadingSpeciesRaster
      })
      names(bothYears) <- paste0("Year", years)

      leadingStackChange <- raster::calc(raster::stack(bothYears[[2]], -bothYears[[1]]),
                                         fun = sum, na.rm = TRUE)

      stopifnot(all(minValue(leadingStackChange) >= -1, maxValue(leadingStackChange) <= 1))

      leadingStackChange[is.na(rasterToMatch)] <- NA
      names(leadingStackChange) <- paste("leadingMapChange", studyAreaName, climateScenario, rep, sep = "_")
      leadingStackChange
    })
    names(allReps) <- paste0("rep", 1:Nreps)

    fmeanLeadingChange <- file.path(outputDir, studyAreaName,
                                    paste0("leadingChange_", studyAreaName, "_", climateScenario, ".tif"))
    if (length(allReps) > 1) {
      meanLeadingChange <- raster::calc(raster::stack(allReps), mean, na.rm = TRUE)
    } else {
      meanLeadingChange <- allReps[[1]]
    }
    meanLeadingChange <- mask(crop(meanLeadingChange, rasterToMatch), rasterToMatch)
    writeRaster(meanLeadingChange, filename = fmeanLeadingChange, overwrite = TRUE)

    maxV <- max(abs(round(minValue(meanLeadingChange), 1)), abs(round(maxValue(meanLeadingChange), 1)))
    AT <- seq(-maxV, maxV, length.out = 12)

    pal <- RColorBrewer::brewer.pal(11, "RdYlBu")
    pal[6] <- "#f7f4f2"

    fmeanLeadingChange_gg <- file.path(outputDir, studyAreaName, "figures",
                                       paste0("leadingChange_", studyAreaName, "_", climateScenario, ".png"))

    fig <- rasterVis::levelplot(
      meanLeadingChange,
      main = paste("Proportional change in leading species under", climateScenario),
      sub = list(
        paste0(" Red: conversion to conifer\n",
               " Blue: conversion to deciduous."),
        cex = 2
      ),
      margin = FALSE,
      maxpixels = 7e6,
      at = AT,
      colorkey = list(
        space = "bottom",
        axis.line = list(col = "black"),
        width = 0.75
      ),
      par.settings = list(
        strip.border = list(col = "transparent"),
        strip.background = list(col = "transparent"),
        axis.line = list(col = "transparent")
      ),
      scales = list(draw = FALSE),
      col.regions = pal,
      par.strip.text = list(cex = 0.8, lines = 1, col = "black")
    )

    ## levelplot (trellis graphics more generally) won't plot correctly inside loop w/o print()
    png(filename = fmeanLeadingChange_gg, width = 1000, height = 1000)
    print(fig)
    dev.off()

    return(list(fmeanLeadingChange, fmeanLeadingChange_gg))
  }
}
