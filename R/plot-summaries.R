.defineLeading <- function(x, leadingPercentage = 0.8, totalCol) {
  colID <- which(x[-length(x)] > (leadingPercentage * x[[totalCol]]))
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
  if (requireNamespace("qs2", quietly = TRUE)) {
    if (is.null(treeType)) {
      treeType <- data.frame(
        leading = as.integer(c(
          seq_len(length(treeSpecies[["Species"]])),
          paste0(length(treeSpecies[["Species"]]) + 1, seq_len(length(treeSpecies[["Species"]])))
        )),
        landcover = c(treeSpecies[["Species"]], paste0("Mixed_", treeSpecies[["Species"]])),
        leadingType = c(
          tolower(treeSpecies[["Type"]]),
          rep("mixed", length(treeSpecies[["Species"]]))
        ),
        stringsAsFactors = FALSE
      )
      treeType$newClass <- ifelse(treeType$leadingType == "deciduous", 1,
        ifelse(treeType$leadingType == "conifer", 0, 0.5)
      )
    }

    ## 1. for each rep within a scenario, calculate difference -->
    ##    if conifer to decid = 1, if decid to conifer = -1, otherwise 0
    ## 2. Create one single map of "proportion net conversion" sum of difference / Nreps
    allReps <- parallel::mclapply(1:Nreps, function(rep) {
      resultsDir <- file.path(outputDir, sprintf("rep%02d", rep))

      bothYears <- lapply(years, function(year) {
        cohortData <- resultsDir |>
          file.path(paste0("cohortData_year", year, ".qs2")) |>
          qs2::qs_read()
        pixelGroupMap <- resultsDir |>
          file.path(paste0("pixelGroupMap_year", year, ".tif")) |>
          rasterRead()

        cohortDataReduced <- cohortData[, list(sumBio = sum(B, na.rm = TRUE)),
          by = c("speciesCode", "pixelGroup")
        ]

        biomassStack <- .stack(lapply(treeSpecies[["Species"]], function(tSp) {
          message(paste0(
            "[", studyAreaName, "_", climateScenario, "]: creating biomass map for ",
            tSp, " in year ", year, " [rep ", rep, "]"
          ))
          r <- SpaDES.tools::rasterizeReduced(
            reduced = cohortDataReduced[speciesCode == tSp, ],
            fullRaster = pixelGroupMap,
            newRasterCols = "sumBio",
            mapcode = "pixelGroup"
          )
          r[is.na(r[])] <- 0
          r[is.na(pixelGroupMap)] <- NA
          return(r)
        }))
        names(biomassStack) <- treeSpecies[["Species"]]

        biomassDT <- data.table(pixelID = 1:ncell(biomassStack), biomassStack[])
        biomassDT[, totalBiomass := rowSums(.SD, na.rm = TRUE),
                  .SDcols = names(biomassDT)[names(biomassDT) != "pixelID"]]
        biomassDT <- biomassDT[totalBiomass != 0, ]
        biomassDT[, leading := apply(.SD, 1, defineLeading,
                                     leadingPercentage = leadingPercentage,
                                     totalCol = "totalBiomass"),
                  .SDcols = names(biomassDT)[names(biomassDT) != "pixelID"]]
        biomassDT <- merge(biomassDT, treeType[, c("leading", "newClass")])
        allPixels <- data.table(pixelID = 1:ncell(biomassStack))
        biomassDTfilled <- merge(allPixels, biomassDT, all.x = TRUE, by = "pixelID")
        leadingSpeciesRaster <- rasterRead(biomassStack)
        leadingSpeciesRaster[] <- biomassDTfilled[["newClass"]]

        leadingSpeciesRaster
      })
      names(bothYears) <- paste0("Year", years)

      # bothYearsStk <- .stack(c(bothYears[[2]], -bothYears[[1]]))
      bothYearsStk <- terra::sds(bothYears[[2]], -bothYears[[1]])

      # leadingStackChange <- sum(bothYearsStk, na.rm = TRUE)
      leadingStackChange <- terra::app(bothYearsStk, fun = "sum", na.rm = TRUE)

      stopifnot(all(
        min(leadingStackChange[], na.rm = TRUE) >= -1,
        max(leadingStackChange[], na.rm = TRUE) <= 1
      ))
      leadingStackChange[is.na(rasterToMatch)] <- NA

      leadingStackChange
    })
    names(allReps) <- paste0("rep", 1:Nreps)

    if (length(allReps) > 1) {
      allRepsStk <- .stack(allReps)
      meanLeadingChange <- mean(allRepsStk, na.rm = TRUE)
    } else {
      meanLeadingChange <- allReps[[1]]
    }
    meanLeadingChange <- mask(crop(meanLeadingChange, rasterToMatch), rasterToMatch)

    f_meanLeadingChange <- file.path(
      outputDir, paste0("leadingChange_", studyAreaName, "_", climateScenario, ".tif")
    )
    writeRaster(meanLeadingChange, filename = f_meanLeadingChange, overwrite = TRUE)

    f_meanLeadingChange_gg <- file.path(outputDir, "figures") |>
      reproducible::checkPath(create = TRUE) |>
      file.path(paste0("leadingChange_", studyAreaName, "_", climateScenario, ".png"))

    fig <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = meanLeadingChange, maxcell = 1e+06) +
      tidyterra::scale_fill_whitebox_c("bl_yl_rd", direction = -1) +
      ggplot2::facet_wrap(~lyr) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.direction = "horizontal",
        legend.position = "top"
      ) +
      ggplot2::labs(
        title = paste("Proportional change in leading species under", climateScenario),
        subtitle = "Red: conversion to conifer; Blue: conversion to deciduous."
      )

    ggplot2::ggsave(filename = f_meanLeadingChange_gg, fig, width = 12, height = 12)

    return(list(f_meanLeadingChange, f_meanLeadingChange_gg))
  }
}
