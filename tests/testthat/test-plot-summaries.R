## One replicate where every pixel goes from conifer- to deciduous-leading.
local_leadingSpeciesRun <- function(env = parent.frame()) {
  outputDir <- withr::local_tempdir("plotLeadingSpecies_", .local_envir = env)
  repDir <- file.path(outputDir, "rep01")
  dir.create(repDir)

  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, vals = 1)
  pixelGroupMap <- terra::rast(rtm, vals = 1:4)
  leading <- c("2011" = "Pice_mar", "2100" = "Popu_tre")
  for (year in names(leading)) {
    cohortData <- data.table::data.table(pixelGroup = 1:4, speciesCode = leading[[year]], B = 100L)
    qs2::qs_save(cohortData, file.path(repDir, paste0("cohortData_year", year, ".qs2")))
    terra::writeRaster(pixelGroupMap, file.path(repDir, paste0("pixelGroupMap_year", year, ".tif")))
  }

  list(
    studyAreaName = "test", climateScenario = "CanESM5_SSP370", Nreps = 1L,
    years = c(2011L, 2100L), outputDir = outputDir, rasterToMatch = rtm,
    treeSpecies = data.table::data.table(
      Species = c("Pice_mar", "Popu_tre"), Type = c("Conifer", "Deciduous")
    )
  )
}

test_that("plotLeadingSpecies() saves the figure to figurePath", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("tidyterra")
  withr::local_options(mc.cores = 1L)

  args <- local_leadingSpeciesRun()
  figDir <- file.path(withr::local_tempdir(), "figs")
  out <- suppressMessages(do.call(plotLeadingSpecies, c(args, figurePath = figDir)))

  expect_identical(file.exists(out), c(TRUE, TRUE))
  ## checkPath() normalizes the directory (`/private/var` on macOS, `/` on Windows)
  expect_identical(
    normalizePath(out[[2]], winslash = "/"),
    normalizePath(file.path(figDir, "leadingChange_test_CanESM5_SSP370.png"), winslash = "/")
  )
  expect_false(dir.exists(file.path(args$outputDir, "figures")))
})

test_that("plotLeadingSpecies() writes one leading-change layer, not one per species", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("tidyterra")
  withr::local_options(mc.cores = 1L)

  args <- local_leadingSpeciesRun()
  out <- suppressMessages(do.call(plotLeadingSpecies, args))

  leadingChange <- terra::rast(out[[1]])
  expect_identical(names(leadingChange), "leadingChange")
  expect_identical(as.vector(terra::values(leadingChange)), c(1, 1, 1, 1))
})

test_that("plotLeadingSpecies() saves the figure under outputDir by default", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("tidyterra")
  withr::local_options(mc.cores = 1L)

  args <- local_leadingSpeciesRun()
  out <- suppressMessages(do.call(plotLeadingSpecies, args))

  expect_identical(file.exists(out), c(TRUE, TRUE))
  expect_identical(
    normalizePath(out[[2]], winslash = "/"),
    normalizePath(file.path(args$outputDir, "figures", "leadingChange_test_CanESM5_SSP370.png"),
                  winslash = "/")
  )
})
