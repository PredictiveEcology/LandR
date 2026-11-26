## Replicate folders in which every pixel's leading species goes from the 2011 one to the 2100 one;
## by default one replicate, from conifer- to deciduous-leading.
local_leadingSpeciesRun <- function(env = parent.frame(), leadingByRep = NULL) {
  if (is.null(leadingByRep)) {
    leadingByRep <- list(rep01 = c("2011" = "Pice_mar", "2100" = "Popu_tre"))
  }
  outputDir <- withr::local_tempdir("plotLeadingSpecies_", .local_envir = env)

  rtm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, vals = 1)
  pixelGroupMap <- terra::rast(rtm, vals = 1:4)
  for (repName in names(leadingByRep)) {
    repDir <- file.path(outputDir, repName)
    dir.create(repDir)
    leading <- leadingByRep[[repName]]
    for (year in names(leading)) {
      cohortData <- data.table::data.table(pixelGroup = 1:4, speciesCode = leading[[year]],
                                           B = 100L)
      qs2::qs_save(cohortData, file.path(repDir, paste0("cohortData_year", year, ".qs2")))
      terra::writeRaster(pixelGroupMap,
                         file.path(repDir, paste0("pixelGroupMap_year", year, ".tif")))
    }
  }

  list(
    studyAreaName = "test", climateScenario = "CanESM5_SSP370", Nreps = length(leadingByRep),
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

test_that("plotLeadingSpecies() averages over the replicates in `reps` only", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("tidyterra")
  withr::local_options(mc.cores = 1L)

  ## rep02 and rep04 go from conifer- to deciduous-leading (+1); rep01 and rep03 the other way (-1)
  toDeciduous <- c("2011" = "Pice_mar", "2100" = "Popu_tre")
  toConifer <- c("2011" = "Popu_tre", "2100" = "Pice_mar")
  args <- local_leadingSpeciesRun(leadingByRep = list(
    rep01 = toConifer, rep02 = toDeciduous, rep03 = toConifer, rep04 = toDeciduous
  ))
  expect_identical(args$Nreps, 4L) # Biomass_summary passed max(reps) for reps = c(2, 4)

  ## all four replicates: the two directions cancel out
  out <- suppressMessages(do.call(plotLeadingSpecies, args))
  expect_identical(as.vector(terra::values(terra::rast(out[[1]]))), c(0, 0, 0, 0))

  ## rep01 and rep03 left out (a separate file name, so the first map is not overwritten)
  args24 <- utils::modifyList(args, list(studyAreaName = "reps24", reps = c(2L, 4L)))
  out24 <- suppressMessages(do.call(plotLeadingSpecies, args24))
  expect_identical(basename(out24[[1]]), "leadingChange_reps24_CanESM5_SSP370.tif")
  expect_identical(as.vector(terra::values(terra::rast(out24[[1]]))), c(1, 1, 1, 1))
})

test_that("plotLeadingSpecies() checks `reps` and that their folders exist", {
  skip_if_not_installed("qs2")

  args <- local_leadingSpeciesRun()
  for (badReps in list(0L, c(1L, 1L), 1.5, NA_integer_, Inf, integer(0), "1")) {
    expect_error(
      do.call(plotLeadingSpecies, c(args, list(reps = badReps))),
      "unique positive whole numbers",
      info = deparse(badReps)
    )
  }
  expect_error(do.call(plotLeadingSpecies, c(args, list(reps = 1:2))), "folders not found: .*rep02")
})
