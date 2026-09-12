## Regenerate committed test fixtures under `tests/testthat/fixtures/`.
##
## Run manually; not part of the build. These exist so tests exercising the
## Canadian ecostratification polygons can run on CI: sis.agr.gc.ca throttles
## connections from CI runners (connection-level failure, not a cert or
## User-Agent problem), so downloading them at test time is not an option.

library(sf)

out_dir <- file.path("tests", "testthat", "fixtures")
stopifnot(dir.exists(out_dir))

td <- tempfile("ecodistricts_")
dir.create(td)

ecodistricts <- reproducible::prepInputs(
  url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
  destinationPath = td
) |>
  sf::st_make_valid()

## Ecoregion 88 (prairie/boreal SK) intersects the synthetic test raster;
## ecoregion 115 (Avalon Peninsula, NL) exists but does not -- keeping both lets
## the test cover the "no spatial intersection" branch. Retain the un-dissolved
## parts so `prep_polygons()` is genuinely exercised (the group/union step is
## part of the code path under test).
subset_ids <- c(88L, 115L)

## Simplify to a 1 km tolerance: the test raster has ~10 km cells, so this is
## well below the resolution anything under test can see, and it takes the
## committed fixture from ~143 KB to ~22 KB.
ecodistricts |>
  subset(ECOREGION %in% subset_ids, select = c("ECODISTRIC", "ECOREGION")) |>
  sf::st_transform(4326) |>
  sf::st_simplify(dTolerance = 1000) |>
  sf::st_make_valid() |>
  saveRDS(file.path(out_dir, "ecodistricts_88_115.rds"), compress = "xz")

message("wrote ", file.path(out_dir, "ecodistricts_88_115.rds"))
