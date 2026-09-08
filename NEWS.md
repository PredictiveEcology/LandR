# LandR (development version)

## Breaking changes

* drop support for R 4.2 due to changes in dependency packages;
* remove `rasterRead()` to use version from `reproducible`;
* `columnsForPixelGroups` is now a function (i.e., use `columnsForPixelGroups()` for consistent `pixelGroup` definitions);

## Dependency changes

* **now requires `reproducible (>= 3.1.1.9063)`** following non-backwards-compatible
  changes to the `reproducible` API: `options("reproducible.gdalwarp")` was removed,
  and `reproducible.inputPaths` was renamed to `reproducible.destinationPathShared`
  (the old name remains as a deprecated alias);
* remove deprecated package `crayon` in favour of `cli` instead;
* remove deprecated package `qs` in favour of `qs2` instead;
* add `knitr` and `rmarkdown` to Suggests for the new vignette, and declare both in
  `VignetteBuilder`. Listing `rmarkdown` there is what lets the vignette be rebuilt
  under `_R_CHECK_DEPENDS_ONLY_=true` (the `--as-cran` "no suggests" check), which
  otherwise hides it from the build even when installed;
* move `ggpubr` to Suggests;
* add `arrow` to Suggests;

## New features

* add functions to visualize vegetation type transitions;
* new function `cohortDefinitionCols()` to ensure consistent cohort definitions;
* new function `lccMapGenerator()` to calculate landcover classes from `cohortData` and `pixelGroupMap`;
* add new `prepInputs_NTEMS_DominantSpecies` function for importing dominant species layers from NTEMS website;
* add new `speciesPresentFromNTEMS` function to import dominant species layer from NTEMS and create factor raster to be hosted on Google drive;
* add `loadSCANFISpeciesLayers` and `prepSpeciesLayers_SCANFI` functions for loading SCANFI species data from Google drive;
* add `adjustAgeToLongevity` to adjust initial cohort ages based on `longevity` for each species;
* add `studyAreaEco` function to extend `studyArea` to ecological boundaries;
  - `studyAreaEco` allows `studyArea = NULL`, uses `type = "ecozone"` by default;
* add `plot_raster_stats` and `calc_raster_counts` for generating summaries of numeric rasters in Canada;
* new vignette `scanfi-summary`, summarizing a national stand age map by ecozone
  (closes #177). The rasters these summaries are designed for cannot be processed
  during a build -- SCANFI at 30 m is 2.1e10 pixels, a single pass reads tens of GB,
  CI runners have ~14 GB of disk, and the ecostratification polygon host throttles CI
  IP ranges -- so the vignette is split at the one expensive step. The value-frequency
  table from `calc_raster_counts()` is precomputed by `data-raw/precompute-vignettes.R`
  and committed under `inst/extdata/` (~640 KB, including a coarsened raster for maps
  and simplified polygons); every statistic, histogram and map is then recomputed from
  it when the vignette is built -- offline, in ~15 s. The committed table is pinned to
  a hash of the normalized source of the functions that produced it, so it cannot go
  stale unnoticed: `tests/testthat/test-vignette-artifacts.R` fails when the two
  diverge and names the script to re-run (that check is skipped under `covr`, which
  rewrites function bodies to insert trace counters, so the source it would hash is
  not the source that produced the artifacts). The summary statistics are deliberately
  excluded from that hash -- they are an exact function of the counts table, so
  changing the quantile definition shows up in the next build rather than forcing a
  multi-hour regeneration of a table that is still valid;

## Enhancements

* SCANFI download failures now explain themselves (closes #163). SCANFI is distributed
  through a Google Drive folder shared with collaborators, so a user without access saw only
  `reproducible`'s generic `Could not access the Google Drive resource ... (404) Not Found`.
  That reads like a dead link and invites a hunt for a public mirror that does not exist.
  `prepInputsStandAgeMap()`, `prepRawBiomassMap()`, `loadSCANFISpeciesLayers()`,
  `convert_SCANFI_LCC_codes()` and `prepInputs_SCANFI_LCC_FAO()` now report it as a
  permissions problem, point at <https://opendata.nfis.org/> to request access, and name the
  ways out -- supplying your own copy, or `dataSource = "KNN"` / `"NTEMS"` where the function
  offers them. The underlying error is still shown in full, and failures that are *not* access
  problems pass through untouched;

* `convertUnwantedLCC()` no longer uses the iterative `spread2()` search, whose run time grew
  with the square of the radius of the largest contiguous block of `classesToReplace`. On
  study areas containing large lakes/burns masked to an irregular boundary that search could
  take **hours, or never finish** (a 3.8 M-cell Western-Alberta study area: >2.8 h and
  unfinished). Each unwanted pixel now takes an available class drawn with probability
  proportional to that class's abundance within the smallest window reaching the pixel's
  nearest available class — the window at which `spread2()` would have stopped. Distances
  come from one vectorized `terra::distance()` transform per candidate class and the counts
  from a summed-area table, so the cost is independent of blob geometry: the same study area
  now completes in seconds. Per-ecoregion availability constraints are preserved exactly,
  and the abundance weighting reproduces the former search's class mix to within its own
  seed-to-seed variance. The exact output is *not* reproducible — the window radius now
  comes from a distance transform rather than from counting spread iterations — so runs
  needing pre-1.2.0.9004 output bit-for-bit must pin `LandR (<= 1.2.0.9003)`.
* `convertUnwantedLCC()` gains a `method` argument, which `overlayLCCs()` also accepts and
  passes through. Both options allocate identically and differ only in where the draw comes
  from. `method = "nearestWeighted"` (the default) keys it on the pixel's ground position, so
  it needs no `set.seed()`, is stable under `Cache()` (which does not key on RNG state), and
  — because the key is the cell centre rather than the cell index — a grid-aligned crop
  reproduces its parent raster cell for cell, so a small development subset agrees with the
  scaled-up run. `method = "nearestRandom"` draws from the RNG instead, for when replicates
  should differ.
* **`convertUnwantedLCC()`'s deterministic nearest-class rule, briefly present in
  1.2.0.9004, has been removed.** It broke distance ties to the lowest land-cover class, and
  ties turn out to be common — 35–41% of unwanted pixels on real landscapes — so it pulled
  systematically toward low-numbered classes. Under the Canada LCC coding those are the
  sparse, non-forest types: pooled over four real landscapes it assigned shrubs **1.69×** as
  often as the previous implementation, broadleaf **0.58×** and mixedwood **0.28×**, moving
  roughly one in fourteen unwanted pixels out of forest altogether. `"nearestWeighted"`
  gives the same determinism without that bias.
* `LANDISDisp()` spiral seed dispersal loop ported to C++ via `Rcpp`
  (~3.5–5.7× faster end-to-end depending on input size; ~5× on landscape-scale
  fixtures of 9 M cells). Memory use also drops dramatically: the
  per-cell-by-species source matrix is replaced with a per-`pixelGroup`
  species bitmask. Bit-identical to the previous R implementation under a
  fixed seed — guarded by 196 seed-locked / parity expectations (209 with
  `LANDR_SLOW_TESTS=1`). Default behaviour is unchanged for callers; the new
  path is on by default and can be opted out with
  `LANDISDisp(..., useCpp = FALSE)` or
  `options(LandR.LANDISDisp.useCpp = FALSE)`. `Rcpp (>= 1.0.10)` added to
  `Imports` and `LinkingTo`; `digest` added to `Suggests` (used by the
  golden-output hash manifest in tests);
* `calc_raster_stats()` now derives min/mean/max/quantiles/proportion-zero from
  the value-frequency table produced by `calc_raster_counts()` instead of
  reading every pixel via `exactextractr`, keeping memory bounded on very large
  rasters; it gains a `counts_df` argument to reuse a precomputed table.
  `plot_raster_stats()` now honours supplied `counts_df`/`stats_df` instead of
  recomputing them, and fetches the Canada inset once, omitting it gracefully
  when the source is unavailable. `zonal` removed from `Suggests`;
* `makeAndCleanInitialCohortData()` now `stop()`s with an informative message
  naming the species that cannot be age-imputed, instead of the cryptic
  `predict.merMod()` "non-conformable arguments": a species with cohorts needing
  age imputation but no usable known-age rows to fit the age model (all dropped
  for zero biomass and/or cover) is absent from the model fit yet present in the
  prediction set, so its fixed-effect `speciesCode` has no coefficient. A `TODO`
  in the code outlines fix options (see #195);
* SCANFI and 2020 now default data source and year for stand age and biomass functions;
* use `writeTo` instead of `filename2` in `prepInputs()` and related calls, following changes in `reproducible`;
* `minRelativeB` defaults updated based on discussion surrounding over-representation of shade tolerant species establishing and generating unreasonably high levels of understory cohorts;
* update `speciesInStudyArea` function to create `dataSource` parameter to direct function to download KNN or NTEMS factor raster from google drive and create associated species list;
* update `prepRawBiomassMap` function to allow for incorporation of NTEMS or SCANFI biomass;
* update `prepInputsStandAgeMap` function to allow for incorporation of SCANFI age map;
* update documentation and citations for `prepSpeciesLayers_*` functions;
* standardized `sppEquivalencies_CA` naming convention for provincial forestry columns with `<province>_forestry` ;
* remove undifferentiated tree species variants from provincial forestry columns in `sppEquivalencies_CA`;
* `plotVTM` now does not use `Plot` internally (with #140);
* several minor updates to `loadSCANFISpeciesLayers`, `prepSpeciesLayers_SCANFI` to address more edge cases;
* `prepSpeciesLayers_SCANFI` updates to improve join `sppEquiv` so "multiple - to - one" can be used;
* improved transition plots, use `arrow` datasets to minimize memory use (important for large study areas);
* use `.scanfi_v1_years` and `.scanfi_v2_years` instead of harcoded years in multiple places;
* `convert_SCANFI_LCC_codes()` now writes an `INT1U`, 256x256-tiled, LZW-compressed
  GeoTIFF (`NAflag = 255`) and gains a `writeTo` argument to name the output. The recode is
  a bijection whose every output code lies in 20-230, so one unsigned byte suffices --- the
  same footprint as the `Byte` SCANFI source. `terra`'s default of `Float32` in full-width
  strips instead quadrupled the per-pixel cost and forced every crop of these national
  178400 x 119100 rasters to read entire 178400-pixel rows. All 12 precomputed layers on
  Google Drive (V2 1985-2025, V1 2000/2010/2020) were rebuilt and verified against their
  distributed copies: **12/12 identical**, with zero value mismatches and zero NA-pattern
  mismatches over all 21,247,440,000 pixels each, and per-class histograms equal (which for
  V2 2020 also reproduce the counts published in NRCan's own `.tif.aux.xml` sidecar). Total
  storage falls from 32.98 GB to 15.05 GB (~2.19x). Values are unchanged, so a rebuilt layer
  is a drop-in replacement for the precomputed copy it replaces.
* `plot_raster_stats()` gains `inset_polygons`, the outline drawn behind the inset map.
  It previously always downloaded the Canada boundary through `gadm_canada()`, so the
  inset could not be drawn offline -- passing e.g. `sf::st_union()` of the polygons
  being summarized keeps the call local;
* `plot_raster_stats()` gains `fig_width`, `fig_height` and `fig_dpi`. The saved figure
  size was hard-coded at 12 x 7.5 in and 300 dpi, which is ~20x larger than a web page
  or vignette needs;
* `plot_raster_stats()`'s per-polygon figure assembly is factored into internal panel
  builders (`.hist_panel()`, `.map_panel()`, `.inset_panel()`, `.stats_panel()`,
  `.write_polygon_figure()`), taking a 120-line loop body down to one call;
* the `calc_raster_stats()`/`plot_raster_stats()` tests no longer skip on CI. They
  previously downloaded the ecodistrict polygons from sis.agr.gc.ca, which throttles CI
  runners; they now read a committed 23 KB fixture (`data-raw/make-test-fixtures.R`
  regenerates it), so the full plotting path is exercised on every commit. A separate
  test, still skipped on CI, checks the upstream source is reachable;

## Bug fixes

* `prepSpeciesLayers_SCANFI()`: the Google Drive fallback (taken when `RCurl::url.exists()`
  fails, e.g. during a network blip) referenced `year`, which is not a formal, so it
  resolved to `data.table::year` and failed with "cannot coerce type 'closure' to
  vector of type 'character'". It now uses `dataYear`.

* `assertERGs()` now gives an informative error when `ecoregionMap` carries no
  `ecoregionGroup` values (e.g. a GeoTIFF read without its companion `.aux.xml`,
  so the raster attribute table is missing), instead of a cryptic
  "subscript out of bounds" (#190, @SAY-5);
* `convertUnwantedLCC()` again returns the `newPossLCC` column it returned prior to
  1.2.0.9004 (the assigned land-cover class itself, i.e. `ecoregionGroup` without its
  ecoregion prefix). Callers use it to write the replacement classes back into the LCC
  raster — `Biomass_borealDataPrep` does so behind an `is.null()` guard, which silently
  stopped firing when the column disappeared, leaving `rstLCCAdj` (and hence
  `ecoregionMap`) still showing the replaced classes;
* `dropTerm` now can deal with random effects better (#105);
* `prepRawBiomassMap` - needed `overwrite = TRUE` for cases where download was corrupt;
* `prepRawBiomassMap` needs `httr2` package as remote site is failing with `download.file`;
* don't delete `CA_forest_VLCE2` raster in `prepInputs_NTEMS_LCC_FAO()` (#110);
* corrected some BC forestry tree species entries;
* minor bug fixes to `prepInputsFireYear` pertaining to file structure of NFDB data;
* `stats_from_counts()` (and so `calc_raster_stats()`) no longer overflows on
  national-scale pixel counts. `value` and `count` are both integer, so at 30 m
  resolution `value * count` exceeds the 2^31^ - 1 integer limit and the mean came back
  `NA` -- one ecozone of SCANFI holds ~1.5e9 pixels, and the Boreal Shield mean was
  silently lost this way. `sum(count)` and `cumsum(count)` in `weighted_quantile()`
  were within a factor of two of overflowing as well. The arithmetic now promotes to
  double. Only visible at national extent: the same code on kNN at 250 m is ~70x below
  the threshold;
* `prep_polygons()` no longer assumes the geometry column is named `geometry`. It is,
  for a shapefile, but a GeoPackage names it `geom`, and the hard-coded
  `st_union(geometry)` failed with `object 'geometry' not found`. The dissolve now goes
  through `summarise()`'s own `do_union`, so the column name is irrelevant;
* `prep_polygons()` now accepts the `SpatVector` its documentation promises; previously
  only `sf` input worked, since `sf::st_transform()` was called on the object directly;
* `plot_raster_stats()`'s statistics panel no longer runs labels together. Label and
  value were pasted into one centred string, so a wide value collided with its
  neighbour (`Min : 0Mean : 60.68Max : 395`); they are now drawn as right- and
  left-aligned columns about a shared centre;
* `plot_raster_stats()` no longer requires `purrr`, which it has not used for some time
  but still gated on via `stopifnot(requireNamespace(...))`;

# LandR 1.1.5

* use INT2U instead of INT1U when writing rasters in `.overlay()` to avoid warning with larger values;

# LandR 1.1.4

* fix bug in `vegTypeMapGenerator()` when `mixedType = 1`;
* allow `mixedType = 0` in `vegTypeMapGenerator()`;

# LandR 1.1.3

* fixed an assertion;

# LandR 1.1.2

* delete NTEMS file (24 GB) after use in `prepInputs_NTEMS_LCC_FAO()`;
* update Quebec PSP column in `speciesEquivalencies_CA`;

# LandR 1.1.1

* more conversion from `raster` to `terra` throughout;
* remove `gdalUtilities` dependency;
* fix bug in `prepInputLCC`: `orig` argument no longer accepted by `terra::compareGeom`;
* fix bug in `calcSeverityB`: output table was missing the proportion of B killed;
* new functions used to estimate maximum biomass (`maxB`) and species establishment probabilities (`SEP`);
* new function to update the `speciesEcoregion` table (brought over from `Biomass_speciesParameters` module), using estimated `inflationFactor` and `mANPPproportion` to adjust `maxB` and `maxANPP`, respectively; 
* new functions to simulate disturbances - `FireDisturbance` and `FireDisturbancePM` pulled from;
`Biomass_regeneration` and `Biomass_regenerationPM` modules, respectively;
* `overlayLCCs()` now works correctly with `terra` (#99);
* fixed partial argument match warnings (#100);
* new function `standAgeMapGenerator()` to produce `standAgeMap` from `cohortData`;
* new functions `prepInputs_NTEMS_Nonforest()` and `prepInputs_NTEMS_LCC_FAO()`;
* add new assertions: `assertSpeciesTable()` and `assertSpeciesTableRaw()`;

# LandR 1.1.0

* move LandWeb-specific functions to `LandWebUtils` package (#86)
* `Colors` -- new function to help with `terra` migration.
* `prepInputsFireYear` can now handle `rasterToMatch` that is a `SpatRaster`
* drop support for R < 4.2 due to change in dependency package `MuMIn` (Sept 2022)
* new functions to download a set of default biogeoclimatic variables used across
several modules, and subset data layers from different time periods according to 
a year. These functions may only be here temporarily.
* new `raster`/`terra` utility and plotting functions
* new functions related to producing default permafrost input data for modules.
  These functions may only be here temporarily.

# LandR 1.0.9

* new function: `nonForestedPixels` used to detect pixels without species cover or a non-forested land-cover class;
* new function: `prepRawBiomassMap` used to create `rawBiomassMap`;
* new function: `prepRasterToMatch` used to create `RasterToMatch` and `RasterToMatchLarge`;

# LandR 1.0.8

* drop support for R 3.6
* `prepInputsStandAgeMap`can now accept `firePerimeters` layer, avoiding inner download if layer is present;
* new assertion (`assertSppVectors`) to check that species match between vectors (e.g. tables, colours and species list vector);
* new function `sppHarmonize` that deals with the 3 potential ways for a user to input the `sim$sppEquiv`, `P(sim)$sppEquivCol`, and `sim$sppNameVector`;
* update Eliot's email address;
* new functions: `speciesInStudyArea` and `species;
* remove undeclared dependency package `Require`;
* age imputation in `makeAndCleanInitialCohortData` can now be turned off;
* fix bug in `LANDISDisp()`: skip dispersal when `src` or `rcv` data.tables are empty;
* initial cohort biomass (calculated in `.initiateNewCohorts`) can now be a fixed integer or, as before, calculated using the LANDIS-II Biomass Succession Extension v3.2.1 approach (Scheller & Miranda 2015);
* fix bug in `prepInputsFireYear()`;
* drop single-level factor terms in `statsModel`;

# LandR 1.0.7

* Several changes to accommodate the tracking and optional removal of pixels;
data suffered data imputation in `Biomass_borealDataPrep`;
* `assertColumns` gives better message for which columns are incorrect/missing;
* `minRelativeBDefaults` is now a function so they are more easily accessible ;
* `statsModel` was pulling along with it all the data, 5x. Now it does not. The Caching of this should be fast and small now.
 
# LandR 1.0.5

* Support for refitting `modelBiomass` (see `Biomass_boreaDataPrep`) with scaled data or different optimizer;
* Changes to `loadkNNSpeciesLayers` and `prepSpeciesLayers_KNN` prevent issues when default URL is down, or working offline (but layers are present locally);
* Several changes to accommodate LCC 2010 dataset;
* New columns to `sppEquivalenciesCA` (`PSP`, `BC_Forestry` and `FuelClass`);
* Lowered values of dummy `rawBiomassMap`;
* passing `fireURL = NULL` to `prepStandAgeMap()` bypasses age imputation;

# LandR 1.0.4

* New assertion for validation data;
* New function `sppEquivCheck`;
* `loadKNNSpeciesLayers` can accept a `sppEquiv` table with one column
* Improved documentation for `speciesEquivalencies_CA` data;

# LandR 1.0.3

* new function to calculate fire severity as biomass loss;
* bug fixes and improvements to to `speciesTableUpdate`;

# LandR 1.0.2

* Fixes and further speed improvements to seed dispersal functions and general code cleaning;

# LandR 1.0.1

* Complete rewrite of `LANDISDisp` now (back to) native R. It is about 15x faster than the Rcpp implementation, and much simpler, with about 30% of the number of lines of code. It was inspired by the "spiral" approach as was used in the Rcpp in the pre-1.0.0 version of `LandR`, but much more efficiently as it is now correctly identifies *every* pixel outward from a centre pixel using `raster::focalWeight`, with the maximum of the `seeddispersal_max` across all species. RAM use appears under control, even for large problems (tested on 50M pixel Raster with 8M potential Source pixels and 500,000 Receiving pixels, with a peak additional RAM of 3 GB during `LANDISDisp`);
`LANDISDisp` now accommodates sub-`cellSize` dispersal distances, using the original Ward Dispersal equation. Previously, the sub-pixel dispersal was treated as if it was starting from the centre of the pixel. So, if less than a full pixel, then very little horizontal transfer. This has the effect that there will be a large increase in horizontal transfer for the species that have small `seeddistance_max` (i.e., less than cell size);
* add new function `prepSpeciesLayers_ONFRI`;

## Bugfixes
* `LANDISDisp` did not correctly handle `speciesCode` when it is a factor. This is a common possibility. It now handles these correctly.

# LandR 0.0.5

* new function `updateSpeciesTable` (moved from `LandWebUtils`) to allow user to update species parameters by passing a named list.

# LandR 0.0.4

* `assignLightProb` now allows interpolating germination probability between species shade tolerance levels for any given stand shade value. This allows for for decimal values in species shade tolerance traits and greater fine tuning of shade-related germination probabilities.

# LandR 0.0.3

* rounding of age classes and biomass now occurs only inside `makeCohortDataFiles`, as it is the last thing to do before making `cohortData`
* rounding of age and biomass has been taken out of 3 other functions -- `.createCohortData`, `makeAndCleanInitialCohortData` , and a hard coded bit in `Biomass_borealDataPrep`
* updated source for kNN databases (now using v1, instead of v0) - this involved changes in the URLs and how the data is downloaded;
* `statsModel` function has a new argument to improve caching with `reproducible::Cache` (i.e. not used internally);
* function arguments that where previously called `time`, are now called `currentTime` - these changes are matched in `LandR` Biomass modules
* minor code clean-ups/bugfixes and improved clarity

# LandR 0.0.2

* Bug fixes in imports (DESCRIPTION)
* new function `overlayLCCs` which will help with overlaying more than one land cover classification raster.
* major revisions to `convertUnwantedLCC` to accommodate more cases and eliminate redundant arguments.
* `prepInputsLCC()` now works with "LCC10"
