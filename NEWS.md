# LandR (development version)

* **New `prepInputs_SCANFI_structure()`**: fetches SCANFI's canopy height and canopy closure
  layers (V2, 1985-2025 in 5-year steps), the two structural attributes published alongside the
  biomass layer `prepRawBiomassMap()` already serves. They say how much structure a pixel carries
  independent of which species carry it, so they can be used as controls when comparing stands:
  at equal height and closure a biomass difference is composition, not site quality.

* **Two leading/mixedwood thresholds replace three options.** `LandR.mixedwoodProp` (0.75) is the
  GROUP threshold -- all conifers, or all broadleaves -- and is the definition the national
  products use (NTEMS/EOSD, NFI photo plots: 75% of total basal area or volume).
  `LandR.leadingSpeciesProp` is the SINGLE-SPECIES threshold, and takes the mixedwood value unless
  set, so `options(LandR.leadingSpeciesProp = 0.51)` gives "just a majority" without moving the
  mixedwood definition. Read them with the new `mixedwoodProp()` and `leadingSpeciesProp()`.
  `NTEMS.mixedwoodProp`, `LandR.vegLeadingProportion` and `LandR.lccLeadingProportion` are gone.
  The number itself is written down only in `LandROptions()`; no function carries its own default.
  **This changes results:** the single-species threshold was 0.8 and is now 0.75, everywhere.
* **`mixedType = 2` now sums the broadleaf group**, which is what it always claimed to do. It
  tested each deciduous species separately, so three broadleaf species at 15% each -- 45%
  broadleaf, mixedwood by the definition -- was called pure conifer. Deciduous conifers stay
  conifers: `Larix` is `Type == "Conifer"` in `sppEquivalencies_CA`, so tamarack never makes a
  stand mixedwood. Calling `vegTypeMapGenerator()` or `vegTypeGenerator()` with `mixedType = 2`
  and a `LandR.leadingSpeciesProp` that differs from `LandR.mixedwoodProp` now warns, because
  `mixedType = 2` asks the mixedwood question and uses the mixedwood threshold.
* **`subsetDT()`'s default subsample is 500, was 50** (new `subsetDataSize()`, option
  `LandR.subsetDataSize`). 50 was chosen when these fits were expensive; it was small enough that
  repeated runs of the same simulation gave visibly different `maxB` -- a median coefficient of
  variation of 10% across ecoregion x species, up to 62%, on a 60 km boreal test window. The
  `LandR` modules take their `subsetData*Model` defaults from `subsetDataSize()`, so the number
  lives in one place.
* SCANFI files are now fetched from the PredictiveEcology arbutus mirror by default. LandR
  addresses SCANFI v2 by Google Drive id, and some of those ids 404 for anonymous users -- the
  2020 land cover and 2020 stand age among them, which stopped `Biomass_borealDataPrep`'s default
  SCANFI path. LandR now ships the mirror manifest and, when loaded, sets
  `options(reproducible.urlRemap = scanfiUrlRemap())` -- only if no remap is set and
  `LandR.scanfiMirror` is `TRUE` (the default). Species folders are remapped too, so listing them
  needs no Google login. A remap you set yourself is never replaced; `scanfiUrlRemap()` is
  exported so it can be combined with one.

* new `prepInputs_CWIM()` builds a wetland *site* layer from the Canadian Wetland Inventory Map
  v3A (10 m, national, public cloud-optimised GeoTIFF), reading only the study window. SCANFI's
  land cover has no wetland classes, so without it a SCANFI-based map cannot tell treed wetland
  from upland forest. Bog, fen, marsh and swamp count as wet; shallow water and NoData do not.
  A target cell is wet when at least `wetThreshold` (0.5) of it is.
* new `wetlandToLCC()` adds the NTEMS wetland codes to a land-cover map from such a layer: wet
  and treed (210, 220, 230, and 240) becomes 81, wet otherwise becomes 80; water and existing
  wetland codes are left alone.

* `prepInputs_NTEMS_LCC_FAO()` and `prepInputs_SCANFI_LCC_FAO()` now decide *forest land* --
  ground that grows trees, whether or not it carries any in the year being prepared -- from
  the new `forestLandFrom` argument, and share one implementation of the rule (#221).
  Previously both used the 2019 FAO layer's code 2 alone, i.e. "an opening in 2019", so a
  stand that was open in the year being prepared but had grown back by 2019 was code 1 and
  was left as shrubland, dropping it from the simulated forest. Now:
    - `"fao"` uses FAO codes 1 and 2, from `faoYear` (2022 by default, was fixed at 2019);
    - `"lccYears"` calls a pixel forest land if it is treed in any of `forestLandYears`,
      which also sees openings whose disturbance predates the 1984 start of the fire and
      harvest record;
    - `"both"` (default) takes the union. Each scanned year is one more layer to read.
  New `forestLandMask()` and `prepInputs_FAO_forest()` are exported. `convertibleClasses`
  controls which classes may be relabelled; the default, every non-treed class, is
  unchanged behaviour. The NTEMS year range is now 1984-2022: 2023 was accepted although
  NFIS publishes no 2023 land cover. The SCANFI path also gains the fast `terra::ifel`
  implementation, which the NTEMS path already had.
* new `LandROptions()`, which lists the `LandR` options and their defaults, following
  `reproducible::reproducibleOptions()` and `SpaDES.core::spadesOptions()`. `?LandROptions`
  (or `?opts.LandR`) documents each one, and `.onLoad()` now sets the options from it instead
  of from its own inline list. `NTEMS.mixedwoodProp` is a full member with a `NULL` default,
  so it is documented without being set and the
  `getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", <default>))`
  fallthrough still reaches the inner default. The package-level help now points at
  `LandROptions()` rather than repeating a two-option list that said `LandR.assertions`
  defaults to `FALSE`, when `.onLoad()` has always set it to `TRUE`.

* `speciesInStudyArea()` also returns `sppEquiv`: the rows of `sppEquivalencies_CA` for the
  species in the study area, without `_Spp` genus entries, only species with LANDIS traits,
  and with the hybrid white x Engelmann spruce (`Pice_eng_gla`) merged into Engelmann spruce
  (`Pice_eng`). This is the table fireSense modules built for themselves. The new argument
  `mergeHybridSpruce` (default `getOption("LandR.mergeHybridSpruce", "engelmann")`) merges it
  into white spruce (`"white"`, `Pice_gla`) instead, or leaves it as its own species (`NA`).
  Only the hybrid being on the raster triggers the merge, and its rows take the target's
  `LandR`, `LANDIS_traits` and `sppEquivCol` names. Rows are matched on `LandR` whatever
  naming the raster uses (SCANFI/NFI `PICE_ENG_GLA` or KNN `Pice_Eng_Gla`), so the table has
  the same rows for any `sppEquivCol`. This merged into `development` at 1.2.0.9020, the
  version already there, so **1.2.0.9021 is the first version a caller can require** for it:
  a `reqdPkgs` floor of `>= 1.2.0.9020` is also met by a 1.2.0.9020 from before the merge,
  which returns no `sppEquiv` and fails at run time instead of at install time.
* `?sppEquiv` (an alias of `?sppEquivalencies_CA`) now describes the `sppEquiv` table in one
  place: its naming conventions, how rows and `sppEquivCol` work, the helpers that use it,
  and which columns `LandR` functions read. The column list now matches the data (30 columns,
  not 27; `*_forestry` names; `SK_forestry`, `ON_forestry` and `NB_forestry` added), and the
  `sppEquiv`/`sppEquivCol` argument docs link to it. The documented `SCANFINamesCol` default of
  `loadSCANFISpeciesLayers()` is corrected to `"SCANFI"`.
* `speciesInStudyArea()` no longer stops with "object 'bb' not found" when `speciesPresentRas`
  is supplied, and uses a supplied `url` instead of ignoring it.
* `speciesTableUpdate()` no longer fails when `sppEquiv` is `NULL`. It built its default from
  `data.table(utils::data("sppEquivalencies_CA", ...))`, which holds the *name* of the dataset
  rather than the dataset, so the call died in `data.table` with "Column or expression 1 of
  'by' ... is type 'list'". It now `get()`s the table, as `prepSpeciesTable()` does.
* `sppColors()`: the test for whether `sppEquiv` has enough distinct `colorHex` values read
  `length(unique(sppEquiv[[sppEquivCol]] <= length(unique(sppEquiv$colorHex))))`, which
  compares species names to a number and takes the length of the result (1 or 2, both
  truthy), so it always passed. Two species sharing one `colorHex` were both given that
  colour instead of falling back to the palette. Also `length(newVals == 1)` is now
  `length(newVals) == 1`.
* `sppEquivalencies_CA`: the `KNN` column was shifted up by one row across the `Ulmus` block,
  so *U. pumila* carried `Ulmu_Rub`, *U. rubra* carried `Ulmu_Spp` and *Ulmus* spp. carried
  `Ulmu_Tho`. `equivalentName("Ulmu_Tho", column = "LandR")` returned the elm genus and
  `"Ulmu_Rub"` returned Siberian elm. Each name now sits on its own species.
* `sppEquivalencies_CA`: rock elm (`ULMU_THO`) and pagoda dogwood (`CORN_ALT`) now have the
  `LandR` names `Ulmu_tho` and `Corn_alt`. Both were blank, and `LandR` is the column rows
  are keyed on, so neither species could be matched.

* the "leading" threshold is no longer hard-coded in each function. Every site now reads a
  nested pair of options,
  `getOption("NTEMS.mixedwoodProp", getOption("LandR.<which>LeadingProportion", <default>))`,
  where `LandR.vegLeadingProportion` (0.8) serves `vegTypeMapGenerator()`, `vegTypeGenerator()`
  and `plotVTM()`, and `LandR.lccLeadingProportion` (0.75) serves `lccMapGenerator()`.
  Setting `NTEMS.mixedwoodProp` moves all of them at once; leaving it unset (the default --
  it is not set at load) leaves each on the value it has always had, so no existing result
  changes. The two inner defaults differ by history rather than by concept: both are the same
  purity threshold on a biomass-like share, and 0.75 is the NTEMS/EOSD value (Wulder & Nelson
  2003: coniferous or broadleaf at 75% or more of total basal area, mixed wood below that).
* **`loadSCANFISpeciesLayers()` and `prepSpeciesLayers_SCANFI()` now take the `*to` family
  (`to`, `cropTo`, `projectTo`, `maskTo`) as formals**, as a first step in retiring
  `rasterToMatch`/`studyArea`. Both still accept the legacy pair -- it arrives through `...`
  and is translated by a new internal `.legacyToTo()`, which implements the table documented
  in `?reproducible::postProcess`: a `rasterToMatch` on its own is `to`; a `studyArea` on its
  own crops and masks but does not reproject (unless `useSAcrs`); and when both are supplied
  the raster gives extent, resolution, projection and alignment while the polygon gives the
  mask. An explicitly passed `*to` argument always wins.
  **This changes output for callers that supplied both.** The previous shims mapped
  `to` -> `studyArea` and `projectTo` -> `rasterToMatch`, then called
  `prepInputs(to = rasterToMatch)`, so the mask was taken from the *raster* rather than from
  the study area -- the reverse of the documented behaviour. Callers that supplied only one
  of the two are unaffected.
  The other `prepSpeciesLayers_*()` functions are unchanged and still take the legacy formals.
* `prepSpeciesLayers_SCANFI()` passed `projectTo = rasterToMatch` twice to
  `loadSCANFISpeciesLayers()`. R accepts duplicate names in `...`, so this was silent
  rather than an error; the duplicate is removed. Its cache entry was also tagged
  `"KNN"`, which is now `"SCANFI"` -- the tag is what `Cache()` searches on, so SCANFI
  species layers were indistinguishable from kNN ones in the cache.
* `sppEquivalencies_CA`: coastal Douglas-fir (`PSEU_MEN_MEN`) now has `FuelClass`
  "DgFrPoPine", like the other two `Pseu_men` rows. It had "CedrMplOther", so any study
  area containing Douglas-fir got two fuel classes for `Pseu_men` and
  `fireSenseUtils::cohortsToFuelClasses()` stopped.
* `makePickellStack()` no longer leaves `terra::terraOptions(memmax)` and
  `raster::rasterOptions(maxmemory)` changed after it returns. Both are global,
  process-wide settings, so a caller that had set its own memory ceiling silently
  kept LandR's for the rest of the session -- visible in a long-lived worker that
  set `memmax = 4` and later found it at 1.

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
* `loadSCANFISpeciesLayers()` now lists the SCANFI species-layer folder via
  `reproducible::listGoogleDriveFolder()` instead of `googledrive::drive_ls()`.
  When a directory-remap manifest is set (`options(reproducible.urlRemap = ...)`,
  e.g. one built with `buckethost::makeMirrorManifest(directories = TRUE)`), the
  SCANFI files are fetched from a public mirror with **no Google
  authentication** — the immediate use case is training/workshops, where
  participants can pull SCANFI layers without a Google account or `drive_auth()`.
  Behaviour is unchanged when no manifest is set: it falls back to `drive_ls()`
  and authenticates as before. (Requires the companion `reproducible` change that
  adds `listGoogleDriveFolder()` and directory remaps.)
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
