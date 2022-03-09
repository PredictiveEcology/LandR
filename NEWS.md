bug reports https://github.com/PredictiveEcology/LandR/issues

version 1.0.7.9000
=============
* drop support for R 3.6
* update Eliot's email address
* new functions: `speciesInStudyArea` and `species
* remove undeclared dependency package `Require`
* fix bug in `LANDISDisp()`: skip dispersal when src or rcv data.tables are empty

version 1.0.7
=============
* Several changes to accommodate the tracking and optional removal of pixels
data suffered data imputation in `Biomass_borealDataPrep`
* `assertColumns` gives better message for which columns are incorrect/missing
* `minRelativeBDefaults` is now a function so they are more easily accessible 
* `statsModel` was pulling along with it all the data, 5x. Now it does not. The Caching of this should be fast and small now.
 
version 1.0.5
=============
* Support for refitting `modelBiomass` (see `Biomass_boreaDataPrep`) with scaled data or different optimizer
* Changes to `loadkNNSpeciesLayers` and `prepSpeciesLayers_KNN` prevent issues when default URL is down, or working offline (but layers are present locally)
* Several changes to accommodate LCC2010
* New columns to `sppEquivalenciesCA` (`PSP`, `BC_Forestry` and `FuelClass`)
* Lowered values of dummy `rawBiomassMap`
* passing `fireURL = NULL` to `prepStandAgeMap()` bypasses age imputation

version 1.0.4
=============
* New assertion for validation data
* New function `sppEquivCheck`
* `loadKNNSpeciesLayers` can accept a `sppEquiv` table with one column
* Improved documentation for `speciesEquivalencies_CA` data

version 1.0.3
=============
* new function to calculate fire severity as biomass loss
* bug fixes and improvements to to `speciesTableUpdate`

version 1.0.2
=============
* Fixes and further speed improvements to seed dispersal functions and general code cleaning

version 1.0.1
=============
* Complete rewrite of `LANDISDisp` now (back to) native R. It is about 15x faster than the Rcpp implementation, and much simpler, with about 30% of the number of lines of code. It was inspired by the "spiral" approach as was used in the Rcpp in the pre-1.0.0 version of `LandR`, but much more efficiently as it is now correctly identifies *every* pixel outward from a center pixel using `raster::focalWeight`, with the maximum of the `seeddispersal_max` across all species. RAM use appears under control, even for large problems (tested on 50M pixel Raster with 8M potential Source pixels and 500,000 Receiving pixels, with a peak additional RAM of 3 GB during `LANDISDisp`)
`LANDISDisp` now accommodates sub-cellSize dispersal distances, using the original Ward Dispersal equation. Previously, the sub-pixel dispersal was treated as if it was starting from the centre of the pixel. So, if less than a full pixel, then very little horizontal transfer. This has the effect that there will be a large increase in horizontal transfer for the species that have small `seeddistance_max` (i.e., less than cell size)
* add new function `prepSpeciesLayers_ONFRI`

## Bugfixes
* `LANDISDisp` did not correctly handle `speciesCode` when it is a factor. This is a common possibility. It now handles these correctly.

version 0.0.5
=============
* new function `updateSpeciesTable` (moved from `LandWebUtils`) to allow user to update species parameters by passing a named list.

version 0.0.4
=============
* `assignLightProb` now allows interpolating germination probability between species shade tolerance levels for any given stand shade value. This allows for for decimal values in species shade tolerance traits and greater fine tuning of shade-related germination probabilities.

version 0.0.3
=============
* rounding of age classes and biomass now occurs only inside `makeCohortDataFiles`, as it is the last thing to do before making `cohortData`
* rounding of age and biomass has been taken out of 3 other functions -- `.createCohortData`, `makeAndCleanInitialCohortData` , and a hard coded bit in `Biomass_borealDataPrep`
* updated source for kNN databases (now using v1, instead of v0) - this involved changes in the URLs and how the data is downloaded;
* `statsModel` function has a new argument to improve caching with `reproducible::Cache` (i.e. not used internally);
* function arguments that where previously called `time`, are now called `currentTime` - these changes are matched in `LandR` Biomass modules
* minor code clean-ups/bugfixes and improved clarity

version 0.0.2
=============
* Bug fixes in imports (DESCRIPTION)
* new function `overlayLCCs` which will help with overlaying more than one land cover classification raster.
* major revisions to `convertUnwantedLCC` to accommodate more cases and eliminate redundant arguments.
* `prepInputsLCC()` now works with "LCC10"
