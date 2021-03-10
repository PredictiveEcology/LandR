bug reports https://github.com/PredictiveEcology/LandR/issues

version 1.0.0
=============
* `LANDISDisp` now accommodates sub-cellSize dispersal distances, using the original Ward Dispersal equation. Previously, the sub-pixel dispersal was treated as if it was starting from the centre of the pixel. So, if less than a full pixel, then very little horizontal transfer. This has the effect that there will be a large increase in horizontal transfer for the species that have small `seeddistance_max` (i.e., less than cell size)
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
