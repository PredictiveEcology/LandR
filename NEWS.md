bug reports https://github.com/PredictiveEcology/LandR/issues

version 0.0.3
=============
* rounding of age classes and biomass now occurs only inside `makeCohortDataFiles`, as it is the last thing to do before making `cohortData`
* rounding of age and biomass has been taken out of 3 other functions -- .createCohortData, makeAndCleanInitialCohortData , and a hard coded bit in Biomass_borealDataPrep 
* updated source for kNN databases (now using v1, instead of v0) - this involved changes in the URLs and how the data is downloaded;
* `statsModel` function has a new argument to improve caching with `reproducible::Cache` (i.e. not used internally);
* function arguments that where previously called `time`, are now called `currentTime` - these changes are matched in LandR Biomass modules
* minor code clean-ups/bugfixes and improved clarity

version 0.0.2
=============

* Bug fixes in imports (DESCRIPTION)
* new function `overlayLCCs` which will help with overlaying more than one land cover classification raster.
* major revisions to `convertUnwantedLCC` to accommodate more cases and eliminate redundant arguments.
* `prepInputsLCC()` now works with "LCC10"
