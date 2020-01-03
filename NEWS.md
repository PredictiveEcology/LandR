bug reports https://github.com/PredictiveEcology/LandR/issues

version 0.0.3
=============

* updated source for kNN databases (now using v1, instead of v0) - this involved changes in the URLs and how the data is downloaded;
* `statsModel` function has a new argument to improve caching with `reproducible::Cache` (i.e. not used internally);
* function arguments that where previously called `time`, are now called `currentTime` - these changes are matched in LandR Biomass modules
* minor code clean-ups/bugfixes and improved clarity
