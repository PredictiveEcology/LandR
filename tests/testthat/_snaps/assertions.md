# assertERGs errors informatively when ecoregionMap has no categories

    Code
      assertERGs(f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
      doAssertion = TRUE)
    Condition
      Error in `assertERGs()`:
      ! `ecoregionMap` has no 'ecoregionGroup' values. Either its raster attribute
      table is missing, it has no 'ecoregionGroup' column, or all its values are NA.
      Note that a GeoTIFF stores its categories in a companion '.aux.xml' file,
      which must be retrieved alongside the '.tif'.

# assertERGs errors informatively when ecoregionMap has no ecoregionGroup column

    Code
      assertERGs(f$ecoregionMap, f$cohortData, f$speciesEcoregion, f$minRelativeB,
      doAssertion = TRUE)
    Condition
      Error in `assertERGs()`:
      ! `ecoregionMap` has no 'ecoregionGroup' values. Either its raster attribute
      table is missing, it has no 'ecoregionGroup' column, or all its values are NA.
      Note that a GeoTIFF stores its categories in a companion '.aux.xml' file,
      which must be retrieved alongside the '.tif'.

