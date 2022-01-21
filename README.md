# LandR

<!-- badges: start -->
[![R build status](https://github.com/PredictiveEcology/LandR/workflows/R-CMD-check/badge.svg)](https://github.com/PredictiveEcology/LandR/actions)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LandR)](https://cran.r-project.org/package=LandR)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/LandR)](https://cran.r-project.org/package=LandR)
[![Codecov test coverage](https://codecov.io/gh/PredictiveEcology/LandR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PredictiveEcology/LandR?branch=master)
[![Gitter](https://badges.gitter.im/PredictiveEcology/LandR_Biomass.svg)](https://gitter.im/PredictiveEcology/LandR_Biomass?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
<!-- badges: end -->

**Landscape Ecosystem Modelling in R**

Utilities for 'LandR' suite of landscape simulation models.
These models simulate forest vegetation dynamics based on LANDIS-II, and incorporate fire and insect disturbance, as well as other important ecological processes.
Models are implemented as `SpaDES` modules.

## Installation

### Current release

[![R build status](https://github.com/PredictiveEcology/LandR/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/PredictiveEcology/LandR/actions)
[![Codecov test coverage](https://codecov.io/gh/PredictiveEcology/LandR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PredictiveEcology/LandR?branch=master)

**Install from CRAN:**

```r
#install.packages("LandR") ## not yet on CRAN
```

**Install from GitHub:**
    
```r
#install.packages("devtools")
library("devtools")
install_github("PredictiveEcology/LandR", dependencies = TRUE) 
```

### Development version

[![R build status](https://github.com/PredictiveEcology/LandR/workflows/R-CMD-check/badge.svg?branch=development)](https://github.com/PredictiveEcology/LandR/actions)
[![Codecov test coverage](https://codecov.io/gh/PredictiveEcology/LandR/branch/development/graph/badge.svg)](https://app.codecov.io/gh/PredictiveEcology/LandR?branch=development)

**Install from GitHub:**

```r
#install.packages("devtools")
library("devtools")
install_github("PredictiveEcology/LandR", ref = "development", dependencies = TRUE) 
```

## Getting help

- Discussion: https://gitter.im/PredictiveEcology/LandR_Biomass
- Bug reports: https://github.com/PredictiveEcology/LandR/issues
