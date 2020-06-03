# LandR

<!-- badges: start -->
[![R build status](https://github.com/PredictiveEcology/LandR/workflows/R-CMD-check/badge.svg)](https://github.com/PredictiveEcology/LandR/actions)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LandR)](https://cran.r-project.org/package=LandR)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/LandR)](https://cran.r-project.org/package=LandR)
<!-- badges: end -->

**Landscape Ecosystem Modelling in R**

Utilities for 'LandR' suite of landscape simulation models.
These models simulate forest vegetation dynamics based on LANDIS-II, and incorporate fire and insect disturbance, as well as other important ecological processes.
Models are implemented as `SpaDES` modules.

## Installation

### Current release

[![Build Status](https://travis-ci.org/PredictiveEcology/LandR.svg?branch=master)](https://travis-ci.org/PredictiveEcology/LandR)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/2fxqhgk6miv2fytd/branch/master?svg=true)](https://ci.appveyor.com/project/achubaty/LandR/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/PredictiveEcology/LandR/badge.svg?branch=master)](https://coveralls.io/github/PredictiveEcology/LandR?branch=master)

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

[![Build Status](https://travis-ci.org/PredictiveEcology/LandR.svg?branch=development)](https://travis-ci.org/PredictiveEcology/LandR)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/2fxqhgk6miv2fytd/branch/development?svg=true)](https://ci.appveyor.com/project/achubaty/LandR/branch/development)
[![Coverage Status](https://coveralls.io/repos/github/PredictiveEcology/LandR/badge.svg?branch=development)](https://coveralls.io/github/PredictiveEcology/LandR?branch=development)

**Install from GitHub:**

```r
#install.packages("devtools")
library("devtools")
install_github("PredictiveEcology/LandR", ref = "development", dependencies = TRUE) 
```
