## Makes regroup-expected.rds. Run it on LandR before the group-level re-grouping in updateCohortData(),
## the single-pass calcSiteShade() and the vector rmMissingCohorts() (1.2.0.9054, 04bfdcd8), from the
## package root, with the package loaded: source("tests/testthat/fixtures/make-regroup-expected.R")
library(data.table)
source("tests/testthat/helper-fastCohortOps.R")
source("tests/testthat/helper-regroup.R")
saveRDS(regroupRun(), "tests/testthat/fixtures/regroup-expected.rds")
