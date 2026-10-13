## Makes fastCohortOps-expected.rds. Run it on LandR before the leadingByPixelGroup() /
## .initiateNewCohorts() / calcSeverityB() / rmMissingCohorts() speed-ups (1.2.0.9052, f3c5a3e1),
## from the package root, with the package loaded: source("tests/testthat/fixtures/make-fastCohortOps-expected.R")
library(data.table)
source("tests/testthat/helper-fastCohortOps.R")
saveRDS(fastCohortOpsRun(), "tests/testthat/fixtures/fastCohortOps-expected.rds")
