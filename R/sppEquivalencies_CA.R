#' Table of species name equivalencies for Canadian trees
#'
#' A table containing the different species names used across
#' diferent sources (e.g., LANDIS-II test parameter files and trait tables,
#' the LandR standard, kNN species biomass layers, ...).
#' Each column refers to a diferent source or species naming approach.
#' Presently only containing Canadian native tree species, with name equivalencies
#' coming from:
#' \itemize{
#'   \item LANDIS-II test parameter tables (column \code{LANDIStest_names};
#'         source: \url{https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/species.txt});
#'   \item LANDIS-II Canada-wide trait table (column \code{LANDIStraits_names};
#'         source: \url{https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv});
#'   \item LandR family of SpaDES modules (column \code{LandR_names};
#'   \item CFS kNN species biomass layers (column \code{KNN_names};
#'         source: \url{http://tree.pfc.forestry.ca/kNN-Species.tar});
#'   \item Canadian Common Attribute Schema for Forest Resource Inventories (column \code{CASFRI_names};
#'         source \url{http://www.borealbirds.ca/files/CAS_Document_Final_Mar_2010_ALL_APPENDICES.pdf}).
#' }
#' Remaining columns have been filled with some other useful ways to name species (e.g., for plotting).
#'
#' This table is currently used as the default equivalencies table in LandR SpaDES modules,
#' but can also serve as a template to costumise species names equivalencies by the user.
#'
#' @format A \code{data.frame} with 271 rows and 10 variables:
#' \describe{
#'   \item{LANDIS_test}{species names from LANDIS-II test parameter table}
#'   \item{LANDIS_traits}{species names from LANDIS-II traits parameter table}
#'   \item{LandR}{species names from LandR modules}
#'   \item{KNN}{species names from Knn datasets}
#'   \item{CASFRI}{species names from CASFRI database}
#'   \item{Latin_full}{accepted species latin names as in \url{http://theplantlist.org}}
#'   \item{EN_generic_short}{Short version of species' common names in English}
#'   \item{EN_generic_full}{Full species common names in English}
#'   \item{Leading}{Simple common english names used for leading species}
#'   \item{Notes}{additional notes and information}
#' }
#'
"sppEquivalencies_CA"
