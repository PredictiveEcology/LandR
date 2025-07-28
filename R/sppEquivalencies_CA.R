#' Table of species name equivalencies for Canadian trees
#'
#' A table containing the different species names used across different sources (e.g., LANDIS-II
#' test parameter files and trait tables, the LandR standard, kNN species biomass layers, etc.).
#' Each column refers to a different source or species naming approach.
#' Presently only containing Canadian native tree species, with name equivalencies coming from:
#'
#' This table is currently used as the default equivalencies table in LandR SpaDES modules,
#' but can also serve as a template to customize species names equivalencies by the user.
#'
#' @format A `data.frame` with 204 rows and 27 variables:
#' \describe{
#'   \item{LANDIS_test}{species names from LANDIS-II test parameter table;
#'         source: <https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/species.txt>}
#'   \item{LANDIS_traits}{species names from LANDIS-II traits parameter table;
#'         source: <https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv>}
#'   \item{LandR}{species names from LandR modules}
#'   \item{KNN}{species names from CFS kNN datasets; source: <http://tree.pfc.forestry.ca/kNN-Species.tar>}
#'   \item{CASFRI}{species names from Canadian Common Attribute Schema for Forest Resource Inventories;
#'         source <http://www.borealbirds.ca/files/CAS_Document_Final_Mar_2010_ALL_APPENDICES.pdf>}
#'   \item{Latin_full}{accepted species latin names as in <http://theplantlist.org>}
#'   \item{EN_generic_short}{Short version of species' common names in English}
#'   \item{EN_generic_full}{Full species common names in English}
#'   \item{Leading}{Simple common English names used for leading species}
#'   \item{Notes}{additional notes and information}
#'   \item{Boreal}{Species present in the Boreal forests of Canada}
#'   \item{Broadleaf}{logical indicating whether the species is broad leaf}
#'   \item{Type}{Whether the species is a deciduous or conifer species}
#'   \item{PSP}{Species name from the module `ianmseddy/PSP_Clean`}
#'   \item{ApproxFBP}{Species groups roughly corresponding to Canadian Forest Fire Behavior Prediction (FBP) System}
#'   \item{FuelClass}{The fuel class used by the module PredictiveEcology/fireSense}#'
#'   \item{BC_Forestry}{Species code adopted by the Government of British Columbia}
#'   \item{AB_Forestry}{Species code adopted by the Government of Alberta}
#'   \item{MB_Forestry}{Species code adopted by the Government of Manitoba}
#'   \item{NFI}{Species code used by the National Forest Inventory}
#'   \item{QCPSP}{Species code used in Québec PSP data}
#'   \item{CanfiCode}{Species code used by Canada's Forest Inventory (CANFI)}
#'   \item{CanfiNote}{Additional notes for CANFI species data}
#'   \item{NTEMS_Species_Code}{Species code used by the National Terrestrial Ecosystem Monitoring System (NTEMS)}
#'   \item{CBM_speciesID}{Species code used by the CFS Carbon Budget Model}
#'   \item{SCANFI}{Species names from the Spatialized Canadian National Forest Inventory (SCANFI)}
#'   \item{colorHex}{hexadecimal colour codes for use with plotting}
#' }
#'
"sppEquivalencies_CA"
