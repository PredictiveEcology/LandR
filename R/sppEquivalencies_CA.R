#' Species equivalencies table (`sppEquiv`) and `sppEquivalencies_CA`
#'
#' @description
#' `sppEquiv` is the species name equivalencies table that `LandR` functions and
#' `LandR` `SpaDES` modules (e.g., `Biomass_borealDataPrep`, `Biomass_core`) use to
#' translate species names between data sources and naming conventions (e.g., kNN, SCANFI
#' and NTEMS species layers, LANDIS-II trait tables, CASFRI, common names).
#' Each row is one taxon (a species, variety, hybrid or genus); each column is one
#' naming convention or attribute.
#'
#' `sppEquivalencies_CA` is the default `sppEquiv` shipped with `LandR`, covering Canadian
#' tree species. It is used whenever a function's `sppEquiv` argument is not supplied, and
#' is the template for building a custom `sppEquiv`: a `data.table` with the same column
#' names (e.g., a subset of its rows) will work.
#'
#' @section Naming conventions:
#' The same species (balsam fir) in each species name column:
#'
#' | Column             | Format                     | Example              |
#' |--------------------|----------------------------|----------------------|
#' | `LandR`            | `Genu_spe`                 | `Abie_bal`           |
#' | `KNN`, `Boreal`    | `Genu_Spe`                 | `Abie_Bal`           |
#' | `NFI`, `SCANFI`    | `GENU_SPE`                 | `ABIE_BAL`           |
#' | `LANDIS_traits`    | `GENU.SPE`                 | `ABIE.BAL`           |
#' | `LANDIS_test`      | `genuspec`                 | `abiebals`           |
#' | `CASFRI`           | `Genu spec`                | `Abie bals`          |
#' | `EN_generic_short` | short common name          | `Bal fir`            |
#' | `Leading`          | common name + `" leading"` | `Balsam fir leading` |
#'
#' Varieties and hybrids add a third part (e.g., `Pice_eng_gla`, `PINU_CON_LAT`).
#' Genus-level entries use `spp` (e.g., `LandR` `Pice_spp`, `KNN` `Pice_Spp`, `NFI` `PICE_SP`).
#'
#' @section Structure of the table:
#' - A blank string (`""`) means the source has no equivalent for that taxon.
#'   Most columns are blank for most rows: only 115 rows have a `KNN` name, 37 a `SCANFI`
#'   name, 40 a `LANDIS_traits` name, and 42 a `Boreal` name.
#' - Several rows can share one name in a column. Varieties map onto their species (e.g.,
#'   the three *Pseudotsuga menziesii* rows all have `LandR` `Pseu_men`), and `Boreal`
#'   `Pinu_Con` covers both lodgepole pine varieties. Choosing such a column as
#'   `sppEquivCol` therefore merges those taxa into one species.
#' - `sppEquivCol` (in modules, usually `P(sim)$sppEquivCol`) names the column whose names
#'   are the species used in the simulation. Rows with a blank or `NA` in that column are
#'   not simulated. [sppHarmonize()] defaults it to `"Boreal"` when it cannot be determined
#'   from `sppNameVector`.
#' - To simulate a different set of species, or to merge species, subset the rows of
#'   `sppEquivalencies_CA`, or edit (or add) a naming column, and pass the result as
#'   `sppEquiv` with that column as `sppEquivCol`.
#'
#' @section Helpers:
#' - [equivalentName()]: translate names from one column to another; it finds the column
#'   the input names come from. [equivalentNameColumn()] returns that column's name.
#' - [sppEquivCheck()]: turn a character vector of names into a one-column `sppEquiv`, and
#'   if it lacks any `ensureColumns`, join it with `sppEquivalencies_CA` to get all columns.
#' - [sppHarmonize()]: make `sppEquiv`, `sppEquivCol`, `sppNameVector` and `sppColorVect`
#'   consistent with each other, subsetting `sppEquiv` to the species being simulated.
#' - [speciesInStudyArea()]: the species present in a study area, and the rows of
#'   `sppEquivalencies_CA` for them (element `sppEquiv`); see also the option
#'   `LandR.mergeHybridSpruce` in [LandROptions()].
#' - [sppColors()]: a named colour vector for the species in `sppEquivCol`, using `colorHex`
#'   when every row has one.
#' - [assertSppVectors()], [assertSpeciesPlotLabels()]: check `sppEquiv` against
#'   `sppNameVector` and `sppColorVect`, and that species plot labels are unique.
#'
#' @section Columns used by `LandR` functions:
#' - `LandR`: [speciesInStudyArea()] matches rows on it; default `sppEquivCol` of
#'   [prepSpeciesTable()].
#' - `LANDIS_traits`: [prepSpeciesTable()], [speciesTableUpdate()] and
#'   [speciesInStudyArea()], to find species traits.
#' - `KNN`: [loadkNNSpeciesLayers()] (argument `knnNamesCol`).
#' - `SCANFI`: [loadSCANFISpeciesLayers()] (argument `SCANFINamesCol`).
#' - `NTEMS_Species_Code`: [prepInputs_NTEMS_DominantSpecies()] and [speciesInStudyArea()].
#' - `CASFRI`: [loadCASFRI()].
#' - `EN_generic_short`, `Leading`: plot labels; see [assertSpeciesPlotLabels()].
#' - `Type`, `Broadleaf`: conifer vs. deciduous, e.g., in [vegTypeMapGenerator()] with
#'   `mixedType = 2`, and in [partitionBiomass()].
#' - `colorHex`: [sppColors()].
#' - `FuelClass`: fuel classes in `fireSense` modules.
#'
#' @format A `data.table` with 204 rows and 30 columns (all `character` unless stated):
#' \describe{
#'   \item{LANDIS_test}{species names from the LANDIS-II test parameter table;
#'         source: <https://raw.githubusercontent.com/LANDIS-II-Foundation/Extensions-Succession/master/biomass-succession-archive/trunk/tests/v6.0-2.0/species.txt>}
#'   \item{LANDIS_traits}{species names from the LANDIS-II traits parameter table;
#'         source: <https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv>}
#'   \item{LandR}{species names used by `LandR`; the common key for matching rows}
#'   \item{KNN}{species names from the CFS kNN species layers; source:
#'         <https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/>}
#'   \item{CASFRI}{species names from the Common Attribute Schema for Forest Resource
#'         Inventories (CASFRI); source: <https://github.com/CASFRI/CASFRI>}
#'   \item{Latin_full}{full Latin names, following The Plant List}
#'   \item{EN_generic_short}{short English common names, used as plot labels}
#'   \item{EN_generic_full}{full English common names}
#'   \item{Leading}{English common names followed by `" leading"`, used as leading-species labels}
#'   \item{Notes}{additional notes}
#'   \item{Boreal}{`KNN`-style names for a default set of 39 species (42 rows), mostly
#'         from the boreal forest; the default `sppEquivCol` of [sppHarmonize()] and
#'         [loadkNNSpeciesLayers()]}
#'   \item{Broadleaf}{`logical`: whether the species is a broadleaf}
#'   \item{Type}{`"Deciduous"` or `"Conifer"`}
#'   \item{PSP}{species names used by the `ianmseddy/PSP_Clean` module}
#'   \item{ApproxFBP}{approximate Canadian Forest Fire Behaviour Prediction (FBP) System fuel
#'         type (`"C2"`, `"C3"`, `"C5"` or `"C7"`), for 7 pine, spruce and Douglas-fir rows}
#'   \item{FuelClass}{fuel class used by `PredictiveEcology/fireSense` modules; one of
#'         `"BlkSprc"`, `"CedrMplOther"`, `"DgFrPoPine"`, `"LdgJkPine"`, `"PopBrch"`,
#'         `"RdWhPine"` or `"SprcFrLrch"`}
#'   \item{BC_forestry}{species codes used by the Government of British Columbia}
#'   \item{AB_forestry}{species codes used by the Government of Alberta}
#'   \item{SK_forestry}{species codes used by the Government of Saskatchewan}
#'   \item{MB_forestry}{species codes used by the Government of Manitoba}
#'   \item{ON_forestry}{species codes used by the Government of Ontario}
#'   \item{QCPSP}{species codes used in Québec permanent sample plot (PSP) data}
#'   \item{NB_forestry}{species codes used by the Government of New Brunswick}
#'   \item{NFI}{species codes used by the National Forest Inventory}
#'   \item{CanfiCode}{`integer`: species codes used by Canada's Forest Inventory (CanFI)}
#'   \item{CanfiNote}{notes on `CanfiCode`, e.g., generic genus codes}
#'   \item{NTEMS_Species_Code}{`integer`: species codes in the National Terrestrial Ecosystem
#'         Monitoring System (NTEMS) dominant species layer}
#'   \item{CBM_speciesID}{`integer`: species IDs used by the CFS Carbon Budget Model (CBM)}
#'   \item{SCANFI}{species names used by the Spatialized Canadian National Forest Inventory
#'         (SCANFI) species layers}
#'   \item{colorHex}{hexadecimal colour codes for plotting, used by [sppColors()]}
#' }
#'
#' @source `data-raw/sppEquivalencies_CA.csv`, built by `data-raw/sppEquivalencies_CA.R`
#'   in the package source.
#'
#' @aliases sppEquiv
#' @examples
#' library(data.table)
#' sppEquiv <- LandR::sppEquivalencies_CA
#'
#' ## translate names between conventions
#' equivalentName(c("Abie_Bal", "Pice_Mar"), sppEquiv, column = "LandR")
#' equivalentName(c("Abie_Bal", "Pice_Mar"), sppEquiv, column = "LANDIS_traits")
#' equivalentNameColumn(c("ABIE_BAL", "PICE_MAR"), sppEquiv)
#'
#' ## a custom sppEquiv with three species, using the "Boreal" names
#' mySppEquiv <- sppEquiv[Boreal %in% c("Abie_Bal", "Pice_Gla", "Popu_Tre")]
#' sppColors(mySppEquiv, sppEquivCol = "Boreal", newVals = "Mixed")
#'
#' ## expand a vector of names into a full sppEquiv
#' sppEquivCheck(c("Abie_bal", "Pice_gla"), sppEquivCol = "LandR", ensureColumns = "KNN")
"sppEquivalencies_CA"
