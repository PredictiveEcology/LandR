#' Customize species trait table values
#'
#' @template speciesTable
#' @param params       A named list (of parameters) of named lists (by species), with species
#'                     traits overrides (e.g., `list(seeddistance_eff = list(Abie_sp = 25))`).
#'
#' @author Alex Chubaty and Ceres Barros
#' @export
updateSpeciesTable <- function(speciesTable, params) {
  ## checks:
  traits <- names(params)
  missingTraits <- traits[!traits %in% names(speciesTable)]
  if (length(missingTraits)) {
    stop("The traits: ", paste(missingTraits, collapse = ", "), "\ndo not exist in `speciesTable`")
  }

  for (trt in traits) {
    subParams <- params[[trt]]

    ## checks:
    spp <- names(subParams)
    missingSpp <- spp[!spp %in% speciesTable$species]
    if (length(missingSpp)) {
      stop(
        "The species: ", paste(missingSpp, collapse = ", "),
        "\ndo not exist in `speciesTable$species`"
      )
    }

    ## this is sub ideal to convert classes, but the only solution I found.
    ## neither class(...) <- class(..) nor as(..., class(...)) were changing integer to numeric.
    ## only as.numeric(...) worked
    if (class(unlist(subParams)) != class(speciesTable[[trt]])) {
      asFUN <- get(paste0("as.", class(unlist(subParams))))
      if (!is.function(asFUN)) {
        stop("Can't find an `as.` method to convert class to '", class(unlist(subParams)), "'")
      }

      speciesTable <- speciesTable[, (trt) := lapply(.SD, asFUN), .SDcols = trt]
    }

    for (sp in spp) {
      set(speciesTable, which(speciesTable$species == sp), trt, subParams[[sp]])
    }
  }

  return(speciesTable)
}
