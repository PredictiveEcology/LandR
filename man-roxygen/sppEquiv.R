#' @param sppEquiv table with species name equivalencies between the kNN and final naming formats.
#'     See \code{data("sppEquivalencies_CA", "LandR")}.
#'     For functions that have \code{mixedType} e.g., \code{vegTypeMapGenerator},
#'     this only necessary if \code{mixedType == 2}.
#'     If not provided and \code{mixedType == 2}, will attempt to use
#'     \code{data("sppEquivalencies_CA", "LandR")}.
