#' Create default study areas for use with LandR modules
#'
#' This simply re-exports \code{SpaDES.tools::randomStudyArea}
#'
#' @inheritParams SpaDES.tools::randomStudyArea
#'
#' @export
#' @importFrom utils getFromNamespace
randomStudyArea <- getFromNamespace("randomStudyArea", "SpaDES.tools")
