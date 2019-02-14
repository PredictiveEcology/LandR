#' Create default study areas for use with LandR modules
#'
#' This simply re-exports \code{SpaDES.tools::randomStudyArea}
#'
#' @export
#' @importFrom utils getFromNamespace
#' @inheritParams SpaDES.tools::randomStudyArea
randomStudyArea <- getFromNamespace("randomStudyArea", "SpaDES.tools")
