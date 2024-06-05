#' Create default study areas for use with LandR modules
#'
#' This simply re-exports `SpaDES.tools::randomStudyArea`
#'
#' @inheritParams SpaDES.tools::randomStudyArea
#'
#' @export
randomStudyArea <- utils::getFromNamespace("randomStudyArea", "SpaDES.tools")
