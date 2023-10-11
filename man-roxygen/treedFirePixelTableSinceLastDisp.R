#' @param treedFirePixelTableSinceLastDisp data.table with 3 columns: `pixelIndex`,
#'  `pixelGroup`, and `burnTime`. Each row represents a forested pixel that was
#'  burned up to and including this year, since last dispersal event, with its
#'  corresponding `pixelGroup` and time it occurred Pixel group IDs correspond to
#'  the last year's `pixelGroupMap` and not necessarily the pixelGroupMap of
#'  the `burnTime` year.
