utils::globalVariables(c(
  "possLCC"
))
#' Convert Land Cover Classes (LCC) to another value in its neighbourhood
#'
#' This will search around the pixels on \code{rstLCC} that have
#' \code{classesToReplace}, and search in iteratively increasing
#' radii outwards for other Land Cover Classes than the those indicated in
#' \code{classesToReplace}.
#'.
#' @param classesToReplace Integer vector of classes that are are to be replaced,
#'     e.g., 34, 35, 36 on LCC2005, which are burned young, burned 10 year, and cities.
#' @param rstLCC LCC raster, e.g., LCC2005
#' @param nIterations the number of iterations to perform
#' @param defaultNewValue the value to assign a pixel in \code{classesToReplace} if no valid  pixel
#'  is closer after \code{nIterations}
#' @param invalidClasses classes that are not valid options
#'
#'
#'
#' @return
#' A \code{rasterLayer} with values in classesToReplace converted to adjacent values or NA.
#'
#' @author Eliot McIntire Ian Eddy
#' @export
#' @importFrom data.table as.data.table is.data.table rbindlist setnames
#' @importFrom raster raster ncell getValues
#' @importFrom SpaDES.tools spread2
convertUnwantedLCC2 <- function(classesToReplace = 34:36, rstLCC, nIterations = 6, defaultNewValue = NA,
                                invalidClasses = NA) {

  lccDT <- data.table(pixelIndex = 1:ncell(rstLCC), lcc = getValues(rstLCC))

  theUnwantedPixels <- lccDT[!is.na(lcc) & lcc %in% classesToReplace]$pixelIndex
  # lccShort <- lccDT[pixelIndex %in% theUnwantedPixels] #keep for joining
  lccOut <- data.table("pixelIndex" = integer(0), lcc = integer(0),
                       pixels = integer(0), possLCC = integer(0))
  iterations <- 1
  # remove the lines that have the code "classesToReplace"
  repeatsOnSameUnwanted <- 0

  while (length(theUnwantedPixels) > 0) {
    message("Converting unwanted LCCs: ", length(theUnwantedPixels), " pixels remaining.")
    out <- spread2(rstLCC,
                   start = theUnwantedPixels, asRaster = FALSE,
                   iterations = iterations, allowOverlap = TRUE, spreadProb = 1
    )
    out <- out[initialPixels != pixels] # rm pixels which are same as initialPixels --> these are known wrong
    iterations <- iterations + 1
    out[, possLCC := rstLCC[][pixels]]
    out[possLCC %in% c(classesToReplace, invalidClasses), possLCC := NA]
    out <- na.omit(out)
    out2 <- lccDT[out[, state := NULL],
                   allow.cartesian = TRUE,
                   on = c("pixelIndex" = "initialPixels"), nomatch = NA]

    # sanity check -- don't let an infinite loop

    if (repeatsOnSameUnwanted > nIterations) {
            message(
        "  reclassifying ", NROW(theUnwantedPixels), " pixels of class ",
        paste(unique(rstLCC[theUnwantedPixels]), collapse = ", "), " as class ", defaultNewValue,
        " because no suitable replacement was found"
      )
      pixelsToNA <- theUnwantedPixels
      theUnwantedPixels <- integer()
    }

    lccOut <- rbind(out2, lccOut)

    theUnwantedPixels <- theUnwantedPixels[!theUnwantedPixels %in% out2$pixelIndex]
    repeatsOnSameUnwanted <- repeatsOnSameUnwanted + 1
  }

  ## take random sample of available, weighted by abundance
  rowsToKeep <- lccOut[, list(keep = .resample(.I, 1)), by = c("pixelIndex")]
  lccOut <- lccOut[rowsToKeep$keep]
  rstLCC[lccOut$pixelIndex] <- lccOut$possLCC
  rstLCC[pixelsToNA] <- defaultNewValue

 return(rstLCC)

}
