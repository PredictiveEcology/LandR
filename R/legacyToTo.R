## `$` on a list partial-matches, so `dots$studyArea` silently returns a `studyAreaName` that the
## caller passed through `...`, and the legacy translation below then fires on a character string.
## Legacy names must therefore be read exactly. An explicitly-NULL entry counts as absent.
.legacyDot <- function(dots, name) {
  if (name %in% names(dots)) dots[[name]] else NULL
}

#' Translate the legacy `rasterToMatch` / `studyArea` pair into the `*to` family
#'
#' `rasterToMatch` and `studyArea` are being retired in favour of
#' [reproducible::postProcessTo()]'s `to`, `cropTo`, `projectTo` and `maskTo`. The old pair is
#' ambiguous -- what it asks for depends on *which* of the two was supplied -- and that
#' ambiguity is the reason for the migration. This function is the single place that resolves
#' it, following the table documented in [reproducible::postProcess()]:
#'
#' |            | `rasterToMatch` | `studyArea` | both            |
#' |------------|-----------------|-------------|-----------------|
#' | extent     | yes             | yes         | `rasterToMatch` |
#' | resolution | yes             | no          | `rasterToMatch` |
#' | projection | yes             | no*         | `rasterToMatch` |
#' | alignment  | yes             | no          | `rasterToMatch` |
#' | mask       | no**            | yes         | `studyArea`     |
#'
#' Notes: `*` overridden by `useSAcrs`; `**` masks with `rasterToMatch`'s own `NA`s if
#' `maskWithRTM`.
#'
#' So: a `rasterToMatch` on its own defines every geometry property and is simply `to`; a
#' `studyArea` on its own crops and masks but does **not** reproject; and when both are given
#' the raster supplies the geometry while the polygon supplies the mask.
#'
#' An explicitly supplied `*to` argument always wins -- this only fills what the caller left
#' `NULL` -- so a caller that has already migrated is never overridden by a legacy argument
#' arriving through `...`.
#'
#' @param to,cropTo,projectTo,maskTo the `*to` family as the caller supplied them.
#' @param rasterToMatch,studyArea the legacy pair, usually pulled out of `...`.
#' @param useSAcrs logical. A `studyArea` on its own does not supply the CRS; set this to take
#'   the CRS from it. Ignored unless `studyArea` is the only one supplied.
#' @param maskWithRTM logical. A `rasterToMatch` on its own masks with its own `NA`s; set
#'   `FALSE` to omit masking. Ignored unless `rasterToMatch` is the only one supplied.
#'
#' @return a named list with elements `to`, `cropTo`, `projectTo` and `maskTo`, ready to pass
#'   to [reproducible::prepInputs()]. `NA` means "omit this step", per `reproducible`.
#'
#' @keywords internal
#' @rdname legacyToTo
.legacyToTo <- function(to = NULL, cropTo = NULL, projectTo = NULL, maskTo = NULL,
                        rasterToMatch = NULL, studyArea = NULL,
                        useSAcrs = FALSE, maskWithRTM = TRUE) {
  out <- list(to = to, cropTo = cropTo, projectTo = projectTo, maskTo = maskTo)

  ## nothing legacy to translate, or the caller already said everything explicitly
  if (is.null(rasterToMatch) && is.null(studyArea)) {
    return(out)
  }
  if (!any(vapply(out, is.null, logical(1)))) {
    return(out)
  }

  legacy <- if (!is.null(rasterToMatch) && !is.null(studyArea)) {
    ## the raster gives extent, resolution, projection and alignment; the polygon masks
    list(to = NULL, cropTo = rasterToMatch, projectTo = rasterToMatch, maskTo = studyArea)
  } else if (!is.null(rasterToMatch)) {
    ## on its own it defines every geometry property, so it is simply `to`
    list(to = rasterToMatch, cropTo = NULL, projectTo = NULL,
         maskTo = if (isFALSE(maskWithRTM)) NA else NULL)
  } else {
    ## on its own it crops and masks, but reprojects only when asked
    list(to = NULL, cropTo = studyArea, maskTo = studyArea,
         projectTo = if (isTRUE(useSAcrs)) studyArea else NA)
  }

  for (nm in names(out)) {
    if (is.null(out[[nm]])) {
      out[[nm]] <- legacy[[nm]]
    }
  }
  out
}
