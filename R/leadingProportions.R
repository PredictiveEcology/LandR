#' Leading-species and mixedwood proportion thresholds
#'
#' Two proportions decide how a stand's composition is labelled, and they answer different
#' questions. Both are read from options so a project sets them once.
#'
#' `mixedwoodProp()` is the **group** threshold: the share of the stand held by all conifers, or
#' by all broadleaves, at or above which the stand stops being mixedwood. It is the definition
#' the national products are written on -- NTEMS/EOSD (Wulder & Nelson 2003) and the NFI photo
#' plot dictionary call a stand coniferous or broadleaf at 75% or more of total basal area (NFI
#' photo plots: total tree volume), and mixed wood when neither group reaches it.
#'
#' `leadingSpeciesProp()` is the **single-species** threshold: the share one species must hold to
#' be called leading. It is a different question and a user may reasonably want a different
#' answer -- 0.51, "just a majority", is a sensible setting -- so it has its own option. When it
#' is not set it takes the mixedwood value, which is what makes the two agree by default.
#'
#' Set them with, e.g., `options(LandR.mixedwoodProp = 0.7)` or
#' `options(LandR.leadingSpeciesProp = 0.51)`. The shipped default for `LandR.mixedwoodProp` lives
#' in [LandROptions()] and is the only place the number itself is written down, so changing it
#' later changes every function at once.
#'
#' Deciduous conifers are conifers here. `Larix` is `Type == "Conifer"` in
#' [sppEquivalencies_CA], so it never counts toward the broadleaf share.
#'
#' @return A numeric proportion between 0 and 1.
#'
#' @seealso [vegTypeMapGenerator()], [lccMapGenerator()]
#'
#' @export
#' @rdname leadingProportions
#' @aliases leadingProportions
mixedwoodProp <- function() {
  getOption("LandR.mixedwoodProp", LandROptions()[["LandR.mixedwoodProp"]])
}

#' @export
#' @rdname leadingProportions
leadingSpeciesProp <- function() {
  getOption("LandR.leadingSpeciesProp", mixedwoodProp())
}

## The threshold a caller that did not supply one should use, given what it is asking.
##
## `mixedType = 2` asks the mixedwood question (conifer group vs broadleaf group), so it takes
## `mixedwoodProp()`. Every other mixedType asks about a single species and takes
## `leadingSpeciesProp()`. Because `mixedType = 2` is the default in every caller, a user who
## sets `LandR.leadingSpeciesProp` and changes nothing else would silently see no effect; that
## is worth a warning rather than a surprise.
.leadingProp <- function(vegLeadingProportion = NULL, mixedType = 2) {
  if (!is.null(vegLeadingProportion)) {
    return(vegLeadingProportion)   ## an explicit argument always wins, and is never ambiguous
  }

  if (isTRUE(mixedType == 2)) {
    set <- getOption("LandR.leadingSpeciesProp")
    if (!is.null(set) && !isTRUE(all.equal(set, mixedwoodProp()))) {
      warning(
        "'LandR.leadingSpeciesProp' is set to ", set, ", but mixedType = 2 asks whether the ",
        "stand is mixedwood (all conifers vs all broadleaves), so it uses ",
        "'LandR.mixedwoodProp' (", mixedwoodProp(), ") instead. Set 'LandR.mixedwoodProp' to ",
        "change this call, or use mixedType = 1 for a single-species leading rule."
      )
    }
    mixedwoodProp()
  } else {
    leadingSpeciesProp()
  }
}

## Is this stand mixedwood? The broadleaf share is summed over the WHOLE broadleaf group and
## compared with the threshold; the conifer share is its complement. Called once per
## pixelGroup, so `speciesProportion` and `Type` are that group's rows.
##
## Deciduous conifers are conifers: `Type` comes from `sppEquiv`, where `Larix` is "Conifer",
## so tamarack biomass counts toward the conifer share and never makes a stand mixedwood.
## `Type %in% "Deciduous"` (not `==`) so a species missing from `sppEquiv`, whose Type is NA
## after the merge, is simply not counted as broadleaf rather than poisoning the sum.
##
## Previously this tested each deciduous species SEPARATELY -- a stand was mixedwood if any ONE
## deciduous species sat inside the band. Three broadleaf species at 0.15 each is 45% broadleaf,
## which is mixedwood by the group definition, but no single species reached the band, so the
## stand was called pure conifer.
.isMixedwood <- function(speciesProportion, Type, vegLeadingProportion) {
  broadleaf <- sum(speciesProportion[Type %in% "Deciduous"], na.rm = TRUE)
  broadleaf < vegLeadingProportion & broadleaf > 1 - vegLeadingProportion
}

#' Leading species of each `pixelGroup`
#'
#' The one place where the leading species (or "Mixed") of every `pixelGroup` in a `cohortData`
#' table is worked out. [vegTypeGenerator()] and [vegTypeMapGenerator()] use it, and so can
#' a caller that only needs the result per `pixelGroup`, which avoids the join back onto the
#' (much longer) `cohortData` that [vegTypeGenerator()] does.
#'
#' With `mixedType = 2`, `B`-based leading type and a factor `speciesCode` (what `cohortData`
#' normally has), the species sums per `pixelGroup` are made in one pass on an integer key and
#' the broadleaf proportion of the (few) `pixelGroup`s with two or more deciduous species is
#' summed per group with [base::sum()], exactly as `.isMixedwood()` does. Any other case
#' (`mixedType` 0 or 1, `cover`, a character `speciesCode`, `NA` in `B`, a species with two
#' `Type`s in `sppEquiv`) takes the original, slower algorithm. Both give identical results.
#'
#' @param x A `cohortData` object, with `B` (or `cover`) and `speciesCode`.
#'
#' @template vegLeadingProportion
#'
#' @param mixedType An integer defining whether mixed stands are:
#'                  not differentiated (0), any kind of species admixture (1),
#'                  or deciduous mixed with conifer (2; default).
#'
#' @template sppEquiv
#'
#' @template sppEquivCol
#'
#' @param pixelGroupColName Name of the column in `x` identifying the `pixelGroup`.
#'
#' @param leadingBasedOn `"B"` (default if `x` has a `B` column) or `"cover"`.
#'
#' @return A `data.table` with one row per `pixelGroup`: the `pixelGroupColName` column,
#'   `leading` (a factor, or `speciesCode` when `mixedType = 0`) and, for `mixedType` 1 and 2,
#'   `mixed`.
#'
#' @author Eliot McIntire, Ceres Barros, Alex Chubaty
#' @export
#' @rdname leadingByPixelGroup
#'
#' @examples
#' library(data.table)
#' x <- data.table(
#'   pixelGroup = rep(1:2, each = 2), B = c(100, 200, 20, 400),
#'   speciesCode = factor(rep(c("Pice_gla", "Popu_tre"), 2))
#' )
#' sppEquiv <- data.table(LandR = c("Pice_gla", "Popu_tre"), Type = c("Conifer", "Deciduous"))
#' leadingByPixelGroup(x, sppEquiv = sppEquiv, sppEquivCol = "LandR")
leadingByPixelGroup <- function(
  x,
  vegLeadingProportion = NULL,
  mixedType = 2,
  sppEquiv = NULL,
  sppEquivCol,
  pixelGroupColName = "pixelGroup",
  leadingBasedOn = if ("B" %in% colnames(x)) "B" else "cover"
) {
  stopifnot(mixedType %in% 0:2, length(mixedType) == 1)
  vegLeadingProportion <- .leadingProp(vegLeadingProportion, mixedType)
  if (mixedType == 2 && !sppEquivCol %in% colnames(sppEquiv)) {
    stop(sppEquivCol, " is not in sppEquiv. Please pass an existing sppEquivCol")
  }

  if (.canSumLeadingByIntegerKey(x, mixedType, sppEquiv, sppEquivCol, pixelGroupColName, leadingBasedOn)) {
    .leadingByIntegerKey(x, vegLeadingProportion, sppEquiv, sppEquivCol, pixelGroupColName)
  } else {
    .leadingByPixelGroupSlow(
      x, vegLeadingProportion, mixedType, sppEquiv, sppEquivCol, pixelGroupColName, leadingBasedOn
    )
  }
}

## The fast algorithm builds one integer key from (pixelGroup, speciesCode), so it needs a
## positive whole-number pixelGroup, a factor speciesCode without NA, no NA in B (the slow
## algorithm's total for a one-cohort group is then NA, not 0), and one Type per species.
.canSumLeadingByIntegerKey <- function(x, mixedType, sppEquiv, sppEquivCol, pixelGroupColName,
                                       leadingBasedOn) {
  if (mixedType != 2 || !identical(leadingBasedOn, "B")) return(FALSE)
  pg <- x[[pixelGroupColName]]
  sc <- x[["speciesCode"]]
  if (!is.factor(sc) || !is.numeric(pg) || anyNA(sc) || anyNA(pg) || anyNA(x[["B"]])) return(FALSE)
  if (length(pg) && (min(pg) < 1 || max(pg) * nlevels(sc) >= .Machine$integer.max)) return(FALSE)
  if (is.double(pg) && any(pg != trunc(pg))) return(FALSE)
  sppEq <- unique(data.table(speciesCode = sppEquiv[[sppEquivCol]], Type = sppEquiv[["Type"]]))
  !anyDuplicated(sppEq$speciesCode)
}

## The algorithm the vegTypeGenerator and vegTypeMapGenerator.data.table used to carry one copy each of.
.leadingByPixelGroupSlow <- function(x, vegLeadingProportion, mixedType, sppEquiv, sppEquivCol,
                                     pixelGroupColName, leadingBasedOn) {
  pgdAndSc <- c(pixelGroupColName, "speciesCode")
  pgdAndScAndLeading <- c(pgdAndSc, leadingBasedOn)
  totalOfLeadingBasedOn <- paste0("total", leadingBasedOn)
  speciesOfLeadingBasedOn <- paste0("speciesGroup", leadingBasedOn)

  ## 1. Find length of each pixelGroup -- don't include pixelGroups in "by" that have only 1 cohort: N = 1
  pixelGroupData1 <- x[, list(N = .N), by = pixelGroupColName]

  ## Calculate speciesProportion from cover or B
  pixelGroupData1 <- x[, ..pgdAndScAndLeading][pixelGroupData1, on = pixelGroupColName]
  set(pixelGroupData1, NULL, totalOfLeadingBasedOn, pixelGroupData1[[leadingBasedOn]])

  if (identical(leadingBasedOn, "cover")) {
    pixelGroupData1[
      N != 1,
      (totalOfLeadingBasedOn) := sum(cover, na.rm = TRUE),
      by = pixelGroupColName
    ]
    pixelGroupData1 <- pixelGroupData1[,
      list(sum(cover, na.rm = TRUE), totalcover[1]),
      by = pgdAndSc
    ]
  } else {
    pixelGroupData1[N != 1, (totalOfLeadingBasedOn) := sum(B, na.rm = TRUE), by = pixelGroupColName]
    pixelGroupData1 <- pixelGroupData1[, list(sum(B, na.rm = TRUE), totalB[1]), by = pgdAndSc]
  }
  setnames(
    pixelGroupData1,
    old = c("V1", "V2"),
    new = c(speciesOfLeadingBasedOn, totalOfLeadingBasedOn)
  )

  set(
    pixelGroupData1,
    NULL,
    "speciesProportion",
    pixelGroupData1[[speciesOfLeadingBasedOn]] / pixelGroupData1[[totalOfLeadingBasedOn]]
  )

  setorderv(pixelGroupData1, pixelGroupColName)

  pixelGroupData <- pixelGroupData1
  rm(pixelGroupData1)

  ## Determine "mixed" -----------------------------------
  if (mixedType == 0) {
    ## 1. sort on pixelGroup and speciesProportion, reverse so 1st row of each pixelGroup is the largest
    ## 2. Keep only first row in each pixelGroup
    pixelGroupData3 <- pixelGroupData[, list(
      speciesCode,
      get(pixelGroupColName),
      speciesProportion
    )]
    setnames(pixelGroupData3, "V2", pixelGroupColName)
    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName]
    setnames(pixelGroupData3, "speciesCode", "leading")
  } else if (mixedType == 1) {
    ## create "mixed" class (2019-05: faster than previous below)
    ## 1. anything with >= vegLeadingProportion is "pure"
    ## 2. sort on pixelGroup and speciesProportion, reverse so that 1st row of each pixelGroup is the largest
    ## 3. Keep only first row in each pixelGroup
    ## 4. change column names and convert pure to mixed ==> mixed <- !pure
    pixelGroupData3 <- pixelGroupData[, list(
      pure = speciesProportion >= vegLeadingProportion,
      speciesCode,
      get(pixelGroupColName),
      speciesProportion
    )]
    setnames(pixelGroupData3, "V3", pixelGroupColName)
    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName]
    pixelGroupData3[pure == FALSE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    pixelGroupData3[, pure := !pure]
    setnames(pixelGroupData3, "pure", "mixed")
  } else {
    sppEq <- data.table(sppEquiv[[sppEquivCol]], sppEquiv[["Type"]])

    names(sppEq) <- c("speciesCode", "Type")
    setkey(pixelGroupData, speciesCode)

    # don't need all columns now
    colsToDelete <- c("rasterToMatch", leadingBasedOn, totalOfLeadingBasedOn)
    colsToDelete <- colsToDelete[colsToDelete %in% colnames(pixelGroupData)]
    set(pixelGroupData, NULL, colsToDelete, NULL)
    pixelGroupData3 <- merge(pixelGroupData, sppEq[!duplicated(sppEq)], all.x = TRUE)
    setkeyv(pixelGroupData, pgdAndSc)

    setkeyv(pixelGroupData3, pgdAndSc)
    pixelGroupData3[,
      mixed := .isMixedwood(speciesProportion, Type, vegLeadingProportion),
      by = pixelGroupColName
    ]

    setorderv(pixelGroupData3, cols = c(pixelGroupColName, "speciesProportion"), order = -1L)
    set(pixelGroupData3, NULL, "speciesProportion", NULL)
    set(pixelGroupData3, NULL, "Type", NULL)
    pixelGroupData3 <- pixelGroupData3[, .SD[1], by = pixelGroupColName] ## sp. w/ highest prop. per pixelGroup
    pixelGroupData3[mixed == TRUE, speciesCode := "Mixed"]
    setnames(pixelGroupData3, "speciesCode", "leading")
    set(pixelGroupData3, NULL, "leading", factor(pixelGroupData3[["leading"]]))
  }

  pixelGroupData3
}

## mixedType = 2 on B. Same result as the slow algorithm, without its 10M-row join and copy, its
## grouped passes over every cohort row, or one R call of .isMixedwood() per pixelGroup.
.leadingByIntegerKey <- function(x, vegLeadingProportion, sppEquiv, sppEquivCol, pixelGroupColName) {
  pg <- x[[pixelGroupColName]]
  sc <- x[["speciesCode"]]
  lv <- levels(sc)
  nSp <- length(lv)

  ## species-within-pixelGroup sums: one integer key, one grouped pass
  d <- data.table(k = (pg - 1L) * nSp + as.integer(sc), B = x[["B"]])[,
    list(speciesGroupB = sum(B, na.rm = TRUE)),
    by = k
  ]
  set(d, NULL, "pixelGroup", (d$k - 1L) %/% nSp + 1L)
  ## character, as the merge() on sppEquiv's character speciesCode leaves it in the slow
  ## algorithm: the (pixelGroup, speciesCode) sort and the levels of `leading` are then in
  ## character order, with "Mixed" first
  set(d, NULL, "speciesCode", lv[(d$k - 1L) %% nSp + 1L])
  tot <- d[, list(totalB = sum(speciesGroupB, na.rm = TRUE)), by = pixelGroup]
  set(d, NULL, "totalB", tot$totalB[match(d$pixelGroup, tot$pixelGroup)])
  set(d, NULL, "speciesProportion", d$speciesGroupB / d$totalB)
  set(d, NULL, c("k", "totalB"), NULL)
  setnames(d, "pixelGroup", pixelGroupColName)
  setorderv(d, c(pixelGroupColName, "speciesCode"))

  sppEq <- unique(data.table(speciesCode = sppEquiv[[sppEquivCol]], Type = sppEquiv[["Type"]]))
  set(d, NULL, "Type", sppEq$Type[match(d$speciesCode, sppEq$speciesCode)])

  ## broadleaf proportion per pixelGroup. A single deciduous species: its proportion (the grouped
  ## sum of one value is exact). Two or more: base sum(), as .isMixedwood() does, whose
  ## accumulator is wider than the grouped sum's; the last bit can matter at the threshold.
  dec <- d[Type %in% "Deciduous",
    list(broadleaf = sum(speciesProportion, na.rm = TRUE), nDec = .N),
    by = pixelGroupColName
  ]
  multi <- dec$nDec >= 2L
  if (any(multi)) {
    pgMulti <- dec[[pixelGroupColName]][multi]
    rows <- d[Type %in% "Deciduous" & d[[pixelGroupColName]] %in% pgMulti]
    byGroup <- split(
      rows$speciesProportion,
      factor(match(rows[[pixelGroupColName]], pgMulti), levels = seq_along(pgMulti))
    )
    dec$broadleaf[multi] <- vapply(byGroup, sum, numeric(1), na.rm = TRUE, USE.NAMES = FALSE)
  }
  bl <- dec$broadleaf[match(d[[pixelGroupColName]], dec[[pixelGroupColName]])]
  bl[is.na(bl)] <- 0
  set(d, NULL, "mixed", bl < vegLeadingProportion & bl > 1 - vegLeadingProportion)
  set(d, NULL, "Type", NULL)

  ## leading = highest proportion, ties to the first speciesCode in character order
  setorderv(d, c(pixelGroupColName, "speciesProportion", "speciesCode"), order = c(-1L, -1L, 1L))
  keep <- c(pixelGroupColName, "speciesCode", "speciesGroupB", "mixed")
  pixelGroupData3 <- d[!duplicated(d[[pixelGroupColName]]), keep, with = FALSE]
  pixelGroupData3[mixed == TRUE, speciesCode := "Mixed"]
  setnames(pixelGroupData3, "speciesCode", "leading")
  set(pixelGroupData3, NULL, "leading", factor(pixelGroupData3[["leading"]]))
  pixelGroupData3
}
