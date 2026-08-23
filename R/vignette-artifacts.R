## ---------------------------------------------------------------------------
## Vignette artifact manifest
##
## The SCANFI summary vignette (`vignettes/scanfi-summary.Rmd`) reports exact
## national statistics derived from 30 m rasters that cannot be processed on a
## CI runner: a single pass reads tens of GB, the runners have ~14 GB of disk,
## and the polygon sources throttle CI IPs.
##
## Rather than freezing the rendered output, the expensive step is cut at its
## natural seam. `calc_raster_counts()` is the only thing that touches the full
## raster; its output -- a value-frequency table of a few thousand rows -- is
## precomputed on a workstation and committed. Everything downstream (summary
## statistics, histograms, maps drawn from a coarsened display raster) is
## recomputed when the vignette is built, in seconds, offline. Those parts
## therefore cannot go stale.
##
## What *can* go stale is the committed counts table, if the code that produced
## it changes. So it is pinned to a hash of the normalized source of the
## functions it depends on. `removeSource()` + `deparse()` strips comments and
## normalises whitespace, so roxygen edits and reformatting do not trigger false
## alarms; only a real change to the code does.
##
## Note what is deliberately absent below: `stats_from_counts()` and
## `weighted_quantile()` are an exact function of the counts table, so the
## vignette recomputes them at build time. A change to the quantile definition
## must NOT force a multi-hour regeneration of a counts table that is still
## perfectly valid.
## ---------------------------------------------------------------------------

## Function sets whose source determines each committed artifact.
.artifact_fns <- list(
  counts = c(
    "prep_polygons",
    "calc_raster_counts"
  )
)

## Normalized source of a single function: `removeSource()` drops the srcref
## that carries comments and original formatting, so `deparse()` returns
## canonical code. Two functions that differ only in comments or whitespace
## normalise identically.
.normalize_src <- function(f) {
  paste(deparse(removeSource(f)), collapse = "\n")
}

## Hash the normalized source of each function set. Function *names* are folded
## into the hashed payload, so adding or removing a helper also invalidates the
## artifact even if the remaining sources are untouched.
.artifact_hashes <- function() {
  stopifnot(requireNamespace("digest", quietly = TRUE))

  ns <- asNamespace("LandR")

  vapply(.artifact_fns, function(fns) {
    src <- vapply(sort(fns), function(nm) .normalize_src(get(nm, envir = ns)), character(1))
    digest::digest(paste(names(src), src, sep = "\n", collapse = "\n\n"))
  }, character(1))
}

## Path to the committed manifest; "" when the artifacts have not been generated.
.artifact_manifest_file <- function() {
  system.file("extdata", "vignette-artifacts.dcf", package = "LandR")
}

.read_artifact_manifest <- function(path = .artifact_manifest_file()) {
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  as.data.frame(read.dcf(path), stringsAsFactors = FALSE)
}

## Record the hashes alongside whatever provenance the caller wants to note
## (data source, year, polygon vintage, ...). Called by
## `data-raw/precompute-vignettes.R`; not used at package run time.
.write_artifact_manifest <- function(path, provenance = list()) {
  hashes <- as.list(.artifact_hashes())
  names(hashes) <- paste0(names(hashes), "_hash")

  fields <- c(
    list(
      Generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      LandRVersion = as.character(utils::packageVersion("LandR")),
      RVersion = R.version.string
    ),
    lapply(provenance, as.character),
    hashes
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ## `width` keeps long values (URLs) on one line: read.dcf folds continuation
  ## lines back in with their leading whitespace, which renders badly.
  write.dcf(as.data.frame(fields, stringsAsFactors = FALSE), file = path, width = 10000L)

  invisible(path)
}
