## Helpers for the Google Drive fallback in prepSpeciesLayers_*().

#' Is this a Google Drive URL?
#'
#' A Drive source needs no reachability ping before use: the fallback that the
#' ping would select is also Drive, and under packet loss the ping fails while
#' the listing itself would have succeeded. On 2026-09-07 that ping sent jobs
#' into the fallback repeatedly, two hours into each run.
#'
#' @param url character
#' @return logical(1)
#' @keywords internal
.isGoogleDriveUrl <- function(url) {
  is.character(url) && length(url) == 1L && !is.na(url) &&
    grepl("^https?://(drive|docs)\\.google\\.com/", url)
}

#' Link to a folder found by name in a Drive listing
#'
#' The fallback taken when a primary source cannot be reached. It fails loudly
#' when the folder is not in the listing: `googledrive::drive_link(character(0))`
#' is `character(0)`, and `googledrive::drive_ls(character(0))` reports only
#' "Parent specified via `path` is invalid: Does not exist", which says nothing
#' about which folder was wanted or where.
#'
#' @param driveDT data.table with columns `name` and `id`, from
#'   `googledrive::drive_ls()` on the shared folder.
#' @param driveFolder character; the folder name wanted.
#' @param shared_drive_url character; where it was looked for (for the message).
#' @return character(1) folder URL.
#' @keywords internal
.driveFolderLink <- function(driveDT, driveFolder, shared_drive_url) {
  ids <- as.character(driveDT[["id"]][driveDT[["name"]] == driveFolder])
  if (!length(ids)) {
    stop("Google Drive fallback: no folder named '", driveFolder, "' in ", shared_drive_url,
         ". Folders there: ", paste(driveDT[["name"]], collapse = ", "), call. = FALSE)
  }
  paste0("https://drive.google.com/drive/folders/", ids[1])
}
