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
#' The fallback taken when a primary source cannot be reached. It exists only to
#' fail loudly when the folder is not in the listing: `drive_link(character(0))`
#' is `character(0)`, and the `drive_ls()` that follows reports only "Parent
#' specified via `path` is invalid: Does not exist", which says nothing about
#' which folder was wanted or where. The link itself comes from
#' [googledrive::drive_link()] on the matching row, which reads the `webViewLink`
#' that `drive_ls()` already fetched -- no further Drive call.
#'
#' @param driveLs a `dribble`, from [googledrive::drive_ls()] on the shared folder.
#' @param driveFolder character; the folder name wanted.
#' @param shared_drive_url character; where it was looked for (for the message).
#' @return character(1) folder URL.
#' @keywords internal
.driveFolderLink <- function(driveLs, driveFolder, shared_drive_url) {
  i <- which(driveLs[["name"]] == driveFolder)
  if (!length(i)) {
    stop("Google Drive fallback: no folder named '", driveFolder, "' in ", shared_drive_url,
         ". Folders there: ", paste(driveLs[["name"]], collapse = ", "), call. = FALSE)
  }
  googledrive::with_drive_quiet(googledrive::drive_link(driveLs[i[1], ]))
}
