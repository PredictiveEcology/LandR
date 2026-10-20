## ---------------------------------------------------------------------------
## Messaging for users without SCANFI access
##
## SCANFI is distributed through a Google Drive folder shared with collaborators
## rather than published openly. A user who has not been granted access gets only
## reproducible's generic failure -- "Could not access the Google Drive resource
## ... Client error: (404) Not Found" -- which reads like a dead link and invites
## a hunt for a broken link. It is a permissions problem, so say so, and name the
## ways out. SCANFI v2 is also on a public mirror; a failure there is a network problem.
## With no internet at all (an HPC compute node) Drive and the mirror fail the same
## way, and reproducible wraps the Drive one in the same "Could not access" text, so
## that case is told apart first: a local copy fixes it, access does not.
## ---------------------------------------------------------------------------

## Official SCANFI distribution point; where to download the files yourself.
.scanfiAccessURL <- "https://opendata.nfis.org/"

## Does this error look like Google Drive refusing us a file? Matches both the
## wrapper reproducible raises and the underlying googledrive status text, since
## which of them surfaces depends on where the request failed. A bare status code
## is enough here because this is only ever applied to a known Drive fetch.
.isGoogleDriveAccessError <- function(e) {
  msg <- paste(conditionMessage(e), collapse = "\n")

  grepl("google drive", msg, ignore.case = TRUE) ||
    grepl("(^|[^0-9])(401|403|404)([^0-9]|$)", msg) ||
    grepl("not found|permission|access denied|forbidden|unauthorized", msg, ignore.case = TRUE)
}

.scanfiAccessMessage <- function(
  e,
  what,
  dataYear = NULL,
  dataVersion = NULL,
  urlArg = NULL,
  alternatives = NULL,
  url = NULL,
  offline = FALSE
) {
  qualifier <- paste(c(dataVersion, dataYear), collapse = " ")
  original <- paste(conditionMessage(e), collapse = "\n")

  ## the error says which server refused: Drive, whatever remap was meant to apply, or the mirror
  fromDrive <- grepl("drive.google.com|google drive", original, ignore.case = TRUE)
  fromMirror <- !fromDrive && .remapsToMirror(url, dataVersion)
  ownCopy <- if (!is.null(urlArg)) paste0("pass your own copy with `", urlArg, "`")

  mirrorWay <- if (identical(dataVersion, "V2")) {
    if (!is.null(getOption("reproducible.urlRemap"))) {
      paste0(
        "check that your reproducible.urlRemap covers this file; LandR's list of mirrored ",
        "SCANFI files is system.file(\"extdata\", \"arbutus_manifest_SCANFI_v2.csv\", package = \"LandR\")"
      )
    } else if (!isTRUE(getOption("LandR.scanfiMirror", TRUE))) {
      "set options(LandR.scanfiMirror = TRUE) to use the public SCANFI mirror"
    }
  }
  ways <- if (offline) {
    c(
      paste0("run once with internet: the file is then kept in destinationPath (and in ",
             "getOption(\"reproducible.destinationPathShared\") if set) and read from there"),
      ownCopy
    )
  } else if (!fromMirror) {
    c(
      mirrorWay,
      paste0("ask for access to the SCANFI Google Drive folder, then run ",
             "googledrive::drive_auth() with that account"),
      paste0("download it from ", .scanfiAccessURL,
             if (!is.null(urlArg)) paste0(" and pass it with `", urlArg, "`")),
      if (length(alternatives)) paste0("use ", paste0('dataSource = "', alternatives, '"', collapse = " or "))
    )
  }

  paste0(
    "Could not download ", what,
    if (nzchar(qualifier)) paste0(" (", qualifier, ")"),
    if (offline) {
      ": no internet connection, so the server could not be reached."
    } else if (fromMirror) {
      " from the SCANFI mirror. This is usually a network problem; retry later."
    } else {
      paste0(
        " from Google Drive. SCANFI's Drive folder is not public, so this usually means ",
        "your Google account has no access, not a broken link."
      )
    },
    if (length(ways)) paste0("\nTo get it:\n", paste0("  * ", ways, collapse = "\n")),
    "\nOriginal error:\n",
    paste0("  ", strsplit(original, "\n")[[1]], collapse = "\n")
  )
}

## Wrap a SCANFI download so an access failure explains itself. Anything that is
## neither a connection failure nor a Drive access problem is re-thrown untouched.
## `enabled` lets functions that serve several data sources wrap the call in place,
## without duplicating it for the SCANFI and non-SCANFI branches.
.withSCANFIAccess <- function(
  expr,
  what,
  dataYear = NULL,
  dataVersion = NULL,
  urlArg = NULL,
  alternatives = NULL,
  url = NULL,
  enabled = TRUE
) {
  if (!isTRUE(enabled)) {
    return(expr)
  }
  tryCatch(.withScanfiRemap(expr), error = function(e) {
    offline <- reproducible::isConnectionError(e)
    if (!offline && !.isGoogleDriveAccessError(e)) {
      stop(e)
    }
    stop(.scanfiAccessMessage(e, what, dataYear, dataVersion, urlArg, alternatives, url, offline),
         call. = FALSE)
  })
}

## ---------------------------------------------------------------------------
## Messaging for SCANFI v3 access failures
##
## V3 is served directly by NRCan (no Google Drive layer), so a failure here is
## usually a network hiccup or an outage, not a permissions problem, and the
## guidance differs from .withSCANFIAccess() above.
## ---------------------------------------------------------------------------

## Where else to get v3 if the https server is unreachable.
.scanfiV3FTPURL <- "ftp://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/SCANFI/v3/"

## Does this error look like the v3 https server refusing or dropping us? A bare status
## code is enough here because this is only ever applied to a known SCANFI v3 fetch.
.isSCANFIv3AccessError <- function(e) {
  msg <- paste(conditionMessage(e), collapse = "\n")

  grepl("(^|[^0-9])(403|404)([^0-9]|$)", msg) ||
    grepl(
      "forbidden|could not open|cannot open|curl error|timed? ?out|connection",
      msg,
      ignore.case = TRUE
    )
}

.scanfiV3AccessMessage <- function(e, what, dataYear = NULL) {
  qualifier <- paste(c("V3", dataYear), collapse = " ")

  paste0(
    "Could not download ",
    what,
    if (nzchar(qualifier)) paste0(" (", qualifier, ")") else "",
    ".\n\n",
    "SCANFI v3 is read from NRCan's server (ftp.maps.canada.ca) as one large\n",
    "Cloud-Optimized GeoTIFF per year; the server can be slow or unreachable under load.\n\n",
    "To resolve it:\n",
    "  * retry later\n",
    "  * or fetch the year's file from the FTP mirror instead: ", .scanfiV3FTPURL, "\n\n",
    "Original error:\n",
    paste0("  ", strsplit(paste(conditionMessage(e), collapse = "\n"), "\n")[[1]], collapse = "\n")
  )
}

## Wrap a SCANFI v3 read so an access failure explains itself. Anything else is
## re-thrown untouched.
.withSCANFIv3Access <- function(expr, what, dataYear = NULL, enabled = TRUE) {
  if (!isTRUE(enabled)) {
    return(expr)
  }

  tryCatch(expr, error = function(e) {
    if (!.isSCANFIv3AccessError(e)) {
      stop(e)
    }
    stop(.scanfiV3AccessMessage(e, what, dataYear), call. = FALSE)
  })
}
