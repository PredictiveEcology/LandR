## ---------------------------------------------------------------------------
## Messaging for users without SCANFI access
##
## SCANFI is distributed through a Google Drive folder shared with collaborators
## rather than published openly. A user who has not been granted access gets only
## reproducible's generic failure -- "Could not access the Google Drive resource
## ... Client error: (404) Not Found" -- which reads like a dead link and invites
## a hunt for a mirror that does not exist. It is a permissions problem, so say
## so, and name the ways out.
## ---------------------------------------------------------------------------

## Official SCANFI distribution point; where to ask for access.
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

.scanfiAccessMessage <- function(e, what, dataYear = NULL, dataVersion = NULL,
                                 urlArg = NULL, alternatives = NULL) {
  qualifier <- paste(c(dataVersion, dataYear), collapse = " ")

  ways <- c(
    paste0(
      "request access to the SCANFI data from the SCANFI team (see ",
      .scanfiAccessURL, "), then re-authenticate with googledrive::drive_auth()"
    ),
    if (!is.null(urlArg)) {
      paste0("or, if you already have the data, pass it via `", urlArg, "`")
    },
    if (length(alternatives)) {
      paste0(
        "or use a different data source: ",
        paste0('dataSource = "', alternatives, '"', collapse = " or ")
      )
    }
  )

  paste0(
    "Could not download ", what,
    if (nzchar(qualifier)) paste0(" (", qualifier, ")") else "", ".\n\n",
    "SCANFI is distributed through a restricted Google Drive folder, so this is\n",
    "usually a permissions problem rather than a broken link: the Google account\n",
    "you are authenticated as may not have been granted access.\n\n",
    "To resolve it:\n",
    paste(
      vapply(
        ways,
        function(w) paste(strwrap(w, width = 76, initial = "  * ", prefix = "    "),
                          collapse = "\n"),
        character(1)
      ),
      collapse = "\n"
    ), "\n\n",
    "Original error:\n",
    paste0("  ", strsplit(paste(conditionMessage(e), collapse = "\n"), "\n")[[1]],
           collapse = "\n")
  )
}

## Wrap a SCANFI download so an access failure explains itself. Anything that is
## not a Drive access problem is re-thrown untouched. `enabled` lets functions
## that serve several data sources wrap the call in place, without duplicating it
## for the SCANFI and non-SCANFI branches.
.withSCANFIAccess <- function(expr, what, dataYear = NULL, dataVersion = NULL,
                              urlArg = NULL, alternatives = NULL, enabled = TRUE) {
  if (!isTRUE(enabled)) {
    return(expr)
  }

  tryCatch(
    expr,
    error = function(e) {
      if (!.isGoogleDriveAccessError(e)) {
        stop(e)
      }
      stop(
        .scanfiAccessMessage(e, what, dataYear, dataVersion, urlArg, alternatives),
        call. = FALSE
      )
    }
  )
}
