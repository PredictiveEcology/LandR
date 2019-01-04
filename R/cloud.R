if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("checksumsFilename", "checksumsID", "id"))
}

#' Basic tool for using cloud-based caching
#'
#' Very experimental
#'
#' @param toDigest The R object to consider, e.g., all the arguments to a function
#'
#' @param checksumsFileID A google file ID where the checksums data.table is located,
#'   provided as a character string
#'
#' @param writeChecksumsFolderID The google folder ID where a new checksums file should
#'    be written. This will only be used if \code{checksumsFileID} is not provided
#'   provided as a character string
#'
#' @export
#' @importFrom fastdigest fastdigest
#' @importFrom googledrive as_id drive_download drive_upload
checkCloud <- function(toDigest, checksumsFileID, writeChecksumsFolderID = NULL) {
  dig <- fastdigest::fastdigest(toDigest)
  if (!is.null(writeChecksumsFolderID)) {
    checksums <- data.table(hash = character(), id = character(), time = character())
    checksumsFilename <- file.path(tempdir(), "biomassModel.rds")
    saveRDS(checksums, file = checksumsFilename)
    drive_upload(checksums, path = as_id(writeChecksumsFolderID), name = "checksums")
  }
  checksums <- downloadChecksums(checksumsFileID)
  hashExists <- checksums$hash == dig
  outBiomass <- if (isTRUE(any(hashExists))) {
    biomassModelFilename2 <- tempfile(fileext = ".rds");
    drive_download(as_id(checksums[hashExists, id]), path = biomassModelFilename2)
    readRDS(biomassModelFilename2)
  } else {
    NULL
  }
  return(list(object = outBiomass, digest = dig, checksums = checksums))
}

#' Basic tool for using cloud-based caching
#'
#' Very experimental
#'
#' @param object The R object to write to cloud
#'
#' @param digest The hash of the input arguments, outputted from \code{checkCloud}
#'
#' @param checksums A \code{data.table} that is outputted from \code{checkCloud} that
#'   is the the checksums file
#'
#' @param cloudFolderID The google folder ID where a new object should
#'    be written
#'
#' @export
#' @importFrom googledrive as_id drive_update drive_upload
#' @seealso \code{\link{checkCloud}}
writeCloud <- function(object, digest, cloudFolderID = NULL, checksums) {
  if (!is.null(cloudFolderID)) {
    biomassModelFile <- tempfile(fileext = ".rds")
    saveRDS(outBiomass, file = biomassModelFile)
    # checksums <- downloadChecksums(checksumsFileID)
    uploadRes <- drive_upload(biomassModelFile, path = as_id(cloudFolderID),
                              name = paste0("biomassModel", NROW(checksums) + 1))
    checksums <- rbindlist(
      list(checksums,
           data.table(hash = digest, id = uploadRes$id, time = as.character(Sys.time()))))
    saveRDS(checksums, file = checksumsFilename)
    drive_update(as_id(checksumsID), media = checksumsFilename)
  } else {
    stop("cloudFolder must be provided as a google id string to a folder (not a file)")
  }
}

#' Download checksums file from Google Drive
#'
#' Very experimental
#'
#' @param checksumsFileID the Google Drive file id for the remote checksums file
#'
#' @importFrom googledrive as_id drive_download
#' @keywords internal
downloadChecksums <- function(checksumsFileID) {
  checksumsFilename <- tempfile(fileext = ".rds");
  drive_download(as_id(checksumsFileID), path = checksumsFilename)
  checksums <- readRDS(checksumsFilename)
}
