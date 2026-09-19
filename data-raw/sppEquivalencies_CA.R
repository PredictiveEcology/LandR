library(data.table)

sppEquivalencies_CA <- fread("data-raw/sppEquivalencies_CA.csv",
  stringsAsFactors = FALSE, header = TRUE
)
str(sppEquivalencies_CA)

## subset kNN species names to existing layers
if (FALSE) { ## run this manually
  library(RCurl)
  library(XML)

  fileURLs2001 <- getURL(
    paste0(
      "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
      "canada-forests-attributes_attributs-forests-canada/",
      "2001-attributes_attributs-2001/"
    ),
    dirlistonly = TRUE,
    .opts = list(followlocation = TRUE)
  )
  fileNames2001 <- getHTMLLinks(fileURLs2001)

  fileURLs2011 <- getURL(
    paste0(
      "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
      "canada-forests-attributes_attributs-forests-canada/",
      "2011-attributes_attributs-2011/"
    ),
    dirlistonly = TRUE,
    .opts = list(followlocation = TRUE)
  )
  fileNames2011 <- getHTMLLinks(fileURLs2011)

  fileNames <- union(fileNames2001, fileNames2011)
  fileNames <- grep("(Species|SpeciesGroups)_.*\\.tif$", fileNames, value = TRUE)

  allSpp <- fileNames |>
    sub("_v1\\.tif", "", x = _) |>
    sub(".*(Species|SpeciesGroups)_", "", x = _)
  allSpp <- unique(allSpp)

  ## remove species that have no layers
  sppEquivalencies_CA[!KNN %in% allSpp, KNN := ""]
  ## checks
  setdiff(sppEquivalencies_CA$KNN, allSpp)
  setdiff(allSpp, sppEquivalencies_CA$KNN)
  ## Ceres: not sure what these are: "Genc_Spp", "Genh_Spp"
}

usethis::use_data(sppEquivalencies_CA, overwrite = TRUE)
