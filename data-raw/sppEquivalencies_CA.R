library(data.table)
sppEquivalencies_CA <- as.data.table(read.csv("data-raw/sppEquivalencies_CA.csv",
                                              stringsAsFactors = FALSE, header = TRUE))
str(sppEquivalencies_CA)

usethis::use_data(sppEquivalencies_CA, overwrite = TRUE)
