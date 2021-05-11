## code to prepare `DATASET` dataset goes here
library(lipdR)
library(tidyverse)

#Palmyra coral record

palmyra <- lipdR::readLipd("https://lipdverse.org/Pages2kTemperature/2_1_0/Ocn-Palmyra.Cobb.2013.lpd") %>%
  extractTs() %>%
  ts2tibble()

usethis::use_data(palmyra, overwrite = TRUE)


#Great basin tree ring record

great_basin <- lipdR::readLipd("https://lipdverse.org/Pages2kTemperature/2_1_0/NAm-GreatBasin.Salzer.2013.lpd") %>%
  extractTs() %>%
  ts2tibble()

usethis::use_data(great_basin, overwrite = TRUE)


#Marine record
odp658c <- lipdR::readLipd("https://lipdverse.org/HoloceneAbruptChange/0_10_0/ODP658C.deMenocal.1997.lpd") %>%
  extractTs() %>%
  ts2tibble()

usethis::use_data(odp658c, overwrite = TRUE)


#Ice core record

ngrip <- lipdR::readLipd("https://lipdverse.org/HoloceneAbruptChange/0_10_0/NGRIP.NGRIP.2004.lpd") %>%
  extractTs() %>%
  ts2tibble()

usethis::use_data(ngrip, overwrite = TRUE)


#full ensemble data

geob7926_2 <- lipdR::readLipd("https://github.com/nickmckay/Temperature12k/raw/master/ScientificDataAnalysis/lipdFilesWithEnsembles/GeoB7926_2.Kim.2012.lpd") %>%
  extractTs() %>%
  ts2tibble()


usethis::use_data(geob7926_2, overwrite = TRUE)


