## Libraries
library(tidyverse)
library(lubridate)

## Read in files
files <- list.files("../data/tomst_probe_raw/", 
                    recursive = T,
                    pattern = "\\.csv$", 
                    full.names = T)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))
tomst <- lapply(files, read.csv) %>%
  reshape::merge_all()

## Parameters for quadratic curve (calibration set
## for silt loam from TOMST website)
quad_params <- data.frame(a = 0.000000017,
                          b = 0.000118119,
                          c = -0.101168511)

## Calculate volumetric soil moisture
tomst %>%
  filter(soil_moisture_raw > 1000) %>%
  mutate(soil_moisture_volum = quad_params$a * (soil_moisture_raw)^2 +
           quad_params$b * soil_moisture_raw + quad_params$c,
         date = ymd_hm(date))



