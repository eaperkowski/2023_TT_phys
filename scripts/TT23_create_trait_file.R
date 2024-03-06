#####################################################################
# Load libraries and data files
#####################################################################
# Libraries
library(dplyr)
library(tidyverse)

# Read data files for gas exchange
phys <- read.csv("../data/raw_data/TT23_gasExchange.csv")
multispeq <- read.csv("../data/raw_data/TT23_multispeq_data.csv")

#####################################################################
# Compile resin strip data into single object, then standardize
# for duration in the field
#####################################################################
phos_nit_corrected <- read.csv("../data/raw_data/TT23_resin_strip_corrected_po4_no3.csv")
ammonium_corrected <- read.csv("../data/raw_data/TT23_resin_strip_corrected_ammonium.csv")

nutrients <- phos_nit_corrected %>%
  full_join(ammonium_corrected, by = c("plot", "composite", "canopy")) %>%
  dplyr::select(plot, composite:days_deployed,
                phosphate_ppm, nitrate_ppm, ammonium_ppm) %>%
  mutate(days_deployed = ifelse(canopy == "post_closure" & plot == 3,
                                29, 
                                ifelse(canopy == "post_closure" & plot == 5,
                                       28, 
                                       ifelse(canopy == "post_closure" & plot == 7,
                                              29,
                                              ifelse(canopy == "pre_closure" & plot == 3,
                                                     41,
                                                     ifelse(canopy == "pre_closure" & plot == 5,
                                                            42,
                                                            ifelse(canopy == "pre_closure" & plot == 7,
                                                                   35, NA))))))) %>%
  group_by(canopy, plot, composite) %>%
  summarize(phosphate_ppm = mean(phosphate_ppm, na.rm = TRUE),
            nitrate_ppm = mean(nitrate_ppm, na.rm = TRUE),
            ammonium_ppm = mean(ammonium_ppm, na.rm = TRUE)) %>%
  mutate(phosphate_ppm = ifelse(phosphate_ppm == "NaN", NA, phosphate_ppm),
         nitrate_ppm = ifelse(nitrate_ppm == "NaN", NA, nitrate_ppm),
         ammonium_ppm = ifelse(ammonium_ppm == "NaN", NA, ammonium_ppm)) %>%
  mutate(inorg_n_ppm = ifelse(!is.na(nitrate_ppm) & !is.na(ammonium_ppm),
                             nitrate_ppm + ammonium_ppm, NA))

#####################################################################
# Join multispeq and physiology data set
#####################################################################
phys.total <- phys %>%
  full_join(multispeq) %>%
  full_join(nutrients, by = c("plot", "composite", "canopy")) %>%
  filter(!is.na(id) & !is.na(anet)) %>%
  mutate(stom.lim = ifelse(stom.lim < 0, NA, stom.lim))

#####################################################################
# Write compiled data file
#####################################################################
write.csv(phys.total, "../data/TT23_compiled_datasheet.csv", row.names = FALSE)



