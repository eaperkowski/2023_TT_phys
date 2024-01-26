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
  dplyr::select(-subplot.y, -X, 
                plot, subplot = subplot.x, composite:days_deployed,
                phosphate_ug = phosphate_ppm, nitrate_ug = nitrate_ppm,
                ammonium_ug = ammonium_ppm) %>%
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
  group_by(canopy, plot, composite, days_deployed) %>%
  summarize(phosphate_ug = mean(phosphate_ug, na.rm = TRUE),
            nitrate_ug = mean(nitrate_ug, na.rm = TRUE),
            ammonium_ug = mean(ammonium_ug, na.rm = TRUE)) %>%
  mutate(phosphate_ug = na_if(phosphate_ug, NaN),
         nitrate_ug = na_if(nitrate_ug, NaN),
         ammonium_ug = na_if(ammonium_ug, NaN)) %>%
  mutate(inorg_n_ug = ifelse(!is.na(nitrate_ug) & !is.na(ammonium_ug),
                             nitrate_ug + ammonium_ug, NA),
         phosphate_ug_day = phosphate_ug / days_deployed,
         nitrate_ug_day = nitrate_ug / days_deployed,
         ammonium_ug_day = ammonium_ug / days_deployed,
         inorg_n_ug_day = inorg_n_ug / days_deployed)

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



