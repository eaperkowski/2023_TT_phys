#####################################################################
# Load libraries and data files
#####################################################################
# Libraries
library(dplyr)
library(tidyverse)

# Read data files for gas exchange
phys <- read.csv("../data/raw_data/TT23_gasExchange.csv")
multispeq <- read.csv("../data/raw_data/TT23_multispeq_data.csv")

# Read resin strip data file, then subset by composite and standardize
# values by number of days deployed. Also calculate plant-available N
# concentration as sum of nitrate and ammonium concentration
nutrients <- read.csv("../data/raw_data/TT23_resin_strip_data.csv") %>%
  group_by(round, plot, composite, days_deployed) %>%
  summarize(phosphate_ppm = mean(phosphate_ppm, na.rm = TRUE),
            nitrate_ppm = mean(nitrate_ppm, na.rm = TRUE),
            ammonium_ppm = mean(ammonium_ppm, na.rm = TRUE)) %>%
  mutate(n_plantAvail = nitrate_ppm + ammonium_ppm,
         phosphate_ppm_day = phosphate_ppm / days_deployed,
         nitrate_ppm_day = nitrate_ppm / days_deployed,
         ammonium_ppm_day = ammonium_ppm / days_deployed,
         n_plantAvail_day = n_plantAvail / days_deployed,
         canopy = ifelse(round == 1, 
                         "pre_closure",
                         ifelse(round == 2,
                                "post_closure",
                                NA))) %>%
  ungroup(round) %>%
  dplyr::select(plot, composite, canopy, 
                days_deployed:n_plantAvail_day) %>%
  filter(plot == 3 | plot == 5 | plot == 7)

#####################################################################
# Join multispeq and physiology data set
#####################################################################
phys.total <- phys %>%
  full_join(multispeq) %>%
  full_join(nutrients, by = c("plot", "composite", "canopy")) %>%
  filter(!is.na(id))


#####################################################################
# Write compiled data file
#####################################################################
write.csv(phys.total, "../data/TT23_compiled_datasheet.csv", row.names = FALSE)



