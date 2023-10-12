# Compiling script for CO2 response curves and dark respiration values
# Note that all paths in the first part of this script assume working directory 
# root is the path to this script in the "trillium_trail" git
# repository. Path assigned via Apple operating system, may differ  on Windows 
# operating systems.

###############################################################################
## Load readLicorData package
###############################################################################
library(readLicorData)
library(dplyr)

###############################################################################
## Load co2 response curve data files for initial field campaign
###############################################################################

# Albert (note: using a modified .txt file where leaf area was originally 
# modified from Zinny's experiment. .txt file corrected to 6 cm2, so no need 
# to change calculation here. Albert data also filtered for initial curves 
# where gasket issues were apparent)
alb_0419 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-19_co2resp_albert.txt") %>%
  filter(id != "5488" & id != "5488_light")
# write.csv(alb_0419, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-19_co2resp_albert.csv",
#           row.names = FALSE)

# Ozzie (same leaf area issues reported for albert)
ozz_0419 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-19_co2resp_ozz.txt")
# write.csv(ozz_0419, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-19_co2resp_ozz.csv",
#           row.names = FALSE)

# Stan
stan_0419 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-19_co2resp_stan")
# write.csv(stan_0419, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-19_co2resp_stan.csv",
#           row.names = FALSE)

stan_0420 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-20_co2resp_stan")
# write.csv(stan_0420, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-20_co2resp_stan.csv",
#           row.names = FALSE)

stan_0421 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-21_co2resp_stan")
# write.csv(stan_0421, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-21_co2resp_stan.csv",
#           row.names = FALSE)

stan_0505 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-05-05_co2resp_stan")
# write.csv(stan_0505, "../licor_data/licor_cleaned/pre_canopy_closure/2023-05-05_co2resp_stan.csv",
          # row.names = FALSE)

stan_0506 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-05-06_co2resp_stan")
write.csv(stan_0506, "../licor_data/licor_cleaned/pre_canopy_closure/2023-05-06_co2resp_stan.csv",
           row.names = FALSE)

# Yadi
yadi_0419 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-19_co2resp_yadi")
# write.csv(yadi_0419, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-19_co2resp_yadi.csv",
#           row.names = FALSE)

yadi_0420 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-20_co2resp_yadi")
# write.csv(yadi_0420, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-20_co2resp_yadi.csv",
#           row.names = FALSE)

yadi_0421 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-04-21_co2resp_yadi")
# write.csv(yadi_0421, "../licor_data/licor_cleaned/pre_canopy_closure/2023-04-21_co2resp_yadi.csv",
#           row.names = FALSE)

yadi_0505_a <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-05-05_co2resp_yadi")
write.csv(yadi_0505_a, "../licor_data/licor_cleaned/pre_canopy_closure/2023-05-05_co2resp_yadi_a.csv",
          row.names = FALSE)

yadi_0505_b <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-05-05_co2resp_yadi_b")
write.csv(yadi_0505_b, "../licor_data/licor_cleaned/pre_canopy_closure/2023-05-05_co2resp_yadi_b.csv",
          row.names = FALSE)

yadi_0506 <- licorData("../licor_data/licor_raw/pre_canopy_closure/2023-05-06_co2resp_yadi")
write.csv(yadi_0506, "../licor_data/licor_cleaned/pre_canopy_closure/2023-05-06_co2resp_yadi.csv",
          row.names = FALSE)

###############################################################################
## Curve-fitting prep
###############################################################################
file.list <- list.files(path = "../licor_data/licor_cleaned/pre_canopy_closure",
                             recursive = TRUE, pattern = "\\.csv$",
                             full.names = TRUE)

file.list <- setNames(file.list, stringr::str_extract(basename(file.list), 
                                                      '.*(?=\\.csv)'))

# Merge list of data frames, arrange by machine, measurement type, id, and time elapsed
merged_curves <- lapply(file.list, read.csv) %>%
  reshape::merge_all() %>%
  arrange(machine, id, elapsed)

# Write .csv file that compiles all licor data for Trillium and Maianthemum
# collected before canopy closure
write.csv(merged_curves, "../data/TT23_licor_merged_pre_closure.csv",
          row.names = FALSE)

###############################################################################
## Load co2 response curve data files for second field campaign
###############################################################################
stan_0612 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-12_co2resp_stan.txt")
# write.csv(stan_0612, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-12_co2resp_stan.csv",
#           row.names = FALSE)

yadi_0612 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-12_co2resp_yadi")
# write.csv(yadi_0612, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-12_co2resp_yadi.csv",
#           row.names = FALSE)

alb_0612 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-12_co2resp_albert")
# write.csv(alb_0612, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-12_co2resp_albert.csv",
#           row.names = FALSE)

stan_0613 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-13_co2resp_stan")
# write.csv(stan_0613, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-13_co2resp_stan.csv",
#           row.names = FALSE)

yadi_0613 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-13_co2resp_yadi")
# write.csv(yadi_0613, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-13_co2resp_yadi.csv",
#           row.names = FALSE)

alb_0613 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-13_co2resp_albert")
# write.csv(alb_0613, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-13_co2resp_albert.csv",
#           row.names = FALSE)

stan_0614 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-14_co2resp_stan")
# write.csv(stan_0614, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-14_co2resp_stan.csv",
#           row.names = FALSE)

yadi_0614 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-14_co2resp_yadi")
# write.csv(yadi_0614, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-14_co2resp_yadi.csv",
#           row.names = FALSE)

alb_0614 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-14_co2resp_albert")
# write.csv(alb_0614, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-14_co2resp_albert.csv",
#           row.names = FALSE)

stan_0615 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-15_co2resp_stan")
# write.csv(stan_0615, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-15_co2resp_stan.csv",
#           row.names = FALSE)

yadi_0615 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-15_co2resp_yadi")
# write.csv(yadi_0615, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-15_co2resp_yadi.csv",
#           row.names = FALSE)

alb_0615 <- licorData("../licor_data/licor_raw/post_canopy_closure/2023-06-15_co2resp_albert")
# write.csv(alb_0615, "../licor_data/licor_cleaned/post_canopy_closure/2023-06-15_co2resp_albert.csv",
#           row.names = FALSE)

###############################################################################
## Curve-fitting prep
###############################################################################
file.list2 <- list.files(path = "../licor_data/licor_cleaned/post_canopy_closure",
                        recursive = TRUE, pattern = "\\.csv$",
                        full.names = TRUE)

file.list2 <- setNames(file.list, stringr::str_extract(basename(file.list), 
                                                      '.*(?=\\.csv)'))

# Merge list of data frames, arrange by machine, measurement type, id, and time elapsed
merged_curves2 <- lapply(file.list, read.csv) %>%
  reshape::merge_all() %>%
  mutate(date = lubridate::ymd_hms(date)) %>%
  arrange(machine, date, id)

# Write .csv file that compiles all licor data for Trillium and Maianthemum
# collected before canopy closure
write.csv(merged_curves2, "../data/TT23_licor_merged_post_closure.csv",
          row.names = FALSE)
