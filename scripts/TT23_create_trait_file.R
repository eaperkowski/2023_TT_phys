#####################################################################
# Load libraries and data files
#####################################################################
# Libraries
library(dplyr)
library(tidyverse)

# Data files
phys <- read.csv("../data/TT23_phys_data_long.csv")
multispeq <- read.csv("../data/TT23_multispeq_data.csv") %>%
  dplyr::select(id:canopy, FmPrime:FvP_over_FmP, NPQt, 
                Phi2:PhiNPQ, SPAD, -quad)
nutrients <- read.csv("../data/raw_data/TT23_resin_strip_data.csv")

#####################################################################
# Join multispeq and physiology data set
#####################################################################
phys.total <- phys %>%
  full_join(multispeq) %>%
  mutate(composite = 
           ifelse(subplot == 1 | subplot == 2 | subplot == 7 | subplot == 8,
                  "C1",
                  ifelse(subplot == 3 | subplot == 9,
                         "C2",
                         ifelse(subplot == 4 | subplot == 10, 
                                "C3", 
                                ifelse(subplot == 5 | subplot == 6 | subplot == 11 | subplot == 12,
                                       "C4",
                                       ifelse(subplot == 13 | subplot == 14 | subplot == 19 | subplot == 20,
                                              "C5", 
                                              ifelse(subplot == 15 | subplot == 21,
                                                     "C6", 
                                                     ifelse(subplot == 16 | subplot == 22,
                                                            "C7",
                                                            ifelse(subplot == 17 | subplot == 18 | subplot == 23 | subplot == 24,
                                                                   "C8",
                                                                   ifelse(subplot == 25 | subplot == 26 | subplot == 31 | subplot == 32,
                                                                          "C9", 
                                                                          ifelse(subplot == 27 | subplot == 33,
                                                                                 "C10",
                                                                                 ifelse(subplot == 28 | subplot == 34,
                                                                                        "C11",
                                                                                        ifelse(subplot == 29 | subplot == 30 | subplot == 35 | subplot == 36,
                                                                                               "C12", NA)))))))))))))
