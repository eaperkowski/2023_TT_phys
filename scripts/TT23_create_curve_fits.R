###############################################################################
## Load libraries and relevant datasets
###############################################################################

## Load libraries
library(dplyr)
library(plantecophys)

# Load clean and compiled Licor data sheet
licor.data <- read.csv("../data/TT23_licor_merged_pre_closure.csv") %>%
  dplyr::select(id, machine, date, A, Ci, Ca, CO2_r, gsw, Tleaf, VPDleaf, Qin, 
         S, Flow, Tair, Tleaf, Fan_speed, machine) %>%
  mutate_at(c("A", "Ci", "Ca", "CO2_r", "gsw", "Tleaf", "VPDleaf", "Qin", 
              "S", "Flow", "Tair", "Tleaf", "Fan_speed"), as.numeric)

# Load species metadata
spp.meta <- read.csv("../log/TT_tag_metadata.csv") %>%
  dplyr::select(id = tag_id, everything()) %>%
  mutate(id = ifelse(id == "flag_2" & spp == "Mai", "flag2_mai", id),
         id = ifelse(id == "flag_3" & spp == "Mai", "flag3_mai", id),
         id = ifelse(id == "flag_3" & spp == "Tri", "flag3_tri", id),
         id = ifelse(id == "9412", "9412_actual", id),
         id = ifelse(id == "2131" & spp == "Tri", "2131_tri", id),
         id = ifelse(id == "2131" & spp == "Mai", "2131_mai", id))

spp.meta.closed <- read.csv("../log/TT_tag_metadata_canopy_closed.csv") %>%
  mutate(id = ifelse(id == "flag2" & spp == "Mai", "flag2_mai", id),
         id = ifelse(id == "flag2" & spp == "Tri", "flag2_tri", id),
         id = ifelse(id == "flag3" & spp == "Mai", "flag3_mai", id),
         id = ifelse(id == "flag3" & spp == "Tri", "flag3_tri", id))

# Load garlic mustard density data measured by Jessie
gm.density.quad <- read.csv("../data/TT23_GM_density_counts_quad.csv") %>%
  dplyr::select(plot:quad, rosette_quad = rosette_count, adult_quad = adult_count,
                total_quad = total_count)

gm.density.plot <- gm.density.quad %>%
  group_by(plot, subplot) %>%
  summarize(rosette_subplot = sum(rosette_quad), 
            adult_subplot = sum(adult_quad),
            total_subplot = sum(total_quad))

# Add keep row delimiter to remove points
licor.data$keep.row <- "yes"

# Load custom fxn for stomatal limitation of photosynthesis
source("../../r_functions/stomatal_limitation.R")

###############################################################################
###############################################################################
## PRE CANOPY CLOSURE FITS
###############################################################################
###############################################################################

###############################################################################
## Plot 3: Trillium
###############################################################################
## Tag: 6578
tri.6578.pre <- fitaci(data = subset(licor.data, id == "6578_b" & keep.row == "yes"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.6578.pre)
aci.coefs <- data.frame(id = "6578", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                        t(coef(tri.6578.pre)))

## Tag: 5488
tri.5488.pre <- fitaci(data = subset(licor.data, id == "5488"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.5488.pre)
aci.coefs[2,] <- c(id = "5488", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                        t(coef(tri.5488.pre)))

## Tag: 4714
tri.4714.pre <- fitaci(data = subset(licor.data, id == "4714_a" & keep.row == "yes"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.4714.pre)
aci.coefs[3,] <- c(id = "4714", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.4714.pre)))

## Tag: 902
tri.902.pre <- fitaci(data = subset(licor.data, id == "902" & (A < 1 | A > 16)),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.902.pre)
aci.coefs[4,] <- c(id = "902", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.902.pre)))

## Tag: 3563
tri.3563.pre <- fitaci(data = subset(licor.data, id == "3563" & (A < 1.5 | A > 13)),
                  varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                   Ci = "Ci", PPFD = "Qin"))
plot(tri.3563.pre)
aci.coefs[5,] <- c(id = "3563", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.3563.pre)))

## Tag: 1686
tri.1686.pre <- fitaci(data = subset(licor.data, id == "1686" & (A < 1.25 | A > 15)),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.1686.pre)
aci.coefs[6,] <- c(id = "1686", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.1686.pre)))

## Tag: 425
# tri.425.pre <- fitaci(data = subset(licor.data, id == "425" & keep.row == "yes"),
#                    varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
#                                     Ci = "Ci", PPFD = "Qin"),
#                    fitTPU = TRUE)
# plot(tri.425.pre)
aci.coefs[7,] <- c(id = "425", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   Vcmax = NA, Jmax = NA, Rd = NA)

## Tag: 2416
tri.2416.pre <- fitaci(data = subset(licor.data, id == "2416" & (A < 1.25 | A > 14)),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.2416.pre)
aci.coefs[8,] <- c(id = "2416", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.2416.pre)))

## Tag: 552
tri.552.pre <- fitaci(data = subset(licor.data, id == "552" & keep.row == "yes"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.552.pre)
aci.coefs[9,] <- c(id = "552", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.552.pre)))

## Tag: 5495
tri.5495.pre <- fitaci(data = subset(licor.data, id == "5495" & keep.row == "yes"),
                  varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                   Ci = "Ci", PPFD = "Qin"))
# plot(tri.5495.pre)
aci.coefs[10,] <- c(id = "5495", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.5495.pre)))

## Tag: 916
tri.916.pre <- fitaci(data = subset(licor.data, id == "916" & keep.row == "yes"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.916.pre)
aci.coefs[11,] <- c(id = "916", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                   t(coef(tri.916.pre)))


## Tag: 583
tri.583.pre <- fitaci(data = subset(licor.data, id == "583a_redo" & (A < 2 | A > 11.3)),
                  varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                   Ci = "Ci", PPFD = "Qin"))
plot(tri.583.pre)
aci.coefs[12,] <- c(id = "583", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.583.pre)))

## Tag: 5479
tri.5479.pre <- fitaci(data = subset(licor.data, id == "5479a" & (A > 0 & A < 1 | A > 13.5)),
                  varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                   Ci = "Ci", PPFD = "Qin"))
# plot(tri.5479.pre)
aci.coefs[13,] <- c(id = "5479", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.5479.pre)))

## Tag: 4934
tri.4934.pre <- fitaci(data = subset(licor.data, id == "4934a_redo" & keep.row == "yes"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.4934.pre)
aci.coefs[14,] <- c(id = "4934", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.4934.pre)))

## Tag: 6558
tri.6558.pre <- fitaci(data = subset(licor.data, id == "6558" & (A < 1.5 | A > 16)),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.6558.pre)
aci.coefs[15,] <- c(id = "6558", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.6558.pre)))
## Tag: 2329
tri.2329.pre <- fitaci(data = subset(licor.data, id == "2329" & (A < 1.5 | A > 12.1)),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"),
                   citransition = 250)
# plot(tri.2329.pre)
aci.coefs[16,] <- c(id = "2329", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.2329.pre)))

## Tag: 774
tri.774.pre <- fitaci(data = subset(licor.data, id == "774" & (A < 1 | A > 14)),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"),
                   citransition = 300)
# plot(tri.774.pre)
aci.coefs[17,] <- c(id = "774", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.774.pre)))

## Tag: 4942
tri.4942.pre <- fitaci(data = subset(licor.data, id == "4942" & (A < 1.5 | A > 15)),
                  varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                   Ci = "Ci", PPFD = "Qin"))
# plot(tri.4942.pre)
aci.coefs[18,] <- c(id = "4942", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.4942.pre)))


## Tag: 4942
tri.1926.pre <- fitaci(data = subset(licor.data, id == "1926"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.1926.pre)
aci.coefs[136,] <- c(id = "1926", spp = "Tri", plot = "3", timepoint = "pre_closure", 
                    t(coef(tri.4942.pre)))


###############################################################################
## Plot 5: Trillium
###############################################################################

## Tag: 2547
tri.2547.pre <- fitaci(data = subset(licor.data, id == "2547b" & keep.row == "yes"),
                   varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                    Ci = "Ci", PPFD = "Qin"))
# plot(tri.2547.pre)
aci.coefs[19,] <- c(id = "2547", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.2547.pre)))

## Tag: 5115
tri.5115.pre <- fitaci(data = subset(licor.data, id == "5115" & (A < 1 | A > 15)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 300)
# plot(tri.5115.pre)
aci.coefs[20,] <- c(id = "5115", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.5115.pre)))

## Tag: 5228
tri.5228.pre <- fitaci(data = subset(licor.data, id == "5228" & (A > 0 & A < 1 | A > 10)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       fitmethod = 'bilinear', citransition = 250)
# plot(tri.5228.pre)
aci.coefs[21,] <- c(id = "5228", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.5228.pre)))

## Tag: 4000
tri.4000.pre <- fitaci(data = subset(licor.data, id == "4000" & (A > 0 & A < 1.5 | A > 18)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.4000.pre)
aci.coefs[22,] <- c(id = "4000", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4000.pre)))


## Tag: 3004
tri.3004.pre <- fitaci(data = subset(licor.data, id == "3004" & (A < 1 | A > 15)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 300)
# plot(tri.3004.pre)
aci.coefs[23,] <- c(id = "3004", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.3004.pre)))


## Tag: 3077
tri.3077.pre <- fitaci(data = subset(licor.data, id == "3077" & (A > 0 & A < 1 | A > 15)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 250)
# plot(tri.3077.pre)
aci.coefs[24,] <- c(id = "3077", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.3077.pre)))


## Tag: 4576
tri.4576.pre <- fitaci(data = subset(licor.data, id == "4576" & (A > 0 & A < 1 | A > 12.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 290)
# plot(tri.4576.pre)
aci.coefs[25,] <- c(id = "4576", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4576.pre)))


## Tag: 482
tri.482.pre <- fitaci(data = subset(licor.data, id == "482" & (A > 0 & A < 1 | A > 12)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                      citransition = 250)
# plot(tri.482.pre)
aci.coefs[26,] <- c(id = "482", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.482.pre)))


## Tag: 4414
tri.4414.pre <- fitaci(data = subset(licor.data, id == "4414" & (A > 0 & A < 1 | A > 15.7)),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
# plot(tri.4414.pre)
aci.coefs[27,] <- c(id = "4414", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4414.pre)))


## Tag: 4959
tri.4959.pre <- fitaci(data = subset(licor.data, id == "4959" & (A > 0 & A < 1 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 225)
# plot(tri.4959.pre)
aci.coefs[28,] <- c(id = "4959", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4959.pre)))


## Tag: 4990
tri.4990.pre <- fitaci(data = subset(licor.data, id == "4990" & (A > 0 & A < 1.5 | A > 15)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 252)
# plot(tri.4990.pre)
aci.coefs[29,] <- c(id = "4990", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4990.pre)))


## Tag: 4431
tri.4431.pre <- fitaci(data = subset(licor.data, id == "4431" & (A > 0 & A < 2 | A > 11)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.4431.pre)
aci.coefs[30,] <- c(id = "4431", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4431.pre)))


## Tag: 4630 (note: can't find in plot)
aci.coefs[31,] <- c(id = "4630", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    Vcmax = NA, Jmax = NA, Rd = NA)


## Tag: 5877
tri.5877.pre <- fitaci(data = subset(licor.data, id == "5877" & (A > 0 & A < 1 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 240)
# plot(tri.5877.pre)
aci.coefs[32,] <- c(id = "5877", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.5877.pre)))


## Tag: 2988
tri.2988.pre <- fitaci(data = subset(licor.data, id == "2988" & (A > 0 & A < 2 | A > 14.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 240)
# plot(tri.2988.pre)
aci.coefs[33,] <- c(id = "2988", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.2988.pre)))


## Tag: 4109
tri.4109.pre <- fitaci(data = subset(licor.data, id == "4109" & (A > 0 & A < 2 | A > 14.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.4109.pre)
aci.coefs[34,] <- c(id = "4109", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4109.pre)))


## Tag: 2508
tri.2508.pre <- fitaci(data = subset(licor.data, id == "2508" & (A > 0 & A < 2 | A > 14.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2508.pre)
aci.coefs[35,] <- c(id = "2508", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.2508.pre)))


## Tag: 43
tri.43.pre <- fitaci(data = subset(licor.data, id == "43"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.43.pre)
aci.coefs[137,] <- c(id = "43", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.43.pre)))


## Tag: 86
tri.86.pre <- fitaci(data = subset(licor.data, id == "86"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.86.pre)
aci.coefs[36,] <- c(id = "86", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.86.pre)))


## Tag: 5904
tri.5904.pre <- fitaci(data = subset(licor.data, id == "5904"),
                     varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                      Ci = "Ci", PPFD = "Qin"))
# plot(tri.5904.pre)
aci.coefs[37,] <- c(id = "5904", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.5904.pre)))


## Tag: 7147
tri.7147.pre <- fitaci(data = subset(licor.data, id == "7147"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.7147.pre)
aci.coefs[38,] <- c(id = "7147", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.7147.pre)))


## Tag: 1374
tri.1374.pre <- fitaci(data = subset(licor.data, id == "1374" & (A < 4 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.1374.pre)
aci.coefs[39,] <- c(id = "1374", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.1374.pre)))


## Tag: 4105
tri.4105.pre <- fitaci(data = subset(licor.data, id == "4105" & (A < 4 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 250)
# plot(tri.4105.pre)
aci.coefs[40,] <- c(id = "4105", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4105.pre)))


## Tag: 4547
tri.4547.pre <- fitaci(data = subset(licor.data, id == "4547" & (A < 4 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.4547.pre)
aci.coefs[41,] <- c(id = "4547", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4547.pre)))


## Tag: 5381
tri.5381.pre <- fitaci(data = subset(licor.data, id == "5381" & (A < 4 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.5381.pre)
aci.coefs[42,] <- c(id = "5381", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.5381.pre)))


## Tag: 614
tri.614.pre <- fitaci(data = subset(licor.data, id == "614" & (A < 4 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.614.pre)
aci.coefs[43,] <- c(id = "614", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.614.pre)))


## Tag: 2912
tri.2912.pre <- fitaci(data = subset(licor.data, id == "2912" & (A < 6 | A > 16)),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
# plot(tri.2912.pre)
aci.coefs[44,] <- c(id = "2912", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.2912.pre)))


## Tag: 1795
tri.1795.pre <- fitaci(data = subset(licor.data, id == "1795" & (A < 6 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.1795.pre)
aci.coefs[45,] <- c(id = "1795", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.1795.pre)))


## Tag: 2563
tri.2563.pre <- fitaci(data = subset(licor.data, id == "2563" & (A < 4 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2563.pre)
aci.coefs[46,] <- c(id = "2563", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.2563.pre)))


## Tag: flag_8
tri.flag_8.pre <- fitaci(data = subset(licor.data, id == "flag_8" & (A < 4 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.flag_8.pre)
aci.coefs[47,] <- c(id = "flag_8", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.flag_8.pre)))


## Tag: 2573
tri.2573.pre <- fitaci(data = subset(licor.data, id == "2573" & (A < 4 | A > 16)),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
# plot(tri.2573.pre)
aci.coefs[48,] <- c(id = "2573", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.2573.pre)))


## Tag: flag_3
tri.flag_3.pre <- fitaci(data = subset(licor.data, id == "flag3" & (A < 4 | A > 16)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.flag_3.pre)
aci.coefs[49,] <- c(id = "flag3_tri", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.flag_3.pre)))

## Tag: 4177
tri.4177.pre <- fitaci(data = subset(licor.data, id == "4177"),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(tri.4177.pre)
aci.coefs[138,] <- c(id = "4177", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                    t(coef(tri.4177.pre)))

## Tag: 1739
tri.1739.pre <- fitaci(data = subset(licor.data, id == "1739"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.1739.pre)
aci.coefs[139,] <- c(id = "1739", spp = "Tri", plot = "5", timepoint = "pre_closure", 
                     t(coef(tri.4177.pre)))

###############################################################################
## Plot 7: Trillium
###############################################################################

## Tag: 632
tri.632.pre <- fitaci(data = subset(licor.data, id == "632"),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
# plot(tri.632.pre)
aci.coefs[50,] <- c(id = "632", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.632.pre)))


## Tag: 1827
tri.1827.pre <- fitaci(data = subset(licor.data, id == "1827"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(tri.1827.pre)
aci.coefs[51,] <- c(id = "1827", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.1827.pre)))


## Tag: 1851
tri.1851.pre <- fitaci(data = subset(licor.data, id == "1851" & (A < 4 | A > 15)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.1851.pre)
aci.coefs[52,] <- c(id = "1851", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.1851.pre)))


## Tag: 3071
tri.3071.pre <- fitaci(data = subset(licor.data, id == "3071" & (A < 4 | A > 15.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.3071.pre)
aci.coefs[53,] <- c(id = "3071", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.3071.pre)))


## Tag: 2927
tri.2927.pre <- fitaci(data = subset(licor.data, id == "2927" & (A < 4 | A > 15.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2927.pre)
aci.coefs[54,] <- c(id = "2927", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2927.pre)))


## Tag: 2996
tri.2996.pre <- fitaci(data = subset(licor.data, id == "2996" & (A < 4 | A > 15.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2996.pre)
aci.coefs[55,] <- c(id = "2996", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2996.pre)))


## Tag: 972
tri.972.pre <- fitaci(data = subset(licor.data, id == "972" & (A < 4 | A > 13.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.972.pre)
aci.coefs[56,] <- c(id = "972", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.972.pre)))


## Tag: 1647
tri.1647.pre <- fitaci(data = subset(licor.data, id == "1647" & (A < 4 | A > 12)),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
# plot(tri.1647.pre)
aci.coefs[57,] <- c(id = "1647", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.1647.pre)))


## Tag: 1669
tri.1669.pre <- fitaci(data = subset(licor.data, id == "1669" & (A < 4 | A > 15)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.1669.pre)
aci.coefs[58,] <- c(id = "1669", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.1669.pre)))


## Tag: 2131 (two curves). Note: 2131 is tagged twice across plots (Tri in plot
## 7 and Mai in plot 3. Resolving by subsetting by machine to include Tri from
## plot 7 here)
tri.2131.pre <- fitaci(data = subset(licor.data, id == "2131" & is.na(machine) & (A < 6 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.2131.pre)
aci.coefs[59,] <- c(id = "2131_tri", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2131.pre)))


## Tag: 2289
tri.2289.pre <- fitaci(data = subset(licor.data, id == "2289" & (A < 4 | A > 16.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2289.pre)
aci.coefs[60,] <- c(id = "2289", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2289.pre)))


## Tag: 2379
tri.2379.pre <- fitaci(data = subset(licor.data, id == "2379" & (A < 4 | A > 16.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2379.pre)
aci.coefs[61,] <- c(id = "2379", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2379.pre)))


## Tag: 2382
tri.2382.pre <- fitaci(data = subset(licor.data, id == "2382" & (A < 4 | A > 12)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2382.pre)
aci.coefs[62,] <- c(id = "2382", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2382.pre)))


## Tag: 1997
tri.1997.pre <- fitaci(data = subset(licor.data, id == "1997" & (A < 4 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.1997.pre)
aci.coefs[63,] <- c(id = "1997", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.1997.pre)))


## Tag: 2886
tri.2886.pre <- fitaci(data = subset(licor.data, id == "2886" & (A < 4 | A > 11.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.2886.pre)
aci.coefs[64,] <- c(id = "2886", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.2886.pre)))


## Tag: 3043
tri.3043.pre <- fitaci(data = subset(licor.data, id == "3043" & (A < 4 | A > 11.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.3043.pre)
aci.coefs[65,] <- c(id = "3043", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.3043.pre)))


## Tag: 3427
tri.3427.pre <- fitaci(data = subset(licor.data, id == "3427" & (A < 3 | A > 10.5)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.3427.pre)
aci.coefs[66,] <- c(id = "3427", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.3427.pre)))


## Tag: 3526
tri.3526.pre <- fitaci(data = subset(licor.data, id == "3526" & (A < 3 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.3526.pre)
aci.coefs[67,] <- c(id = "3526", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.3526.pre)))


## Tag: 1975
tri.1975.pre <- fitaci(data = subset(licor.data, id == "1975" & (A < 2 | A > 12)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.1975.pre)
aci.coefs[68,] <- c(id = "1975", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.1975.pre)))


## Tag: 4081
tri.4081.pre <- fitaci(data = subset(licor.data, id == "4081" & (A < 3 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.4081.pre)
aci.coefs[69,] <- c(id = "4081", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.4081.pre)))


## Tag: 6507
tri.6507.pre <- fitaci(data = subset(licor.data, id == "6507" & (A < 10 | A > 14)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.6507.pre)
aci.coefs[70,] <- c(id = "6507", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.6507.pre)))


## Tag: 7574
tri.7574.pre <- fitaci(data = subset(licor.data, id == "7574"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
# plot(tri.7574.pre)
aci.coefs[71,] <- c(id = "7574", spp = "Tri", plot = "7", timepoint = "pre_closure", 
                    t(coef(tri.7574.pre)))


###############################################################################
## Plot 3: Maianthemum
###############################################################################

## Tag: 3643
mai.3643.pre <- fitaci(data = subset(licor.data, id == "3643"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.3643.pre)
aci.coefs[72,] <- c(id = "3643", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.3643.pre)))

## Tag: 141
mai.141.pre <- fitaci(data = subset(licor.data, id == "141"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.141.pre)
aci.coefs[73,] <- c(id = "141", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.141.pre)))


## Tag: 4250
mai.4250.pre <- fitaci(data = subset(licor.data, id == "4250"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.4250.pre)
aci.coefs[74,] <- c(id = "4250", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.4250.pre)))


## Tag: 9412_actual
mai.9412act.pre <- fitaci(data = subset(licor.data, id == "9412_actual"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.9412act.pre)
aci.coefs[75,] <- c(id = "9412_actual", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.9412act.pre)))

## Tag: 9412
mai.9412.pre <- fitaci(data = subset(licor.data, id == "9412"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.9412.pre)
aci.coefs[135,] <- c(id = "9412", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.9412.pre)))


## Tag: 2310
mai.2310.pre <- fitaci(data = subset(licor.data, id == "2310"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2310.pre)
aci.coefs[76,] <- c(id = "2310", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.2310.pre)))


## Tag: 1021
mai.1021.pre <- fitaci(data = subset(licor.data, id == "1021"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1021.pre)
aci.coefs[77,] <- c(id = "1021", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.1021.pre)))


## Tag: 2268
mai.2268.pre <- fitaci(data = subset(licor.data, id == "2268"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2268.pre)
aci.coefs[78,] <- c(id = "2268", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.2268.pre)))


## Tag: 5069
mai.5069.pre <- fitaci(data = subset(licor.data, id == "5069"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5069.pre)
aci.coefs[79,] <- c(id = "5069", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.5069.pre)))


## Tag: 5030
mai.5030.pre <- fitaci(data = subset(licor.data, id == "5030"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5030.pre)
aci.coefs[80,] <- c(id = "5030", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.5030.pre)))


## Tag: 2337
mai.2337.pre <- fitaci(data = subset(licor.data, id == "2337"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2337.pre)
aci.coefs[81,] <- c(id = "2337", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.2337.pre)))


## Tag: 179
mai.179.pre <- fitaci(data = subset(licor.data, id == "179"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.179.pre)
aci.coefs[82,] <- c(id = "179", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.179.pre)))

## Tag: 511
mai.511.pre <- fitaci(data = subset(licor.data, id == "511"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.511.pre)
aci.coefs[83,] <- c(id = "511", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.511.pre)))


## Tag: 543
mai.543.pre <- fitaci(data = subset(licor.data, id == "543"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.543.pre)
aci.coefs[84,] <- c(id = "543", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.543.pre)))


## Tag: 1157
mai.1157.pre <- fitaci(data = subset(licor.data, id == "1157"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.1157.pre)
aci.coefs[85,] <- c(id = "1157", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.1157.pre)))


## Tag: 1663
mai.1663.pre <- fitaci(data = subset(licor.data, id == "1663"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1663.pre)
aci.coefs[86,] <- c(id = "1663", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.1663.pre)))


## Tag: 2131 Note: 2131 is tagged twice across plots (Tri in plot
## 7 and Mai in plot 3. Resolving by subsetting by machine to include Mai from
## stan and plot 3 here)
mai.2131.pre <- fitaci(data = subset(licor.data, id == "2131" & machine == "stan"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2131.pre)
aci.coefs[87,] <- c(id = "2131_mai", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                   t(coef(mai.2131.pre)))

## Tag: 3665
mai.3665.pre <- fitaci(data = subset(licor.data, id == "3665"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.3665.pre)
aci.coefs[88,] <- c(id = "3665", spp = "Mai", plot = "3", timepoint = "pre_closure", 
                    t(coef(mai.3665.pre)))

###############################################################################
## Plot 5: Maianthemum
###############################################################################
## Tag: 5052
mai.5052.pre <- fitaci(data = subset(licor.data, id == "5052"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5052.pre)
aci.coefs[89,] <- c(id = "5052", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.5052.pre)))


## Tag: 2616
mai.2616.pre <- fitaci(data = subset(licor.data, id == "2616"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2616.pre)
aci.coefs[90,] <- c(id = "2616", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.2616.pre)))


## Tag: 2662
mai.2662.pre <- fitaci(data = subset(licor.data, id == "2662"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2662.pre)
aci.coefs[91,] <- c(id = "2662", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.2662.pre)))


## Tag: 4781
mai.4781.pre <- fitaci(data = subset(licor.data, id == "4781"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4781.pre)
aci.coefs[92,] <- c(id = "4781", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.4781.pre)))


## Tag: 5657
mai.5657.pre <- fitaci(data = subset(licor.data, id == "5657"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5657.pre)
aci.coefs[93,] <- c(id = "5657", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.5657.pre)))


## Tag: 631
mai.631.pre <- fitaci(data = subset(licor.data, id == "631"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.631.pre)
aci.coefs[94,] <- c(id = "631", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.631.pre)))


## Tag: 3960
mai.3960.pre <- fitaci(data = subset(licor.data, id == "3960"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.3960.pre)
aci.coefs[95,] <- c(id = "3960", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.3960.pre)))


## Tag: flag_2
mai.flag_2.pre <- fitaci(data = subset(licor.data, id == "flag 2"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.flag_2.pre)
aci.coefs[96,] <- c(id = "flag2_mai", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.flag_2.pre)))


## Tag: 1432
mai.1432.pre <- fitaci(data = subset(licor.data, id == "1432"),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(mai.1432.pre)
aci.coefs[97,] <- c(id = "1432", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.1432.pre)))


## Tag: 1499
mai.1499.pre <- fitaci(data = subset(licor.data, id == "1499"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1499.pre)
aci.coefs[98,] <- c(id = "1499", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.1499.pre)))


## Tag: 3993
mai.3993.pre <- fitaci(data = subset(licor.data, id == "3993"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.3993.pre)
aci.coefs[140,] <- c(id = "3993", spp = "Mai", plot = "5", timepoint = "pre_closure", 
                    t(coef(mai.3993.pre)))

###############################################################################
## Plot 7: Maianthemum
###############################################################################
## Tag: 1159
mai.1159.pre <- fitaci(data = subset(licor.data, id == "1159"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1159.pre)
aci.coefs[99,] <- c(id = "1159", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                    t(coef(mai.1159.pre)))


## Tag: 3539
mai.3539.pre <- fitaci(data = subset(licor.data, id == "3539"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.3539.pre)
aci.coefs[100,] <- c(id = "3539", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                    t(coef(mai.3539.pre)))


## Tag: 901
mai.901.pre <- fitaci(data = subset(licor.data, id == "901"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.901.pre)
aci.coefs[101,] <- c(id = "901", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.901.pre)))


## Tag: 1817
mai.1817.pre <- fitaci(data = subset(licor.data, id == "1817"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.1817.pre)
aci.coefs[102,] <- c(id = "1817", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.1817.pre)))


## Tag: 6027
mai.6027.pre <- fitaci(data = subset(licor.data, id == "6027"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.6027.pre)
aci.coefs[103,] <- c(id = "6027", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.6027.pre)))


## Tag: 1153
mai.1153.pre <- fitaci(data = subset(licor.data, id == "1153"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1153.pre)
aci.coefs[104,] <- c(id = "1153", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.1153.pre)))


## Tag: 6000
mai.6000.pre <- fitaci(data = subset(licor.data, id == "6000"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.6000.pre)
aci.coefs[105,] <- c(id = "6000", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.6000.pre)))


## Tag: 6808
mai.6808.pre <- fitaci(data = subset(licor.data, id == "6808"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.6808.pre)
aci.coefs[106,] <- c(id = "6808", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.6808.pre)))


## Tag: 545
mai.545.pre <- fitaci(data = subset(licor.data, id == "545"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.545.pre)
aci.coefs[107,] <- c(id = "545", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.545.pre)))


## Tag: 2546
mai.2546.pre <- fitaci(data = subset(licor.data, id == "2546" & (A < 1 | A > 7.5)),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.2546.pre)
aci.coefs[108,] <- c(id = "2546", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.2546.pre)))


## Tag: 235
mai.235.pre <- fitaci(data = subset(licor.data, id == "235" & (A < 1.5 | A > 10)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.235.pre)
aci.coefs[109,] <- c(id = "235", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.235.pre)))


## Tag: 6100
mai.6100.pre <- fitaci(data = subset(licor.data, id == "6100" & (A <3 | A > 10)),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.6100.pre)
aci.coefs[110,] <- c(id = "6100", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.6100.pre)))


## Tag: 6685
mai.6685.pre <- fitaci(data = subset(licor.data, id == "6685"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.6685.pre)
## Curve directly under maple seedling with leaves fully expanded. Dropped
## from analysis due to noisy curve and low A/gsw under ambient conditions
aci.coefs[111,] <- c(id = "6685", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     Vcmax = NA, Jmax = NA, Rd = NA)


## Tag: 4883
mai.4883.pre <- fitaci(data = subset(licor.data, id == "4883" & (A < 4 | A > 12)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4883.pre)
aci.coefs[112,] <- c(id = "4883", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4883.pre)))


## Tag: 4874
mai.4874.pre <- fitaci(data = subset(licor.data, id == "4874"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4874.pre)
aci.coefs[113,] <- c(id = "4874", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4874.pre)))


## Tag: 42
mai.42.pre <- fitaci(data = subset(licor.data, id == "42"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.42.pre)
aci.coefs[114,] <- c(id = "42", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.42.pre)))


## Tag: 4834
mai.4834.pre <- fitaci(data = subset(licor.data, id == "4834"),
                     varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                      Ci = "Ci", PPFD = "Qin"))
plot(mai.4834.pre)
aci.coefs[115,] <- c(id = "4834", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4834.pre)))


## Tag: 4844
mai.4844.pre <- fitaci(data = subset(licor.data, id == "4844"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4844.pre)
aci.coefs[116,] <- c(id = "4844", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4844.pre)))


## Tag: 4871
mai.4871.pre <- fitaci(data = subset(licor.data, id == "4871"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4871.pre)
aci.coefs[117,] <- c(id = "4871", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4871.pre)))


## Tag: 2666
mai.2666.pre <- fitaci(data = subset(licor.data, id == "2666"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2666.pre)
aci.coefs[118,] <- c(id = "2666", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.2666.pre)))


## Tag: 5073
mai.5073.pre <- fitaci(data = subset(licor.data, id == "5073"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5073.pre)
aci.coefs[119,] <- c(id = "5073", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.5073.pre)))


## Tag: 5015
mai.5015.pre <- fitaci(data = subset(licor.data, id == "5015"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5015.pre)
aci.coefs[120,] <- c(id = "5015", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.5015.pre)))


## Tag: 5077
mai.5077.pre <- fitaci(data = subset(licor.data, id == "5077"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.5077.pre)
aci.coefs[121,] <- c(id = "5077", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.5077.pre)))


## Tag: 1240
mai.1240.pre <- fitaci(data = subset(licor.data, id == "1240"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1240.pre)
aci.coefs[122,] <- c(id = "1240", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.1240.pre)))


## Tag: 1264
mai.1264.pre <- fitaci(data = subset(licor.data, id == "1264"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.1264.pre)
aci.coefs[123,] <- c(id = "1264", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.1264.pre)))


## Tag: 860
mai.860.pre <- fitaci(data = subset(licor.data, id == "860"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.860.pre)
aci.coefs[124,] <- c(id = "860", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.860.pre)))


## Tag: 682
mai.682.pre <- fitaci(data = subset(licor.data, id == "682"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.682.pre)
aci.coefs[125,] <- c(id = "682", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.682.pre)))


## Tag: 3581
mai.3581.pre <- fitaci(data = subset(licor.data, id == "3581"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.3581.pre)
aci.coefs[126,] <- c(id = "3581", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.3581.pre)))


## Tag: 2265
mai.2265.pre <- fitaci(data = subset(licor.data, id == "2265"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2265.pre)
aci.coefs[127,] <- c(id = "2265", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.2265.pre)))


## Tag: 2285
mai.2285.pre <- fitaci(data = subset(licor.data, id == "2285"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2285.pre)
aci.coefs[128,] <- c(id = "2285", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.2285.pre)))


## Tag: 4859
mai.4859.pre <- fitaci(data = subset(licor.data, id == "4859"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4859.pre)
aci.coefs[129,] <- c(id = "4859", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4859.pre)))


## Tag: 4736
mai.4736.pre <- fitaci(data = subset(licor.data, id == "4736"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4736.pre)
## Bad curve fit, removing from data
aci.coefs[130,] <- c(id = "4736", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     Vcmax = NA, Jmax = NA, Rd = NA)


## Tag: 4865
mai.4865.pre <- fitaci(data = subset(licor.data, id == "4865" & (A > 0 & (A < 2 | A > 10))),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 400)
plot(mai.4865.pre)
aci.coefs[131,] <- c(id = "4865", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4865.pre)))


## Tag: 4676
mai.4676.pre <- fitaci(data = subset(licor.data, id == "4676"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4676.pre)
aci.coefs[132,] <- c(id = "4676", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4676.pre)))


## Tag: 4802
mai.4802.pre <- fitaci(data = subset(licor.data, id == "4802" & A > 0 & (A < 1.5 | A > 10)),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"),
                       citransition = 375)
plot(mai.4802.pre)
aci.coefs[133,] <- c(id = "4802", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.4802.pre)))


## Tag: 652
mai.652.pre <- fitaci(data = subset(licor.data, id == "625"), ## note error in ID
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.652.pre)
aci.coefs[134,] <- c(id = "652", spp = "Mai", plot = "7", timepoint = "pre_closure", 
                     t(coef(mai.652.pre)))


###############################################################################
## Snapshot measures at start of each curve
###############################################################################
snapshot <- licor.data %>%
  mutate(id = ifelse(id == "2131" & is.na(machine), "2131_tri", id),
         id = ifelse(id == "2131" & machine == "stan", "2131_mai", id)) %>%
  group_by(id) %>%
  filter(row_number() == 1 & id != "392") %>%
  dplyr::select(id, machine, anet = A, ci = Ci, ca = Ca, co2_ref = CO2_r, gsw, Tleaf) %>%
  mutate(ci.ca = ci / ca,
         iwue = anet / gsw) %>%
  filter(id != "4714_b" & id!= "6578" & id != "6578_a" & id != "5479b" & 
           id != "5479c" & id != "583a" & id != "583b" & id != "4934b" & id != "2927_light" &
           id != "UserDefCon" & id != "") %>%
  mutate(machine = ifelse(is.na(machine), "yadi", machine),
         id = ifelse(id == "4714_a", "4714", id),
         id = ifelse(id == "6578_b", "6578", id),
         id = ifelse(id == "2547b", "2547", id),
         id = ifelse(id == "5479a", "5479", id),
         id = ifelse(id == "583a_redo", "583", id),
         id = ifelse(id == "flag 2", "flag2_mai", id),
         id = ifelse(id == "4934a_redo", "4934", id),
         id = ifelse(id == "flag3", "flag3_tri", id),
         id = ifelse(id == "625", "652", id))
aci.coefs$plot <- as.numeric(aci.coefs$plot)
spp.meta$plot <- as.numeric(spp.meta$plot)

###############################################################################
## Merge snapshot measures with curve coefficients, then add plot and subplot
## information 
###############################################################################
pre.canopy.merged <- snapshot %>%
  full_join(aci.coefs, by = "id") %>%
  mutate(plot = as.numeric(plot),
         jmax.vcmax = as.numeric(Vcmax) / as.numeric(Jmax),
         rd.vcmax = as.numeric(Rd) / as.numeric(Vcmax),
         stomatal_limitation(A_net = as.numeric(anet),
                             Vcmax = as.numeric(Vcmax),
                             Rd.meas = TRUE,
                             Rd = as.numeric(Rd),
                             leaf.temp = Tleaf,
                             temp = "C")[5]) %>%
  filter(id != "425" & id != 4630) %>%
  full_join(spp.meta, by = c("id", "plot")) %>%
  mutate(subplot = as.numeric(subplot),
         quad = as.numeric(quad)) %>%
  full_join(gm.density.plot, by = c("plot", "subplot")) %>%
  full_join(gm.density.quad, by = c("plot", "subplot", "quad")) %>%
  filter(!is.na(anet)) %>%
  mutate(rosette_subplot = ifelse(is.na(rosette_subplot), 0, rosette_subplot),
         adult_subplot = ifelse(is.na(adult_subplot), 0, adult_subplot),
         total_subplot = ifelse(is.na(total_subplot), 0, total_subplot),
         rosette_quad = ifelse(is.na(rosette_quad), 0, rosette_quad),
         adult_quad = ifelse(is.na(adult_quad), 0, rosette_quad),
         total_quad = ifelse(is.na(total_quad), 0, rosette_quad)) %>%
  dplyr::select(id, spp = spp.x, plot, subplot, quad, gm.trt, rosette_subplot:total_quad,
                canopy = timepoint, anet, gsw, ci.ca, 
                iwue, stom.lim = l, vcmax25 = Vcmax, 
                jmax25 = Jmax, rd25 = Rd, jmax.vcmax, 
                rd.vcmax) %>%
  filter(id != "9412") %>%
  mutate(id = ifelse(id == "9412_actual", "9412", id))

write.csv(pre.canopy.merged, "../data/TT23_pre_closure_phys_data.csv",
          row.names = FALSE)

###############################################################################
###############################################################################
## POST CANOPY CLOSURE FITS
###############################################################################
###############################################################################
# Load clean and compiled Licor data sheet
licor.data.post <- read.csv("../data/TT23_licor_merged_post_closure.csv") %>%
  dplyr::select(id, machine, date, A, Ci, Ca, CO2_r, gsw, Tleaf, VPDleaf, Qin, 
                S, Flow, Tair, Tleaf, Fan_speed, machine) %>%
  mutate_at(c("A", "Ci", "Ca", "CO2_r", "gsw", "Tleaf", "VPDleaf", "Qin", 
              "S", "Flow", "Tair", "Tleaf", "Fan_speed"), as.numeric) %>%
  mutate(date = lubridate::ymd_hms(date),
         day = lubridate::date(date))

###############################################################################
## Plot 5 - 06/12/2023
###############################################################################
## Tag: 2671
ari.2671.post <- fitaci(data = subset(licor.data.post, id == "2671"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(ari.2671.post)
aci.coefs.post <- data.frame(id = "2671", spp = "Ari", plot = "5", 
                             timepoint = "post_closure", 
                             t(coef(ari.2671.post)))

## Tag: 5115
tri.5115.post <- fitaci(data = subset(licor.data.post, id == "5115"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.5115.post)
aci.coefs.post[2,] <- c(id = "5115", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.5115.post)))

## Tag: 5337
ari.5337.post <- fitaci(data = subset(licor.data.post, id == "5337"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.5337.post)
aci.coefs.post[3,] <- c(id = "5337", spp = "Ari", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(ari.5337.post)))

## Tag: 2547
tri.2547.post <- fitaci(data = subset(licor.data.post, id == "2547"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.2547.post)
aci.coefs.post[4,] <- c(id = "2547", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.2547.post)))

## Tag: 2508
tri.2508.post <- fitaci(data = subset(licor.data.post, id == "2508"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.2508.post)
aci.coefs.post[5,] <- c(id = "2508", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.2508.post)))

## Tag: 3004
tri.3004.post <- fitaci(data = subset(licor.data.post, id == "3004"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.3004.post)
aci.coefs.post[6,] <- c(id = "3004", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.3004.post)))

## Tag: 4990
tri.4990.post <- fitaci(data = subset(licor.data.post, id == "4990"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.4990.post)
aci.coefs.post[7,] <- c(id = "4990", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.4990.post)))

## Tag: 482
tri.482.post <- fitaci(data = subset(licor.data.post, id == "482"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.482.post)
aci.coefs.post[8,] <- c(id = "482", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.482.post)))

## Tag: 4959
tri.4959.post <- fitaci(data = subset(licor.data.post, id == "4959"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.4959.post)
aci.coefs.post[9,] <- c(id = "4959", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.4959.post)))

## Tag: 4414
tri.4414.post <- fitaci(data = subset(licor.data.post, id == "4414"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.4414.post)
aci.coefs.post[10,] <- c(id = "4414", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.4414.post)))

## Tag: 4431
tri.4431.post <- fitaci(data = subset(licor.data.post, id == "4431"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.4431.post)
aci.coefs.post[11,] <- c(id = "4431", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.4431.post)))

## Tag: 3829
tri.3829.post <- fitaci(data = subset(licor.data.post, id == "3829"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.3829.post)
aci.coefs.post[12,] <- c(id = "3829", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.3829.post)))

## Tag: 5877
tri.5877.post <- fitaci(data = subset(licor.data.post, id == "5877"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.5877.post)
aci.coefs.post[13,] <- c(id = "5877", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.5877.post)))

## Tag: 2988
tri.2988.post <- fitaci(data = subset(licor.data.post, id == "2988"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.2988.post)
aci.coefs.post[14,] <- c(id = "2988", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.2988.post)))

## Tag: 5381
tri.5381.post <- fitaci(data = subset(licor.data.post, id == "5381"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.5381.post)
aci.coefs.post[15,] <- c(id = "5381", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.5381.post)))

## Tag: 4547
tri.4547.post <- fitaci(data = subset(licor.data.post, id == "4547"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.4547.post)
aci.coefs.post[16,] <- c(id = "4547", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.4547.post)))

## Tag: 2563
tri.2563.post <- fitaci(data = subset(licor.data.post, id == "2563"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.2563.post)
aci.coefs.post[17,] <- c(id = "2563", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.2563.post)))

## Tag: flag 8 (polka dot)
tri.flag8.post <- fitaci(data = subset(licor.data.post, id == "polka"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.flag8.post)
aci.coefs.post[18,] <- c(id = "flag8", spp = "Tri", plot = "5", 
                        timepoint = "post_closure", 
                        t(coef(tri.flag8.post)))

## Tag: flag 2
tri.flag2.post <- fitaci(data = subset(licor.data.post, id == "flag2" & is.na(machine)),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(tri.flag2.post)
aci.coefs.post[19,] <- c(id = "flag2_tri", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.flag2.post)))

## Tag: flag 3
tri.flag3.post <- fitaci(data = subset(licor.data.post, id == "flag3" & is.na(machine)),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(tri.flag3.post)
aci.coefs.post[20,] <- c(id = "flag3_tri", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.flag3.post)))

## Tag: 4177
tri.4177.post <- fitaci(data = subset(licor.data.post, id == "4177" & is.na(machine)),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(tri.4177.post)
aci.coefs.post[21,] <- c(id = "4177", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.4177.post)))

## Tag: 7147
tri.7147.post <- fitaci(data = subset(licor.data.post, id == "7147" & is.na(machine)),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.7147.post)
aci.coefs.post[22,] <- c(id = "7147", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.7147.post)))

## Tag: 43
licor.data.post$id[c(7064:7159)] <- 2912

tri.43.post <- fitaci(data = subset(licor.data.post, id == "43" & is.na(machine)),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.43.post)
aci.coefs.post[23,] <- c(id = "43", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.43.post)))

## Tag: 2912
tri.2912.post <- fitaci(data = subset(licor.data.post, id == "2912" & is.na(machine)),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(tri.2912.post)
aci.coefs.post[24,] <- c(id = "2912", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.2912.post)))

## Tag: 1374
tri.1374.post <- fitaci(data = subset(licor.data.post, id == "1374" & is.na(machine)),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.1374.post)
aci.coefs.post[25,] <- c(id = "1374", spp = "Tri", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(tri.1374.post)))

###############################################################################
## Plot 5 - 06/15/2023
###############################################################################
## Tag: 2662
mai.2662.post <- fitaci(data = subset(licor.data.post, id == "2662" & A > 0),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))

plot(mai.2662.post)
aci.coefs.post[26,] <- c(id = "2662", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.2662.post)))

## Tag: 5657
mai.5657.post <- fitaci(data = subset(licor.data.post, id == "5657"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.5657.post)
aci.coefs.post[27,] <- c(id = "5657", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.5637.post)))

## Tag: 6867
licor.data.post$id[c(3128:3223)] <- 5052


ari.6867.post <- fitaci(data = subset(licor.data.post, id == "6867"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.6867.post)
aci.coefs.post[28,] <- c(id = "6867", spp = "Ari", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(ari.6867.post)))

## Tag: 5052
mai.5052.post <- fitaci(data = subset(licor.data.post, id == "5052"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.5052.post)
aci.coefs.post[29,] <- c(id = "5052", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.5052.post)))

## Tag: 631
mai.631.post <- fitaci(data = subset(licor.data.post, id == "631"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.631.post)
aci.coefs.post[30,] <- c(id = "631", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.631.post)))

## Tag: 3993
mai.3993.post <- fitaci(data = subset(licor.data.post, id == "3993"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.3993.post)
aci.coefs.post[31,] <- c(id = "3993", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.3993.post)))

## Tag: flag 2 (Mai)
mai.flag2.post <- fitaci(data = subset(licor.data.post, id == "flag2" & machine == "stan"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.flag2.post)
aci.coefs.post[32,] <- c(id = "flag2_mai", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.flag2.post)))

## Tag: 2215
ari.2215.post <- fitaci(data = subset(licor.data.post, id == "2215"),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(ari.2215.post)
aci.coefs.post[33,] <- c(id = "2215", spp = "Ari", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(ari.2215.post)))

## Tag: 2616
mai.2616.post <- fitaci(data = subset(licor.data.post, id == "2616"),
                         varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                          Ci = "Ci", PPFD = "Qin"))
plot(mai.2616.post)
aci.coefs.post[34,] <- c(id = "2616", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.2616.post)))

## Tag: 1476
mai.1476.post <- fitaci(data = subset(licor.data.post, id == "1476"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"),
                        fitmethod = "bilinear")
plot(mai.1476.post)
aci.coefs.post[35,] <- c(id = "1476", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.1476.post)))

## Tag: 4781
mai.4781.post <- fitaci(data = subset(licor.data.post, id == "4781"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.4781.post)
aci.coefs.post[36,] <- c(id = "4781", spp = "Mai", plot = "5", 
                         timepoint = "post_closure", 
                         t(coef(mai.4781.post)))

###############################################################################
## Plot 7 - 06/13/2023
###############################################################################

## Tag: 901
mai.901.post <- fitaci(data = subset(licor.data.post, id == "901"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.901.post)
aci.coefs.post[37,] <- c(id = "901", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.901.post)))

## Tag: 5151
ari.5151.post <- fitaci(data = subset(licor.data.post, id == "5151"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(ari.5151.post)
aci.coefs.post[38,] <- c(id = "5151", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.5151.post)))

## Tag: 275
mai.275.post <- fitaci(data = subset(licor.data.post, id == "275"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.275.post)
aci.coefs.post[39,] <- c(id = "275", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.275.post)))

## Tag: 2229
mai.2229.post <- fitaci(data = subset(licor.data.post, id == "2229"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2229.post)
aci.coefs.post[40,] <- c(id = "2229", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.2229.post)))

## Tag: 1159
mai.1159.post <- fitaci(data = subset(licor.data.post, id == "1159"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.1159.post)
aci.coefs.post[41,] <- c(id = "1159", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.1159.post)))

## Tag: 6808
mai.6808.post <- fitaci(data = subset(licor.data.post, id == "6808"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.6808.post)
aci.coefs.post[42,] <- c(id = "6808", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.6808.post)))

## Tag: 2546
mai.2546.post <- fitaci(data = subset(licor.data.post, id == "2546"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.2546.post)
aci.coefs.post[43,] <- c(id = "2546", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.2546.post)))

## Tag: 4632
mai.4632.post <- fitaci(data = subset(licor.data.post, id == "4632"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.4632.post)
aci.coefs.post[44,] <- c(id = "4632", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.4632.post)))

## Tag: 3427 (missing 7 on file label)
tri.3427.post <- fitaci(data = subset(licor.data.post, id == "342"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.3427.post)
aci.coefs.post[45,] <- c(id = "3427", spp = "Tri", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(tri.3427.post)))

## Tag: 3526
tri.3526.post <- fitaci(data = subset(licor.data.post, id == "3526"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.3526.post)
aci.coefs.post[46,] <- c(id = "3526", spp = "Tri", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(tri.3526.post)))

## Tag: 42
mai.42.post <- fitaci(data = subset(licor.data.post, id == "42"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.42.post)
aci.coefs.post[47,] <- c(id = "42", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.42.post)))

## Tag: 652
mai.652.post <- fitaci(data = subset(licor.data.post, id == "652_b"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(mai.652.post)
aci.coefs.post[48,] <- c(id = "652", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.652.post)))

## Tag: 1997
tri.1997.post <- fitaci(data = subset(licor.data.post, id == "1997"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.1997.post)
aci.coefs.post[49,] <- c(id = "1997", spp = "Tri", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(tri.1997.post)))

## Tag: 5404
ari.5404.post <- fitaci(data = subset(licor.data.post, id == "5404"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.5404.post)
aci.coefs.post[50,] <- c(id = "5404", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.5404.post)))

## Tag: 5176
ari.5176.post <- fitaci(data = subset(licor.data.post, id == "5176"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"),
                        fitmethod = "bilinear")
plot(ari.5176.post)
aci.coefs.post[51,] <- c(id = "5176", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.5176.post)))

## Tag: 6890
ari.6890.post <- fitaci(data = subset(licor.data.post, id == "6890" & (A > 2 | A < 0)),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.6890.post)
aci.coefs.post[52,] <- c(id = "6890", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.6890.post)))

## Tag: 5253
ari.5253.post <- fitaci(data = subset(licor.data.post, id == "5253"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.5253.post)
aci.coefs.post[53,] <- c(id = "5253", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.5253.post)))

## Tag: 2438
ari.2438.post <- fitaci(data = subset(licor.data.post, id == "2438"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.2438.post)
aci.coefs.post[54,] <- c(id = "2438", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.2438.post)))

## Tag: 2996
tri.2996.post <- fitaci(data = subset(licor.data.post, id == "2996"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.2996.post)
aci.coefs.post[55,] <- c(id = "2996", spp = "Tri", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(tri.2996.post)))

## Tag: 1851
tri.1851.post <- fitaci(data = subset(licor.data.post, id == "1851"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.1851.post)
aci.coefs.post[56,] <- c(id = "1851", spp = "Tri", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(tri.1851.post)))

## Tag: 3071
tri.3071.post <- fitaci(data = subset(licor.data.post, id == "3071"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.3071.post)
aci.coefs.post[57,] <- c(id = "3071", spp = "Tri", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(tri.3071.post)))

## Tag: 784 (note: included as 1704 in data file)
ari.784.post <- fitaci(data = subset(licor.data.post, id == "1704"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.784.post)
aci.coefs.post[58,] <- c(id = "784", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.784.post)))

## Tag: 3785
ari.3785.post <- fitaci(data = subset(licor.data.post, id == "3785"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.3785.post)
aci.coefs.post[59,] <- c(id = "3785", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.3785.post)))

## Tag: 421
licor.data.post$id[c(5336:5431)] <- "421_redo"

ari.421.post <- fitaci(data = subset(licor.data.post, id == "421_redo"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(ari.421.post)
aci.coefs.post[60,] <- c(id = "421", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.421.post)))

## Tag: 818
ari.818.post <- fitaci(data = subset(licor.data.post, id == "818"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(ari.818.post)
aci.coefs.post[61,] <- c(id = "818", spp = "Ari", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(ari.818.post)))

## Tag: 4844
mai.4844.post <- fitaci(data = subset(licor.data.post, id == "4844-redo"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.4844.post)
aci.coefs.post[62,] <- c(id = "4844", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.4844.post)))

## Tag: 4871
mai.4871.post <- fitaci(data = subset(licor.data.post, id == "4871"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.4871.post)
aci.coefs.post[63,] <- c(id = "4871", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.4871.post)))

## Tag: 2666
mai.2666.post <- fitaci(data = subset(licor.data.post, id == "2666"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.2666.post)
aci.coefs.post[64,] <- c(id = "2666", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.2666.post)))

## Tag: 5077
mai.5077.post <- fitaci(data = subset(licor.data.post, id == "5077"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.5077.post)
aci.coefs.post[65,] <- c(id = "5077", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.5077.post)))

## Tag: 4874
mai.4874.post <- fitaci(data = subset(licor.data.post, id == "4874"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.4874.post)
aci.coefs.post[66,] <- c(id = "4874", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.4874.post)))



## Tag: 682 (2302 label used for both 682 and 2265)
licor.data.post$id[c(7801:7896)] <- "682"
licor.data.post$id[c(7897:7992)] <- "2265"

mai.682.post <- fitaci(data = subset(licor.data.post, id == "682"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.682.post)
aci.coefs.post[67,] <- c(id = "682", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.682.post)))

## Tag: 2265
mai.2265.post <- fitaci(data = subset(licor.data.post, id == "2265"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2265.post)
aci.coefs.post[68,] <- c(id = "2265", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.2265.post)))

## Tag: 4859
mai.4859.post <- fitaci(data = subset(licor.data.post, id == "4859"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.4859.post)
aci.coefs.post[69,] <- c(id = "4859", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.4859.post)))

## Tag: 3581
mai.3581.post <- fitaci(data = subset(licor.data.post, id == "3581"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.3581.post)
aci.coefs.post[70,] <- c(id = "3581", spp = "Mai", plot = "7", 
                         timepoint = "post_closure", 
                         t(coef(mai.3581.post)))

###############################################################################
## Plot 3 - 06/14/2023
###############################################################################
## Tag: 5069
mai.5069.post <- fitaci(data = subset(licor.data.post, id == "5069_b"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.5069.post)
aci.coefs.post[71,] <- c(id = "5069", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.5069.post)))

## Tag: 1021
mai.1021.post <- fitaci(data = subset(licor.data.post, id == "1021"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.1021.post)
aci.coefs.post[72,] <- c(id = "1021", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.1021.post)))

## Tag: 5030
mai.5030.post <- fitaci(data = subset(licor.data.post, id == "5030"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.5030.post)
aci.coefs.post[73,] <- c(id = "5030", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.5030.post)))

## Tag: 2268
mai.2268.post <- fitaci(data = subset(licor.data.post, id == "2268"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.2268.post)
aci.coefs.post[74,] <- c(id = "2268", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.2268.post)))

## Tag: 2416
tri.2416.post <- fitaci(data = subset(licor.data.post, id == "2416"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.2416.post)
aci.coefs.post[75,] <- c(id = "2416", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.2416.post)))

## Tag: 54
mai.54.post <- fitaci(data = subset(licor.data.post, id == "54"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.54.post)
aci.coefs.post[76,] <- c(id = "54", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.54.post)))

## Tag: Untagged_p3_sub17_q1_tri
tri.untagged1.post <- fitaci(data = subset(licor.data.post, id == "untagged_3_sub17_q1_tri"),
                      varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                       Ci = "Ci", PPFD = "Qin"))
plot(tri.untagged1.post)
aci.coefs.post[77,] <- c(id = "untagged_p3_sub17_q1", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.untagged1.post)))

## Tag: Untagged_p3_sub17_q1_tri
tri.untagged2.post <- fitaci(data = subset(licor.data.post, id == "untagged_3_sub12_q4_tri"),
                             varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                              Ci = "Ci", PPFD = "Qin"))
plot(tri.untagged2.post)
aci.coefs.post[78,] <- c(id = "untagged_p3_sub12_q4", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.untagged2.post)))

## Tag: 5488
tri.5488.post <- fitaci(data = subset(licor.data.post, id == "548"),
                             varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                              Ci = "Ci", PPFD = "Qin"))
plot(tri.5488.post)
aci.coefs.post[79,] <- c(id = "5488", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.5488.post)))

## Tag: 4250
mai.4250.post <- fitaci(data = subset(licor.data.post, id == "4250"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.4250.post)
aci.coefs.post[80,] <- c(id = "4250", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.4250.post)))

## Tag: 5184
mai.5184.post <- fitaci(data = subset(licor.data.post, id == "5184"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.5184.post)
aci.coefs.post[81,] <- c(id = "5184", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.5184.post)))

## Tag: 141
mai.141.post <- fitaci(data = subset(licor.data.post, id == "141"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.141.post)
aci.coefs.post[82,] <- c(id = "141", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.141.post)))

## Tag: 774
tri.774.post <- fitaci(data = subset(licor.data.post, id == "774"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.774.post)
aci.coefs.post[83,] <- c(id = "774", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.774.post)))

## Tag: 5479
tri.5479.post <- fitaci(data = subset(licor.data.post, id == "5479"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(tri.5479.post)
aci.coefs.post[84,] <- c(id = "5479", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.5479.post)))

## Tag: 3643
mai.3643.post <- fitaci(data = subset(licor.data.post, id == "3643"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.3643.post)
aci.coefs.post[85,] <- c(id = "3643", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.3643.post)))

## Tag: 179
mai.179.post <- fitaci(data = subset(licor.data.post, id == "179"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(mai.179.post)
aci.coefs.post[86,] <- c(id = "179", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.179.post)))

## Tag: 2131
mai.2131.post <- fitaci(data = subset(licor.data.post, id == "2131"),
                       varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                        Ci = "Ci", PPFD = "Qin"))
plot(mai.2131.post)
aci.coefs.post[87,] <- c(id = "2131_mai", spp = "Mai", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(mai.2131.post)))

## Tag: Unflagged trillium
tri.unflagged3.post <- fitaci(data = subset(licor.data.post, id == "unflagged_trillium_p3_sub22_q4"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.unflagged3.post)
aci.coefs.post[88,] <- c(id = "untagged_p3_sub22_q4", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.unflagged3.post)))

## Tag: 6558
tri.6558.post <- fitaci(data = subset(licor.data.post, id == "6558"),
                              varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                               Ci = "Ci", PPFD = "Qin"))
plot(tri.6558.post)
aci.coefs.post[89,] <- c(id = "6558", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.6558.post)))

## Tag: 1686
tri.1686.post <- fitaci(data = subset(licor.data.post, id == "1686"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.1686.post)
aci.coefs.post[90,] <- c(id = "1686", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.1686.post)))

## Tag: 1926
tri.1926.post <- fitaci(data = subset(licor.data.post, id == "1926"),
                        varnames = list (ALEAF = "A", Tleaf = "Tleaf", 
                                         Ci = "Ci", PPFD = "Qin"))
plot(tri.1926.post)
aci.coefs.post[91,] <- c(id = "1926", spp = "Tri", plot = "3", 
                         timepoint = "post_closure", 
                         t(coef(tri.1926.post)))

###############################################################################
## Snapshot - post canopy closure
###############################################################################
snapshot.post <- licor.data.post %>%
  mutate(machine = ifelse(is.na(machine), "yadi", machine)) %>%
  group_by(id, machine) %>%
  filter(row_number() == 1) %>%
  dplyr::select(id, machine, anet = A, ci = Ci, ca = Ca, 
                co2_ref = CO2_r, gsw, Tleaf) %>%
  mutate(ci.ca = ci / ca,
         iwue = anet / gsw) %>%
  filter(id != "5069" & id != "652" & id != "4844") %>%
  mutate(id = ifelse(id == "5069_b", "5069", id),
         id = ifelse(id == "652_b", "652", id),
         id = ifelse(id == "4844-redo", "4844", id),
         id = ifelse(id == "polka", "flag_8", id),
         id = ifelse(id == "unflagged_trillium_p3_sub22_q4", "untagged_p3_sub22_q4", id),
         id = ifelse(id == "untagged_3_sub17_q1_tri", "untagged_p3_sub17_q1", id),
         id = ifelse(id == "untagged_3_sub12_q4_tri", "untagged_p3_sub12_q4", id),
         id = ifelse(id == "548", "5488", id),
         id = ifelse(id == "1704", "784", id),
         id = ifelse(id == "235", "275", id),
         id = ifelse(id == "342", "3427", id),
         id = ifelse(id == "flag2" & machine == "stan", "flag2_mai", id),
         id = ifelse(id == "flag2" & machine == "yadi", "flag2_tri", id),
         id = ifelse(id == "flag3" & machine == "yadi", "flag3_tri", id)) %>%
  dplyr::select(id:iwue) %>% ungroup()

###############################################################################
## Merge snapshot and curve fits into central data frame
###############################################################################
post.canopy.merged <- snapshot.post %>%
  full_join(aci.coefs.post, by = "id") %>%
  mutate(plot = as.numeric(plot),
         jmax.vcmax = as.numeric(Vcmax) / as.numeric(Jmax),
         rd.vcmax = as.numeric(Rd) / as.numeric(Vcmax),
         stomatal_limitation(A_net = as.numeric(anet),
                             Vcmax = as.numeric(Vcmax),
                             Rd.meas = TRUE,
                             Rd = as.numeric(Rd),
                             leaf.temp = Tleaf,
                             temp = "C")[5]) %>%
  filter(!is.na(plot)) %>%
  full_join(spp.meta.closed) %>%
  mutate(subplot = as.numeric(subplot),
         quad = as.numeric(quad)) %>%
  full_join(gm.density.plot, by = c("plot", "subplot")) %>%
  full_join(gm.density.quad, by = c("plot", "subplot", "quad")) %>%
  mutate(rosette_subplot = ifelse(is.na(rosette_subplot), 0, rosette_subplot),
         adult_subplot = ifelse(is.na(adult_subplot), 0, adult_subplot),
         total_subplot = ifelse(is.na(total_subplot), 0, total_subplot),
         rosette_quad = ifelse(is.na(rosette_quad), 0, rosette_quad),
         adult_quad = ifelse(is.na(adult_quad), 0, rosette_quad),
         total_quad = ifelse(is.na(total_quad), 0, rosette_quad)) %>%
  dplyr::select(id, spp, plot, subplot, quad, rosette_subplot:total_quad,
                gm.trt, canopy = timepoint,
                anet, gsw, ci.ca, iwue, stom.lim = l, vcmax25 = Vcmax, 
                jmax25 = Jmax, rd25 = Rd, jmax.vcmax, 
                rd.vcmax = rd.vcmax) %>%
  filter(!is.na(ci.ca))

write.csv(post.canopy.merged, "../data/TT23_post_canopy_closure_phys_data.csv",
          row.names = FALSE)

###############################################################################
## Merge pre-canopy closure and post-canopy closure
###############################################################################
pre.canopy.merged %>%
  full_join(post.canopy.merged) %>%
  write.csv("../data/TT23_phys_data_long.csv", row.names = FALSE)


pre_canopy_long <- pre.canopy.merged %>%
  dplyr::select(id:gm.trt, anet_pre = anet, gsw_pre = gsw, ci.ca_pre = ci.ca,
                iwue_pre = iwue, stom.lim_pre = stom.lim, vcmax25_pre = vcmax25,
                jmax25_pre = jmax25, jmax.vcmax_pre = jmax.vcmax, 
                rd25_pre = rd25, rd.vcmax_pre = rd.vcmax)

post_canopy_long <- post.canopy.merged %>%
  dplyr::select(id:gm.trt, anet_post = anet, gsw_post = gsw, ci.ca_post = ci.ca,
                iwue_post = iwue, stom.lim_post = stom.lim, vcmax25_post = vcmax25,
                jmax25_post = jmax25, jmax.vcmax_post = jmax.vcmax, 
                rd25_post = rd25, rd.vcmax_post = rd.vcmax)

canopy_merged_short <- pre_canopy_long %>%
  full_join(post_canopy_long) %>%
  dplyr::select(id:gm.trt, rosette_subplot:total_quad, anet_pre:rd.vcmax_pre, 
                anet_post:rd.vcmax_post) %>%
  mutate(rosette_subplot = ifelse(is.na(rosette_subplot), 
                                  0, rosette_subplot),
         adult_subplot = ifelse(is.na(adult_subplot), 
                                  0, adult_subplot),
         total_subplot = ifelse(is.na(total_subplot), 
                                  0, total_subplot),
         rosette_quad = ifelse(gm.trt == "invaded" & is.na(rosette_quad), 
                                  0, rosette_quad),
         adult_quad = ifelse(gm.trt == "invaded" & is.na(adult_quad), 
                                  0, adult_quad),
         total_quad = ifelse(gm.trt == "invaded" & is.na(total_quad), 
                                  0, total_quad))


write.csv(canopy_merged_short, "../data/TT23_phys_data_short.csv",
          row.names = FALSE)


