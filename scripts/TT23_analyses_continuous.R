## Load libraries
library(tidyverse)
library(ggpubr)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)

## Read compiled data file
df <- read.csv("../data/TT23_compiled_datasheet.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c("invaded", "weeded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")),
         vcmax.gs = vcmax25 / gsw)
head(df)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, total_subplot, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:n_plantAvail_day),
               .funs = mean) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("invaded", "weeded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")))

##############################################################################
## Anet - Tri
##############################################################################
df$anet[c(35, 80, 86)] <- NA

anet.tri <- lmer(
  log(anet) ~ total_subplot * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri" & gm.trt == "invaded"))

# Check model assumptions
plot(anet.tri)
qqnorm(residuals(anet.tri))
qqline(residuals(anet.tri))
densityPlot(residuals(anet.tri))
shapiro.test(residuals(anet.tri))
outlierTest(anet.tri)

# Model output
summary(anet.tri)
Anova(anet.tri)
r.squaredGLMM(anet.tri)

# Pairwise comparisons
test(emtrends(anet.tri, ~canopy, "total_subplot"))
emmeans(anet.tri, pairwise~canopy)

##############################################################################
## Anet - Mai
##############################################################################
df$anet[49] <- NA

anet.mai <- lmer(
  log(anet) ~ total_subplot * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai" & gm.trt == "invaded"))

# Check model assumptions
plot(anet.mai)
qqnorm(residuals(anet.mai))
qqline(residuals(anet.mai))
densityPlot(residuals(anet.mai))
shapiro.test(residuals(anet.mai))
outlierTest(anet.mai)

# Model output
summary(anet.mai)
Anova(anet.mai)
r.squaredGLMM(anet.mai)

# Pairwise comparisons
test(emtrends(anet.mai, ~canopy, "total_subplot"))
emmeans(anet.tri, pairwise~canopy)

##############################################################################
## gs - Tri
##############################################################################
df$gsw[c(86)] <- NA

gsw.tri <- lmer(
 gsw ~ total_subplot * canopy + n_plantAvail_day + phosphate_ppm_day + 
   (1 | plot), data = subset(df, spp == "Tri" & gm.trt == "invaded"))

# Check model assumptions
plot(gsw.tri)
qqnorm(residuals(gsw.tri))
qqline(residuals(gsw.tri))
densityPlot(residuals(gsw.tri))
shapiro.test(residuals(gsw.tri))
outlierTest(gsw.tri)

# Model output
summary(gsw.tri)
Anova(gsw.tri)
r.squaredGLMM(gsw.tri)

# Pairwise comparisons
emmeans(gsw.tri, pairwise~canopy)
emmeans(gsw.tri, pairwise~gm.trt)
test(emtrends(gsw.tri, ~1, "phosphate_ppm_day"))

##############################################################################
## gs - Mai
##############################################################################
df$gsw[c(120)] <- NA

gsw.mai <- lmer(
  gsw ~ total_subplot * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(gsw.mai)
qqnorm(residuals(gsw.mai))
qqline(residuals(gsw.mai))
densityPlot(residuals(gsw.mai))
shapiro.test(residuals(gsw.mai))
outlierTest(gsw.mai)

# Model output
summary(gsw.mai)
Anova(gsw.mai)
r.squaredGLMM(gsw.mai)

# Pairwise comparisons
test(emtrends(anet.mai, ~canopy, "total_subplot"))
emmeans(anet.tri, pairwise~canopy)

##############################################################################
## stomatal limitation - Tri
##############################################################################
df$stom.lim[c(228, 229)] <- NA

stomlim.tri <- lmer(
  log(stom.lim) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri" & stom.lim > 0))

# Check model assumptions
plot(stomlim.tri)
qqnorm(residuals(stomlim.tri))
qqline(residuals(stomlim.tri))
densityPlot(residuals(stomlim.tri))
shapiro.test(residuals(stomlim.tri))
outlierTest(stomlim.tri)

# Model output
summary(stomlim.tri)
Anova(stomlim.tri)
r.squaredGLMM(stomlim.tri)

# Pairwise comparisons
emmeans(stomlim.tri, pairwise~canopy, type = "response")
emmeans(stomlim.tri, pairwise~gm.trt, type = "response")
test(emtrends(stomlim.tri, ~1, "phosphate_ppm_day"))

##############################################################################
## stomatal limitation - Mai
##############################################################################
stom.lim.mai <- lmer(
  log(stom.lim) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai" & stom.lim > 0))

# Check model assumptions
plot(stom.lim.mai)
qqnorm(residuals(stom.lim.mai))
qqline(residuals(stom.lim.mai))
densityPlot(residuals(stom.lim.mai))
shapiro.test(residuals(stom.lim.mai))
outlierTest(stom.lim.mai)

# Model output
summary(stom.lim.mai)
Anova(stom.lim.mai)
r.squaredGLMM(stom.lim.mai)

# Pairwise comparisons
emmeans(stom.lim.mai, pairwise~gm.trt, type = "response")
cld(emmeans(stom.lim.mai, pairwise~gm.trt*canopy, type = "response"))

##############################################################################
## Vcmax - Tri
##############################################################################
df$vcmax25[183] <- NA

vcmax.tri <- lmer(
  log(vcmax25) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(vcmax.tri)
qqnorm(residuals(vcmax.tri))
qqline(residuals(vcmax.tri))
densityPlot(residuals(vcmax.tri))
shapiro.test(residuals(vcmax.tri))
outlierTest(vcmax.tri)

# Model output
summary(vcmax.tri)
Anova(vcmax.tri)
r.squaredGLMM(vcmax.tri)

# Pairwise comparisons
emmeans(vcmax.tri, pairwise~gm.trt, type = "response")
cld(emmeans(vcmax.tri, pairwise~gm.trt*canopy))

##############################################################################
## Vcmax - Mai
##############################################################################
df$vcmax25[231] <- NA

vcmax.mai <- lmer(
  log(vcmax25) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(vcmax.mai)
qqnorm(residuals(vcmax.mai))
qqline(residuals(vcmax.mai))
densityPlot(residuals(vcmax.mai))
shapiro.test(residuals(vcmax.mai))
outlierTest(vcmax.mai)

# Model output
summary(vcmax.mai)
Anova(vcmax.mai)
r.squaredGLMM(vcmax.mai)

##############################################################################
## Jmax - Tri
##############################################################################
df$jmax25[183] <- NA

jmax.tri <- lmer(
  log(jmax25) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(jmax.tri)
qqnorm(residuals(jmax.tri))
qqline(residuals(jmax.tri))
densityPlot(residuals(jmax.tri))
shapiro.test(residuals(jmax.tri))
outlierTest(jmax.tri)

# Model output
summary(jmax.tri)
Anova(jmax.tri)
r.squaredGLMM(jmax.tri)

# Pairwise comparisons
cld(emmeans(jmax.tri, pairwise~gm.trt*canopy))

##############################################################################
## Jmax - Mai
##############################################################################
jmax.mai <- lmer(
  jmax25 ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(jmax.mai)
qqnorm(residuals(jmax.mai))
qqline(residuals(jmax.mai))
densityPlot(residuals(jmax.mai))
shapiro.test(residuals(jmax.mai))
outlierTest(jmax.mai)

# Model output
summary(jmax.mai)
Anova(jmax.mai)
r.squaredGLMM(jmax.mai)

# Pairwise comparisons
emmeans(jmax.mai, pairwise~canopy)
test(emtrends(jmax.mai, ~1, "n_plantAvail_day"))

##############################################################################
## Jmax : Vcmax - Tri
##############################################################################
jmax.vcmax.tri <- lmer(
  log(jmax.vcmax) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(jmax.vcmax.tri)
qqnorm(residuals(jmax.vcmax.tri))
qqline(residuals(jmax.vcmax.tri))
densityPlot(residuals(jmax.vcmax.tri))
shapiro.test(residuals(jmax.vcmax.tri))
outlierTest(jmax.vcmax.tri)

# Model output
summary(jmax.vcmax.tri)
Anova(jmax.vcmax.tri)
r.squaredGLMM(jmax.vcmax.tri)

##############################################################################
## Jmax : Vcmax - Mai
##############################################################################
df$jmax.vcmax[c(94, 225, 231)] <- NA

jmax.vcmax.mai <- lmer(
  log(jmax.vcmax) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(jmax.vcmax.mai)
qqnorm(residuals(jmax.vcmax.mai))
qqline(residuals(jmax.vcmax.mai))
densityPlot(residuals(jmax.vcmax.mai))
shapiro.test(residuals(jmax.vcmax.mai))
outlierTest(jmax.vcmax.mai)

# Model output
summary(jmax.vcmax.mai)
Anova(jmax.vcmax.mai)
r.squaredGLMM(jmax.vcmax.mai)

# Pairwise comparisons
emmeans(jmax.vcmax.mai, pairwise~canopy, type = "response")
emmeans(jmax.vcmax.mai, pairwise~gm.trt, type = "response")

##############################################################################
## iWUE - Tri
##############################################################################
df$iwue[c(80, 86, 158, 159)] <- NA

iwue.tri <- lmer(
  log(iwue) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(iwue.tri)
qqnorm(residuals(iwue.tri))
qqline(residuals(iwue.tri))
densityPlot(residuals(iwue.tri))
shapiro.test(residuals(iwue.tri))
outlierTest(iwue.tri)

# Model output
summary(iwue.tri)
Anova(iwue.tri)
r.squaredGLMM(iwue.tri)

# Pairwise comparisons
emmeans(iwue.tri, pairwise~canopy)
emmeans(iwue.tri, pairwise~gm.trt, type = "response")
test(emtrends(iwue.tri, ~1, "phosphate_ppm_day"))

##############################################################################
## iWUE - Mai
##############################################################################
df$iwue[c(104)] <- NA

iwue.mai <- lmer(
  log(iwue) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(iwue.mai)
qqnorm(residuals(iwue.mai))
qqline(residuals(iwue.mai))
densityPlot(residuals(iwue.mai))
shapiro.test(residuals(iwue.mai))
outlierTest(iwue.mai)

# Model output
summary(iwue.mai)
Anova(iwue.mai)
r.squaredGLMM(iwue.mai)

# Pairwise comparisons
emmeans(iwue.mai, pairwise~canopy)

##############################################################################
## Vcmax25:gs - Tri
##############################################################################
vcmax.gs.tri <- lmer(
  log(vcmax.gs) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(vcmax.gs.tri)
qqnorm(residuals(vcmax.gs.tri))
qqline(residuals(vcmax.gs.tri))
densityPlot(residuals(vcmax.gs.tri))
shapiro.test(residuals(vcmax.gs.tri))
outlierTest(vcmax.gs.tri)

# Model output
summary(vcmax.gs.tri)
Anova(vcmax.gs.tri)
r.squaredGLMM(vcmax.gs.tri)

# Pairwise comparisons
emmeans(vcmax.gs.tri, pairwise~canopy)
test(emtrends(vcmax.gs.tri, ~1, "phosphate_ppm_day"))

##############################################################################
## Vcmax25:gs - Mai
##############################################################################
vcmax.gs.mai <- lmer(
  log(vcmax.gs) ~ gm.trt * canopy + n_plantAvail_day + phosphate_ppm_day + 
    (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(vcmax.gs.mai)
qqnorm(residuals(vcmax.gs.mai))
qqline(residuals(vcmax.gs.mai))
densityPlot(residuals(vcmax.gs.mai))
shapiro.test(residuals(vcmax.gs.mai))
outlierTest(vcmax.gs.mai)

# Model output
summary(vcmax.gs.mai)
Anova(vcmax.gs.mai)
r.squaredGLMM(vcmax.gs.mai)

# Pairwise comparisons
cld(emmeans(vcmax.gs.mai, pairwise~gm.trt*canopy, type = "response"))

##############################################################################
## Write Table 1: Soil nutrients
##############################################################################
soil.nitrogen <- data.frame(Anova(plant_availableN)) %>%
  mutate(treatment = row.names(.),
         chisq_soilN = Chisq,
         p_soilN = Pr..Chisq.,
         across(chisq_soilN:p_soilN, round, digits = 3),
         chisq_soilN = ifelse(chisq_soilN <0.001 & chisq_soilN >= 0, 
                                "<0.001", chisq_soilN),
         p_soilN = ifelse(p_soilN <0.001 & p_soilN >= 0, 
                               "<0.001", p_soilN)) %>%
  dplyr::select(treatment, Df, chisq_soilN, p_soilN)

soil.nitrate <- data.frame(Anova(nitrate)) %>%
  mutate(treatment = row.names(.),
         chisq_nitrate = Chisq,
         p_nitrate = Pr..Chisq.,
         across(chisq_nitrate:p_nitrate, round, digits = 3),
         chisq_nitrate = ifelse(chisq_nitrate < 0.001 & chisq_nitrate >= 0, 
                              "<0.001", chisq_nitrate),
         p_nitrate = ifelse(p_nitrate <0.001 & p_nitrate >= 0, 
                          "<0.001", p_nitrate)) %>%
  dplyr::select(treatment, chisq_nitrate, p_nitrate)

soil.ammonium <- data.frame(Anova(ammonium)) %>%
  mutate(treatment = row.names(.),
         chisq_ammonium = Chisq,
         p_ammonium = Pr..Chisq.,
         across(chisq_ammonium:p_ammonium, round, digits = 3),
         chisq_ammonium = ifelse(chisq_ammonium < 0.001 & chisq_ammonium >= 0, 
                                "<0.001", chisq_ammonium),
         p_ammonium = ifelse(p_ammonium <0.001 & p_ammonium >= 0, 
                            "<0.001", p_ammonium)) %>%
  dplyr::select(treatment, chisq_ammonium, p_ammonium)

soil.phosphate <- data.frame(Anova(phosphate)) %>%
  mutate(treatment = row.names(.),
         chisq_phosphate = Chisq,
         p_phosphate = Pr..Chisq.,
         across(chisq_phosphate:p_phosphate, round, digits = 3),
         chisq_phosphate = ifelse(chisq_phosphate < 0.001 & chisq_phosphate >= 0, 
                                 "<0.001", chisq_phosphate),
         p_phosphate = ifelse(p_phosphate <0.001 & p_phosphate >= 0, 
                             "<0.001", p_phosphate)) %>%
  dplyr::select(treatment, chisq_phosphate, p_phosphate)

table1 <- soil.nitrogen %>% full_join(soil.nitrate) %>% 
  full_join(soil.ammonium) %>% full_join(soil.phosphate)

write.csv(table1, "../drafts/tables/TT23_table1_soil_nutrients.csv",
          row.names = FALSE)

##############################################################################
## Write Table 2: Gas exchange
##############################################################################
anet.tri <- data.frame(Anova(anet.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_anet.tri = Chisq,
         p_anet.tri = Pr..Chisq.,
         across(chisq_anet.tri:p_anet.tri, round, digits = 3),
         chisq_anet.tri = ifelse(chisq_anet.tri < 0.001 & chisq_anet.tri >= 0, 
                              "<0.001", chisq_anet.tri),
         p_anet.tri = ifelse(p_anet.tri <0.001 & p_anet.tri >= 0, 
                          "<0.001", p_anet.tri)) %>%
  dplyr::select(treatment, Df, chisq_anet.tri, p_anet.tri)

anet.mai <- data.frame(Anova(anet.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_anet.mai = Chisq,
         p_anet.mai = Pr..Chisq.,
         across(chisq_anet.mai:p_anet.mai, round, digits = 3),
         chisq_anet.mai = ifelse(chisq_anet.mai < 0.001 & chisq_anet.mai >= 0, 
                                 "<0.001", chisq_anet.mai),
         p_anet.mai = ifelse(p_anet.mai <0.001 & p_anet.mai >= 0, 
                             "<0.001", p_anet.mai)) %>%
  dplyr::select(treatment, chisq_anet.mai, p_anet.mai)

gsw.tri <- data.frame(Anova(gsw.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_gsw.tri = Chisq,
         p_gsw.tri = Pr..Chisq.,
         across(chisq_gsw.tri:p_gsw.tri, round, digits = 3),
         chisq_gsw.tri = ifelse(chisq_gsw.tri < 0.001 & chisq_gsw.tri >= 0, 
                                 "<0.001", chisq_gsw.tri),
         p_gsw.tri = ifelse(p_gsw.tri <0.001 & p_gsw.tri >= 0, 
                             "<0.001", p_gsw.tri)) %>%
  dplyr::select(treatment, chisq_gsw.tri, p_gsw.tri)

gsw.mai <- data.frame(Anova(gsw.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_gsw.mai = Chisq,
         p_gsw.mai = Pr..Chisq.,
         across(chisq_gsw.mai:p_gsw.mai, round, digits = 3),
         chisq_gsw.mai = ifelse(chisq_gsw.mai < 0.001 & chisq_gsw.mai >= 0, 
                                 "<0.001", chisq_gsw.mai),
         p_gsw.mai = ifelse(p_gsw.mai <0.001 & p_gsw.mai >= 0, 
                             "<0.001", p_gsw.mai)) %>%
  dplyr::select(treatment, chisq_gsw.mai, p_gsw.mai)


stom.lim.tri <- data.frame(Anova(stomlim.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_stom.lim.tri = Chisq,
         p_stom.lim.tri = Pr..Chisq.,
         across(chisq_stom.lim.tri:p_stom.lim.tri, round, digits = 3),
         chisq_stom.lim.tri = ifelse(chisq_stom.lim.tri < 0.001 & 
                                       chisq_stom.lim.tri >= 0, 
                                "<0.001", chisq_stom.lim.tri),
         p_stom.lim.tri = ifelse(p_stom.lim.tri <0.001 & p_stom.lim.tri >= 0, 
                            "<0.001", p_stom.lim.tri)) %>%
  dplyr::select(treatment, chisq_stom.lim.tri, p_stom.lim.tri)

stom.lim.mai <- data.frame(Anova(stom.lim.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_stom.lim.mai = Chisq,
         p_stom.lim.mai = Pr..Chisq.,
         across(chisq_stom.lim.mai:p_stom.lim.mai, round, digits = 3),
         chisq_stom.lim.mai = ifelse(chisq_stom.lim.mai < 0.001 & 
                                       chisq_stom.lim.mai >= 0, 
                                "<0.001", chisq_stom.lim.mai),
         p_stom.lim.mai = ifelse(p_stom.lim.mai <0.001 & p_stom.lim.mai >= 0, 
                            "<0.001", p_stom.lim.mai)) %>%
  dplyr::select(treatment, chisq_stom.lim.mai, p_stom.lim.mai)

table2 <- anet.tri %>% full_join(gsw.tri) %>% full_join(stom.lim.tri) %>% 
  full_join(anet.mai) %>% full_join(gsw.mai) %>% full_join(stom.lim.mai)

write.csv(table2, "../drafts/tables/TT23_table2_gas_exchange.csv",
          row.names = FALSE)









