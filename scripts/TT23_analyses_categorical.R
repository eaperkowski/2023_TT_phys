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

# Helper fxn to change "NaN" to "NA"
NaN_to_NA <- function(x) ifelse(is.nan(x), NA, x)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:inorg_n_ppm),
               .funs = mean, na.rm = TRUE) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("invaded", "weeded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure"))) %>%
  mutate(across(nitrate_ppm:inorg_n_ppm, .fns = NaN_to_NA),
         np.ratio = inorg_n_ppm/phosphate_ppm)

##############################################################################
## Nitrate
##############################################################################
nitrate <- lmer(nitrate_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(nitrate)
qqnorm(residuals(nitrate))
qqline(residuals(nitrate))
densityPlot(residuals(nitrate))
shapiro.test(residuals(nitrate))
outlierTest(nitrate)

# Model output
summary(nitrate)
Anova(nitrate)
r.squaredGLMM(nitrate)

# Pairwise comparisons
emmeans(nitrate, pairwise~canopy)

# % change nitrate with canopy status
(5.32 - 18.57) / 18.57 * 100

##############################################################################
## Ammonium
##############################################################################
df.soil$ammonium_ppm[c(40, 56)] <- NA

ammonium <- lmer(ammonium_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(ammonium)
qqnorm(residuals(ammonium))
qqline(residuals(ammonium))
densityPlot(residuals(ammonium))
shapiro.test(residuals(ammonium))
outlierTest(ammonium)

# Model output
summary(ammonium)
Anova(ammonium)
r.squaredGLMM(ammonium)

# Pairwise comparisons
cld(emmeans(ammonium, pairwise~canopy*gm.trt), Letters = LETTERS)
emmeans(ammonium, pairwise~canopy)

##############################################################################
## Phosphate
##############################################################################
phosphate <- lmer(
  phosphate_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(phosphate)
qqnorm(residuals(phosphate))
qqline(residuals(phosphate))
densityPlot(residuals(phosphate))
shapiro.test(residuals(phosphate))
outlierTest(phosphate)

# Model output
summary(phosphate)
Anova(phosphate)
r.squaredGLMM(phosphate)

# Pairwise comparisons
emmeans(phosphate, pairwise~gm.trt)
emmeans(phosphate, pairwise~canopy)

# % change canopy
(0.81 - 1.09) / 1.09 * 100

##############################################################################
## N availability (nitrate + ammonium)
##############################################################################
plant_availableN <- lmer(
  inorg_n_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(plant_availableN)
qqnorm(residuals(plant_availableN))
qqline(residuals(plant_availableN))
densityPlot(residuals(plant_availableN))
shapiro.test(residuals(plant_availableN))
outlierTest(plant_availableN)

# Model output
summary(plant_availableN)
Anova(plant_availableN)
r.squaredGLMM(plant_availableN)

# Pairwise comparisons
cld(emmeans(plant_availableN, pairwise~gm.trt*canopy))
emmeans(plant_availableN, pairwise~gm.trt)
emmeans(plant_availableN, pairwise~canopy)

# Percent change due to canopy
(5.62 - 18.57) / 18.57 * 100

##############################################################################
## Soil N:P
##############################################################################
df.soil$np.ratio[c(55)] <- NA

n_to_p_ratio <- lmer(
  np.ratio ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(n_to_p_ratio)
qqnorm(residuals(n_to_p_ratio))
qqline(residuals(n_to_p_ratio))
densityPlot(residuals(n_to_p_ratio))
shapiro.test(residuals(n_to_p_ratio))
outlierTest(n_to_p_ratio)

# Model output
summary(n_to_p_ratio)
Anova(n_to_p_ratio)
r.squaredGLMM(n_to_p_ratio)

# Pairwise comparisons
emmeans(n_to_p_ratio, pairwise~gm.trt)
emmeans(n_to_p_ratio, pairwise~canopy)

# Percent change due to GM treatment
(17.1 - 10.2) / 10.2 * 100

# Percent change due to canopy
(9.2 - 18.1) / 18.1 * 100

##############################################################################
## Anet - Tri
##############################################################################
df$anet[c(35, 80, 86, 116)] <- NA

anet.tri <- lmer(
  anet ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
emmeans(anet.tri, pairwise~canopy, type = "response")

# Percent change
(4.58 - 12.4) / 12.4 * 100

##############################################################################
## Anet - Mai
##############################################################################
df$anet[c(71, 120)] <- NA

anet.mai <- lmer(
  anet ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(anet.mai, pairwise~canopy)
emmeans(anet.mai, pairwise~gm.trt)

# % change canopy
(4.12 - 9.28) / 9.28 * 100

# % change gm.trt
(6.08 - 7.31) / 7.31 * 100

##############################################################################
## gs - Tri
##############################################################################
df$gsw[86] <- NA

gsw.tri <- lmer(
  gsw ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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

# Canopy % change
(0.113 - 0.139) / 0.139 * 100

##############################################################################
## gs - Mai
##############################################################################
df$gsw[120] <- NA

gsw.mai <- lmer(
  gsw ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(gsw.mai, pairwise~canopy)
emmeans(gsw.mai, pairwise~gm.trt)

# Canopy % change
(0.0594 - 0.1519) / 0.1519 * 100

# % change gm.trt
(0.0871 - 0.1242) / 0.1242 * 100


##############################################################################
## stomatal limitation - Tri
##############################################################################
df$stom.lim[c(228, 229)] <- NA

stomlim.tri <- lmer(log(stom.lim) ~ gm.trt * canopy  + (1 | plot),
                    data = subset(df, spp == "Tri" & stom.lim > 0))

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

# % change canopy
(0.177 - 0.492) / 0.492 * 100

# % change gm.trt
(0.287 - 0.303) / 0.287

##############################################################################
## stomatal limitation - Mai
##############################################################################
stom.lim.mai <- lmer(log(stom.lim) ~ gm.trt * canopy + (1 | plot), 
                     data = subset(df, spp == "Mai" & stom.lim > 0))

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

# % change gm.trt
(0.356 - 0.273) / 0.273 * 100

##############################################################################
## SPAD - Tri
##############################################################################
spad.tri <- lmer(
  SPAD ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(spad.tri)
qqnorm(residuals(spad.tri))
qqline(residuals(spad.tri))
densityPlot(residuals(spad.tri))
shapiro.test(residuals(spad.tri))
outlierTest(spad.tri)

# Model output
summary(spad.tri)
Anova(spad.tri)
r.squaredGLMM(spad.tri)

# Pairwise comparisons
emmeans(spad.tri, pairwise~canopy)

# % change canopy
(44.7 - 35.5) / 35.5 * 100

##############################################################################
## SPAD - Mai
##############################################################################
df$SPAD[c(122)] <- NA

spad.mai <- lmer(
  SPAD ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(spad.mai)
qqnorm(residuals(spad.mai))
qqline(residuals(spad.mai))
densityPlot(residuals(spad.mai))
shapiro.test(residuals(spad.mai))
outlierTest(spad.mai)

# Model output
summary(spad.mai)
Anova(spad.mai)
r.squaredGLMM(spad.mai)

# Pairwise comparisons
emmeans(spad.mai, pairwise~canopy)

# % change canopy
(39.9 - 26.5) / 26.5 * 100

##############################################################################
## Vcmax - Tri
##############################################################################
df$vcmax25[c(183)] <- NA

vcmax.tri <- lmer(
  log(vcmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
emmeans(vcmax.tri, pairwise~canopy, type = "response")
cld(emmeans(vcmax.tri, pairwise~gm.trt*canopy, type = "response"))

# % change canopy
(23.3 - 101.8) / 101.8 * 100

##############################################################################
## Vcmax - Mai
##############################################################################
df$vcmax25[231] <- NA

vcmax.mai <- lmer(
  log(vcmax25) ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

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

# Pairwise comparisons
emmeans(vcmax.mai, pairwise~canopy, type = "response")
emmeans(vcmax.mai, pairwise~gm.trt, type = "response")

# % change canopy
(23.3 - 56.2) / 56.6 * 100

##############################################################################
## Jmax - Tri
##############################################################################
df$jmax25[183] <- NA

jmax.tri <- lmer(
  log(jmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
emmeans(jmax.tri, pairwise~canopy, type = "response")
cld(emmeans(jmax.tri, pairwise~gm.trt*canopy))

# % change canopy
(43.2 - 179.5) / 179.5 * 100

##############################################################################
## Jmax - Mai
##############################################################################
jmax.mai <- lmer(
  log(jmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(jmax.mai, pairwise~canopy, type = "response")
emmeans(jmax.mai, pairwise~gm.trt, type = "response")

# % change canopy
(38.8 - 99.8) / 99.8 * 100

##############################################################################
## Jmax : Vcmax - Tri
##############################################################################
df$jmax.vcmax[184] <- NA

jmax.vcmax.tri <- lmer(
  log(jmax.vcmax) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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

# Pairwise comparisons
emmeans(jmax.vcmax.tri, pairwise~canopy, type = "response")

# % change canopy
(0.543 - 0.567) / 0.567 * 100

##############################################################################
## Jmax : Vcmax - Mai
##############################################################################
df$jmax.vcmax[c(94, 224, 225, 231)] <- NA

jmax.vcmax.mai <- lmer(
  log(jmax.vcmax) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

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

# % change canopy
(0.6 - 0.559) / 0.559 * 100 

# % change gm.trt
(0.592 - 0.567) / 0.567 * 100

##############################################################################
## iWUE - Tri
##############################################################################
df$iwue[c(80, 86, 116)] <- NA

iwue.tri <- lmer(
  iwue ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
emmeans(iwue.tri, pairwise~gm.trt)

# % change canopy
(44.6 - 91.3) / 91.3 * 100

##############################################################################
## iWUE - Mai
##############################################################################
iwue.mai <- lmer(
  iwue ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(iwue.mai, pairwise~gm.trt)

# % change canopy
(78.9 - 62.6) / 62.6 * 100 

# % change gm.trt
(82.6 - 58.9) / 58.9 * 100

##############################################################################
## SPAD:gs - Tri
##############################################################################
df$spad.gs <- df$SPAD / df$gsw

spad.gs.tri <- lmer(
  log(spad.gs) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(spad.gs.tri)
qqnorm(residuals(spad.gs.tri))
qqline(residuals(spad.gs.tri))
densityPlot(residuals(spad.gs.tri))
shapiro.test(residuals(spad.gs.tri))
outlierTest(spad.gs.tri)

# Model output
summary(spad.gs.tri)
Anova(spad.gs.tri)
r.squaredGLMM(spad.gs.tri)

# Pairwise comparisons
emmeans(spad.gs.tri, pairwise~canopy, type = "response")

# % change canopy
(410 - 251) / 251 * 100


##############################################################################
## SPAD:gs - Mai
##############################################################################
df$spad.gs[c(49, 132)] <- NA

spad.gs.mai <- lmer(
  log(spad.gs) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(spad.gs.mai)
qqnorm(residuals(spad.gs.mai))
qqline(residuals(spad.gs.mai))
densityPlot(residuals(spad.gs.mai))
shapiro.test(residuals(spad.gs.mai))
outlierTest(spad.gs.mai)

# Model output
summary(spad.gs.mai)
Anova(spad.gs.mai)
r.squaredGLMM(spad.gs.mai)

# Pairwise comparisons
emmeans(spad.gs.mai, pairwise~canopy)
emmeans(spad.gs.mai, pairwise~gm.trt, type = "response")
cld(emmeans(spad.gs.mai, pairwise~gm.trt*canopy))

# % change canopy
(6.73 - 5.19) / 5.19 * 100

# % change gm.trt
(497 - 301) / 301 * 100


##############################################################################
## Vcmax25:gs - Tri
##############################################################################
vcmax.gs.tri <- lmer(
  log(vcmax.gs) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
emmeans(vcmax.gs.tri, pairwise~canopy, type = "response")

# % change canopy
(221 - 739) / 739 * 100

##############################################################################
## Vcmax25:gs - Mai
##############################################################################
vcmax.gs.mai <- lmer(
  log(vcmax.gs) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(vcmax.gs.mai, pairwise~canopy, type = "response")
emmeans(vcmax.gs.mai, pairwise~gm.trt, type = "response")
cld(emmeans(vcmax.gs.mai, pairwise~gm.trt*canopy))

# % change canopy
(463 -  370) / 370 * 100

# % change gm.trt
(519 - 330) / 330 * 100


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

##############################################################################
## Write Table 3: Indices of photosynthetic capacity
##############################################################################
vcmax.tri <- data.frame(Anova(vcmax.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_vcmax.tri = Chisq,
         p_vcmax.tri = Pr..Chisq.,
         across(chisq_vcmax.tri:p_vcmax.tri, round, digits = 3),
         chisq_vcmax.tri = ifelse(chisq_vcmax.tri < 0.001 & chisq_vcmax.tri >= 0, 
                                 "<0.001", chisq_vcmax.tri),
         p_vcmax.tri = ifelse(p_vcmax.tri <0.001 & p_vcmax.tri >= 0, 
                             "<0.001", p_vcmax.tri)) %>%
  dplyr::select(treatment, Df, chisq_vcmax.tri, p_vcmax.tri)

vcmax.mai <- data.frame(Anova(vcmax.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_vcmax.mai = Chisq,
         p_vcmax.mai = Pr..Chisq.,
         across(chisq_vcmax.mai:p_vcmax.mai, round, digits = 3),
         chisq_vcmax.mai = ifelse(chisq_vcmax.mai < 0.001 & chisq_vcmax.mai >= 0, 
                                 "<0.001", chisq_vcmax.mai),
         p_vcmax.mai = ifelse(p_vcmax.mai <0.001 & p_vcmax.mai >= 0, 
                             "<0.001", p_vcmax.mai)) %>%
  dplyr::select(treatment, chisq_vcmax.mai, p_vcmax.mai)

jmax.tri <- data.frame(Anova(jmax.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_jmax.tri = Chisq,
         p_jmax.tri = Pr..Chisq.,
         across(chisq_jmax.tri:p_jmax.tri, round, digits = 3),
         chisq_jmax.tri = ifelse(chisq_jmax.tri < 0.001 & chisq_jmax.tri >= 0, 
                                  "<0.001", chisq_jmax.tri),
         p_jmax.tri = ifelse(p_jmax.tri <0.001 & p_jmax.tri >= 0, 
                              "<0.001", p_jmax.tri)) %>%
  dplyr::select(treatment, Df, chisq_jmax.tri, p_jmax.tri)

jmax.mai <- data.frame(Anova(jmax.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_jmax.mai = Chisq,
         p_jmax.mai = Pr..Chisq.,
         across(chisq_jmax.mai:p_jmax.mai, round, digits = 3),
         chisq_jmax.mai = ifelse(chisq_jmax.mai < 0.001 & chisq_jmax.mai >= 0, 
                                  "<0.001", chisq_jmax.mai),
         p_jmax.mai = ifelse(p_jmax.mai <0.001 & p_jmax.mai >= 0, 
                              "<0.001", p_jmax.mai)) %>%
  dplyr::select(treatment, chisq_jmax.mai, p_jmax.mai)

jvmax.tri <- data.frame(Anova(jmax.vcmax.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_jvmax.tri = Chisq,
         p_jvmax.tri = Pr..Chisq.,
         across(chisq_jvmax.tri:p_jvmax.tri, round, digits = 3),
         chisq_jvmax.tri = ifelse(chisq_jvmax.tri < 0.001 & chisq_jvmax.tri >= 0, 
                                 "<0.001", chisq_jvmax.tri),
         p_jvmax.tri = ifelse(p_jvmax.tri <0.001 & p_jvmax.tri >= 0, 
                             "<0.001", p_jvmax.tri)) %>%
  dplyr::select(treatment, Df, chisq_jvmax.tri, p_jvmax.tri)

jvmax.mai <- data.frame(Anova(jmax.vcmax.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_jvmax.mai = Chisq,
         p_jvmax.mai = Pr..Chisq.,
         across(chisq_jvmax.mai:p_jvmax.mai, round, digits = 3),
         chisq_jvmax.mai = ifelse(chisq_jvmax.mai < 0.001 & chisq_jvmax.mai >= 0, 
                                 "<0.001", chisq_jvmax.mai),
         p_jvmax.mai = ifelse(p_jvmax.mai <0.001 & p_jvmax.mai >= 0, 
                             "<0.001", p_jvmax.mai)) %>%
  dplyr::select(treatment, chisq_jvmax.mai, p_jvmax.mai)

table3 <- vcmax.tri %>% full_join(jmax.tri) %>%  full_join(jvmax.tri) %>% 
  full_join(vcmax.mai) %>%  full_join(jmax.mai) %>% full_join(jvmax.mai)
write.csv(table3, "../drafts/tables/TT23_table3_photoCapacity.csv", 
          row.names = FALSE)

##############################################################################
## Write Table 4: SPAD and phi psII
##############################################################################
spad.tri <- data.frame(Anova(spad.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_spad.tri = Chisq,
         p_spad.tri = Pr..Chisq.,
         across(chisq_spad.tri:p_spad.tri, round, digits = 3),
         chisq_spad.tri = ifelse(chisq_spad.tri < 0.001 & chisq_spad.tri >= 0, 
                                 "<0.001", chisq_spad.tri),
         p_spad.tri = ifelse(p_spad.tri <0.001 & p_spad.tri >= 0, 
                             "<0.001", p_spad.tri)) %>%
  dplyr::select(treatment, Df, chisq_spad.tri, p_spad.tri)

spad.mai <- data.frame(Anova(spad.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_spad.mai = Chisq,
         p_spad.mai = Pr..Chisq.,
         across(chisq_spad.mai:p_spad.mai, round, digits = 3),
         chisq_spad.mai = ifelse(chisq_spad.mai < 0.001 & chisq_spad.mai >= 0, 
                                 "<0.001", chisq_spad.mai),
         p_spad.mai = ifelse(p_spad.mai <0.001 & p_spad.mai >= 0, 
                             "<0.001", p_spad.mai)) %>%
  dplyr::select(treatment, chisq_spad.mai, p_spad.mai)

phi2.tri <- data.frame(Anova(phips2.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_phi2.tri = Chisq,
         p_phi2.tri = Pr..Chisq.,
         across(chisq_phi2.tri:p_phi2.tri, round, digits = 3),
         chisq_phi2.tri = ifelse(chisq_phi2.tri < 0.001 & chisq_phi2.tri >= 0, 
                                 "<0.001", chisq_phi2.tri),
         p_phi2.tri = ifelse(p_phi2.tri < 0.001 & p_phi2.tri >= 0, 
                             "<0.001", p_phi2.tri)) %>%
  dplyr::select(treatment, Df, chisq_phi2.tri, p_phi2.tri)

phi2.mai <- data.frame(Anova(phips2.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_phi2.mai = Chisq,
         p_phi2.mai = Pr..Chisq.,
         across(chisq_phi2.mai:p_phi2.mai, round, digits = 3),
         chisq_phi2.mai = ifelse(chisq_phi2.mai < 0.001 & chisq_phi2.mai >= 0, 
                                 "<0.001", chisq_phi2.mai),
         p_phi2.mai = ifelse(p_phi2.mai < 0.001 & p_phi2.mai >= 0, 
                             "<0.001", p_phi2.mai)) %>%
  dplyr::select(treatment, chisq_phi2.mai, p_phi2.mai)

table4 <- spad.tri %>% full_join(phi2.tri) %>% full_join(spad.mai) %>%
  full_join(phi2.mai)
write.csv(table4, "../drafts/tables/TT23_table4_chlor_fluor.csv",
          row.names = FALSE)

##############################################################################
## Write Table 5: iWUE and Vcmax:gsw
##############################################################################
iwue.tri <- data.frame(Anova(iwue.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_iwue.tri = Chisq,
         p_iwue.tri = Pr..Chisq.,
         across(chisq_iwue.tri:p_iwue.tri, round, digits = 3),
         chisq_iwue.tri = ifelse(chisq_iwue.tri < 0.001 & chisq_iwue.tri >= 0, 
                                  "<0.001", chisq_iwue.tri),
         p_iwue.tri = ifelse(p_iwue.tri <0.001 & p_iwue.tri >= 0, 
                              "<0.001", p_iwue.tri)) %>%
  dplyr::select(treatment, Df, chisq_iwue.tri, p_iwue.tri)

iwue.mai <- data.frame(Anova(iwue.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_iwue.mai = Chisq,
         p_iwue.mai = Pr..Chisq.,
         across(chisq_iwue.mai:p_iwue.mai, round, digits = 3),
         chisq_iwue.mai = ifelse(chisq_iwue.mai < 0.001 & chisq_iwue.mai >= 0, 
                                  "<0.001", chisq_iwue.mai),
         p_iwue.mai = ifelse(p_iwue.mai <0.001 & p_iwue.mai >= 0, 
                              "<0.001", p_iwue.mai)) %>%
  dplyr::select(treatment, chisq_iwue.mai, p_iwue.mai)

vcmax.gs.tri <- data.frame(Anova(vcmax.gs.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_vcmax.gs.tri = Chisq,
         p_vcmax.gs.tri = Pr..Chisq.,
         across(chisq_vcmax.gs.tri:p_vcmax.gs.tri, round, digits = 3),
         chisq_vcmax.gs.tri = ifelse(chisq_vcmax.gs.tri < 0.001 & 
                                       chisq_vcmax.gs.tri >= 0, 
                                     "<0.001", chisq_vcmax.gs.tri),
         p_iwue.tri = ifelse(p_vcmax.gs.tri <0.001 & p_vcmax.gs.tri >= 0, 
                             "<0.001", p_vcmax.gs.tri)) %>%
  dplyr::select(treatment, Df, chisq_vcmax.gs.tri, p_vcmax.gs.tri)

vcmax.gs.mai <- data.frame(Anova(vcmax.gs.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_vcmax.gs.mai = Chisq,
         p_vcmax.gs.mai = Pr..Chisq.,
         across(chisq_vcmax.gs.mai:p_vcmax.gs.mai, round, digits = 3),
         chisq_vcmax.gs.mai = ifelse(chisq_vcmax.gs.mai < 0.001 & chisq_vcmax.gs.mai >= 0, 
                                 "<0.001", chisq_vcmax.gs.mai),
         p_vcmax.gs.mai = ifelse(p_vcmax.gs.mai <0.001 & p_vcmax.gs.mai >= 0, 
                             "<0.001", p_vcmax.gs.mai)) %>%
  dplyr::select(treatment, chisq_vcmax.gs.mai, p_vcmax.gs.mai)

table5 <- iwue.tri %>% full_join(vcmax.gs.tri) %>% full_join(iwue.mai) %>%
   full_join(vcmax.gs.mai)
write.csv(table5, "../drafts/tables/TT23_table5_iWUE_vcmaxgs.csv",
          row.names = FALSE)




