#####################################################################
# Libraries and custom functions
#####################################################################
# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(MuMIn)
library(multcomp)

# Read compiled trait dataset
photo_traits <- read.csv("../data/TT23_phys_data.csv") %>%
  mutate(date = mdy(date),
         doy_photo = yday(date)) %>%
  dplyr::select(id:gm.trt, date_photo = date, doy_photo, canopy:inorg_n_ppm)

# How many measurements per ID?
n_measurements <- photo_traits %>%
  group_by(id, spp, plot, subplot, gm.trt) %>%
  summarize(n_measurements = length(id)) %>%
  ungroup() %>%
  dplyr::select(id, spp, n_meas = n_measurements)
head(n_measurements)

# Join n measurements into photo traits
photo_traits_analysis <- photo_traits %>%
  left_join(n_measurements, by = "id") %>%
  filter(n_meas > 1) %>%
  dplyr::select(id, spp = spp.x, plot:inorg_n_ppm, n_meas)

##############################################################################
## Anet - Tri
##############################################################################
anet.tri <- lmer(log(anet) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                 data = subset(photo_traits_analysis, spp == "Tri"))

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
cld(emmeans(anet.tri, pairwise~canopy*gm.trt, type = "response"))
emmeans(anet.tri, pairwise~canopy, type = "response")
emmeans(anet.tri, pairwise~gm.trt, type = "response")

##############################################################################
## gsw - Tri
##############################################################################
gsw.tri <- lmer(log(gsw) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                data = subset(photo_traits_analysis, spp == "Tri"))

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

##############################################################################
## stomatal limitation - Tri
##############################################################################
l.tri <- lmer(l ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                data = subset(photo_traits_analysis, spp == "Tri"))

# Check model assumptions
plot(l.tri)
qqnorm(residuals(l.tri))
qqline(residuals(l.tri))
densityPlot(residuals(l.tri))
shapiro.test(residuals(l.tri))
outlierTest(l.tri)

# Model output
summary(l.tri)
Anova(l.tri)
r.squaredGLMM(l.tri)

##############################################################################
## Vcmax25 - Tri
##############################################################################
photo_traits_analysis$vcmax25[c(67, 96)] <- NA

vcmax25.tri <- lmer(log(vcmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                    data = subset(photo_traits_analysis, spp == "Tri"))

# Check model assumptions
plot(vcmax25.tri)
qqnorm(residuals(vcmax25.tri))
qqline(residuals(vcmax25.tri))
densityPlot(residuals(vcmax25.tri))
shapiro.test(residuals(vcmax25.tri))
outlierTest(vcmax25.tri)

# Model output
summary(vcmax25.tri)
Anova(vcmax25.tri)
r.squaredGLMM(vcmax25.tri)

# Pairwise comparisons
cld(emmeans(vcmax25.tri, pairwise~gm.trt*canopy, type = "response"))
emmeans(vcmax25.tri, pairwise~gm.trt, type = "response")
emmeans(vcmax25.tri, pairwise~canopy, type = "response")

##############################################################################
## Jmax25 - Tri
##############################################################################
photo_traits_analysis$jmax25[c(67, 96)] <- NA

jmax25.tri <- lmer(log(jmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                    data = subset(photo_traits_analysis, spp == "Tri"))

# Check model assumptions
plot(jmax25.tri)
qqnorm(residuals(jmax25.tri))
qqline(residuals(jmax25.tri))
densityPlot(residuals(jmax25.tri))
shapiro.test(residuals(jmax25.tri))
outlierTest(jmax25.tri)

# Model output
summary(jmax25.tri)
Anova(jmax25.tri)
r.squaredGLMM(jmax25.tri)

# Pairwise comparisons
cld(emmeans(jmax25.tri, pairwise~gm.trt*canopy, type = "response"))
emmeans(jmax25.tri, pairwise~gm.trt, type = "response")
emmeans(jmax25.tri, pairwise~canopy, type = "response")

##############################################################################
## Jmax25:Vcmax25 - Tri
##############################################################################
photo_traits_analysis$jmax.vcmax[c(97)] <- NA

jmax25_vcmax25.tri <- lmer(jmax.vcmax ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                           data = subset(photo_traits_analysis, spp == "Tri"))

# Check model assumptions
plot(jmax25_vcmax25.tri)
qqnorm(residuals(jmax25_vcmax25.tri))
qqline(residuals(jmax25_vcmax25.tri))
densityPlot(residuals(jmax25_vcmax25.tri))
shapiro.test(residuals(jmax25_vcmax25.tri))
outlierTest(jmax25_vcmax25.tri)

# Model output
summary(jmax25_vcmax25.tri)
Anova(jmax25_vcmax25.tri)
r.squaredGLMM(jmax25_vcmax25.tri)

# Pairwise comparisons
emmeans(jmax25_vcmax25.tri, pairwise~gm.trt, type = "response")
emmeans(jmax25_vcmax25.tri, pairwise~canopy, type = "response")

##############################################################################
## SPAD - Tri
##############################################################################
spad.tri <- lmer(SPAD ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                  data = subset(photo_traits_analysis, spp == "Tri"))

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
cld(emmeans(spad.tri, pairwise~gm.trt*canopy, type = "response"))

emmeans(spad.tri, pairwise~gm.trt, type = "response")
emmeans(spad.tri, pairwise~canopy, type = "response")



##############################################################################
## Anet - Mai
##############################################################################
photo_traits_analysis$anet[43] <- NA

anet.mai <- lmer(anet ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                 data = subset(photo_traits_analysis, spp == "Mai"))

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

##############################################################################
## gsw - Mai
##############################################################################
photo_traits_analysis$gsw[43] <- NA

gsw.mai <- lmer(gsw ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                data = subset(photo_traits_analysis, spp == "Mai"))

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

##############################################################################
## stomatal limitation - Mai
##############################################################################
photo_traits_analysis$l[68] <- NA

l.mai <- lmer(l ~ gm.trt * canopy + (1 | plot) + (1 | id), 
              data = subset(photo_traits_analysis, spp == "Mai"))

# Check model assumptions
plot(l.mai)
qqnorm(residuals(l.mai))
qqline(residuals(l.mai))
densityPlot(residuals(l.mai))
shapiro.test(residuals(l.mai))
outlierTest(l.mai)

# Model output
summary(l.mai)
Anova(l.mai)
r.squaredGLMM(l.mai)

##############################################################################
## Vcmax25 - Mai
##############################################################################
vcmax25.mai <- lmer(log(vcmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                    data = subset(photo_traits_analysis, spp == "Mai"))

# Check model assumptions
plot(vcmax25.mai)
qqnorm(residuals(vcmax25.mai))
qqline(residuals(vcmax25.mai))
densityPlot(residuals(vcmax25.mai))
shapiro.test(residuals(vcmax25.mai))
outlierTest(vcmax25.mai)

# Model output
summary(vcmax25.mai)
Anova(vcmax25.mai)
r.squaredGLMM(vcmax25.mai)

# Pairwise comparisons
emmeans(vcmax25.mai, pairwise~canopy, type = "response")

##############################################################################
## Jmax25 - Tri
##############################################################################
jmax25.mai <- lmer(log(jmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                   data = subset(photo_traits_analysis, spp == "Mai"))

# Check model assumptions
plot(jmax25.mai)
qqnorm(residuals(jmax25.mai))
qqline(residuals(jmax25.mai))
densityPlot(residuals(jmax25.mai))
shapiro.test(residuals(jmax25.mai))
outlierTest(jmax25.mai)

# Model output
summary(jmax25.mai)
Anova(jmax25.mai)
r.squaredGLMM(jmax25.mai)

# Pairwise comparisons
emmeans(jmax25.mai, pairwise~canopy, type = "response")

##############################################################################
## Jmax25:Vcmax25 - Tri
##############################################################################
jmax25_vcmax25.mai <- lmer(jmax.vcmax ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                           data = subset(photo_traits_analysis, spp == "Mai"))

# Check model assumptions
plot(jmax25_vcmax25.mai)
qqnorm(residuals(jmax25_vcmax25.mai))
qqline(residuals(jmax25_vcmax25.mai))
densityPlot(residuals(jmax25_vcmax25.mai))
shapiro.test(residuals(jmax25_vcmax25.mai))
outlierTest(jmax25_vcmax25.mai)

# Model output
summary(jmax25_vcmax25.mai)
Anova(jmax25_vcmax25.mai)
r.squaredGLMM(jmax25_vcmax25.mai)

##############################################################################
## SPAD - Tri
##############################################################################
spad.mai <- lmer(SPAD ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                 data = subset(photo_traits_analysis, spp == "Mai"))

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
emmeans(spad.mai, pairwise~canopy, type = "response")