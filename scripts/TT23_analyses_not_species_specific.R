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
         vcmax.gs = vcmax25 / gsw,
         spad.gs = SPAD / gsw)
head(df)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

##############################################################################
## Anet
##############################################################################
anet <- lmer(
  sqrt(anet) ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(anet)
qqnorm(residuals(anet))
qqline(residuals(anet))
densityPlot(residuals(anet))
shapiro.test(residuals(anet))
outlierTest(anet)

# Model output
summary(anet)
Anova(anet)
r.squaredGLMM(anet)

# Pairwise comparisons
cld(emmeans(anet, pairwise~gm.trt*spp))
cld(emmeans(anet, pairwise~gm.trt*canopy))
cld(emmeans(anet, pairwise~canopy*spp))

##############################################################################
## Stomatal conductance 
##############################################################################
df$gsw[c(53, 91)] <- NA

gsw <- lmer(
  gsw ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(gsw)
qqnorm(residuals(gsw))
qqline(residuals(gsw))
densityPlot(residuals(gsw))
shapiro.test(residuals(gsw))
outlierTest(gsw)

# Model output
summary(gsw)
Anova(gsw)
r.squaredGLMM(gsw)

# Pairwise comparisons
cld(emmeans(gsw, pairwise~canopy*spp))
cld(emmeans(gsw, pairwise~gm.trt*spp))

##############################################################################
## Stomatal limitation 
##############################################################################
l <- lmer(
 log(l) ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(l)
qqnorm(residuals(l))
qqline(residuals(l))
densityPlot(residuals(l))
shapiro.test(residuals(l))
outlierTest(l)

# Model output
summary(l)
Anova(l)
r.squaredGLMM(l)

# Pairwise comparisons
cld(emmeans(l, pairwise~gm.trt*canopy*spp))
cld(emmeans(l, pairwise~canopy*spp))
cld(emmeans(l, pairwise~gm.trt*spp))

##############################################################################
## SPAD
##############################################################################
spad <- lmer(
  SPAD ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(spad)
qqnorm(residuals(spad))
qqline(residuals(spad))
densityPlot(residuals(spad))
shapiro.test(residuals(spad))
outlierTest(spad)

# Model output
summary(spad)
Anova(spad)
r.squaredGLMM(spad)

# Pairwise comparisons
cld(emmeans(spad, pairwise~canopy*spp))

##############################################################################
## Vcmax
##############################################################################
df$vcmax25[c(225)] <- NA

vcmax <- lmer(
  log(vcmax25) ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(vcmax)
qqnorm(residuals(vcmax))
qqline(residuals(vcmax))
densityPlot(residuals(vcmax))
shapiro.test(residuals(vcmax))
outlierTest(vcmax)

# Model output
summary(vcmax)
Anova(vcmax)
r.squaredGLMM(vcmax)

# Pairwise comparisons
cld(emmeans(vcmax, pairwise~canopy*spp, type = "response"))

##############################################################################
## Jmax
##############################################################################
df$jmax25[c(27, 49)] <- NA

jmax <- lmer(
  jmax25 ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(jmax)
qqnorm(residuals(jmax))
qqline(residuals(jmax))
densityPlot(residuals(jmax))
shapiro.test(residuals(jmax))
outlierTest(jmax)

# Model output
summary(jmax)
Anova(jmax)
r.squaredGLMM(jmax)

# Pairwise comparisons
cld(emmeans(jmax, pairwise~canopy*spp, type = "response"))


##############################################################################
## Jmax:Vcmax
##############################################################################
df$jmax.vcmax[c(181)] <- NA

jmax.vcmax <- lmer(
  jmax.vcmax ~ gm.trt * canopy * spp + (1 | plot), data = subset(df, spp != "Ari"))

# Check model assumptions
plot(jmax.vcmax)
qqnorm(residuals(jmax.vcmax))
qqline(residuals(jmax.vcmax))
densityPlot(residuals(jmax.vcmax))
shapiro.test(residuals(jmax.vcmax))
outlierTest(jmax.vcmax)

# Model output
summary(jmax.vcmax)
Anova(jmax.vcmax)
r.squaredGLMM(jmax.vcmax)

# Pairwise comparisons
cld(emmeans(jmax, pairwise~canopy*spp, type = "response"))

