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
## Nitrate
##############################################################################
nitrate <- lmer(nitrate_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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

##############################################################################
## Ammonium
##############################################################################
ammonium <- lmer(
  ammonium_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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
emmeans(ammonium, pairwise~canopy)

##############################################################################
## Phosphate
##############################################################################
phosphate <- lmer(
  phosphate_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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

# Prep file for figure
phosphate_results <- cld(emmeans(phosphate, pairwise~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both"))

phosphate_plot <- ggplot() +
  geom_point(data = df.soil, 
             aes(x = gm.trt, y = phosphate_ppm_day, fill = canopy),
             position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_errorbar(data = phosphate_results, aes(x = gm.trt, y = emmean,
                                             ymin = lower.CL, ymax = upper.CL,
                                             group = canopy),
                linewidth = 1, position = position_dodge(width = 0.75), width = 0.25) +
  geom_point(data = phosphate_results, aes(x = gm.trt, y = emmean, fill = canopy),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = phosphate_results, aes(x = gm.trt, y = 0.06, group = canopy,
                                         label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_fill_discrete(labels = c("Open", "Closed")) +
  scale_x_discrete(labels = c("Ambient", "Weeded")) +
  scale_y_continuous(limits = c(0, 0.06), breaks = seq(0, 0.06, 0.02)) +
  labs(x = expression(paste(bolditalic("Alliaria"), bold(" treatment"))),
       y = expression(bold("Soil P availability (ppm day"^"-1"*")")),
       fill = "Canopy") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
phosphate_plot

##############################################################################
## N availability (nitrate + ammonium)
##############################################################################
plant_availableN <- lmer(
  n_plantAvail_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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

# Prep file for figure
n_results <- cld(emmeans(plant_availableN, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = c("B", "B", "A", "A"))

n_plantavailable_plot <- ggplot() +
  geom_point(data = df.soil, 
             aes(x = gm.trt, y = n_plantAvail_day, fill = canopy),
             position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_errorbar(data = n_results, aes(x = gm.trt, y = emmean,
                                              ymin = lower.CL, ymax = upper.CL,
                                              group = canopy),
                linewidth = 1, position = position_dodge(width = 0.75), width = 0.25) +
  geom_point(data = n_results, aes(x = gm.trt, y = emmean, fill = canopy),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = n_results, 
            aes(x = gm.trt, y = 1, group = canopy, label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("Ambient", "Weeded")) +
  scale_fill_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(paste(bolditalic("Alliaria"), bold(" treatment"))),
       y = expression(bold("Soil N availability (ppm day"^"-1"*")")),
       fill = "Canopy") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
n_plantavailable_plot

## Write plot for soil phosphate and plant-available N
png("../figs/TT23_soil_nutrients.png", width = 12, height = 5,
    units = "in", res = 600)
ggarrange(n_plantavailable_plot, phosphate_plot, common.legend = TRUE,
          legend = "right", ncol = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 18))
dev.off()


##############################################################################
## Anet - Tri
##############################################################################
df$anet[c(35, 80, 86)] <- NA

anet.tri <- lmer(
  log(anet) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
cld(emmeans(anet.tri, pairwise~gm.trt*canopy, type = "response"))

##############################################################################
## Anet - Mai
##############################################################################
df$anet[c(49, 132, 145)] <- NA

anet.mai <- lmer(
  log(anet) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(anet.mai, pairwise~gm.trt)
cld(emmeans(anet.mai, pairwise~gm.trt*canopy, type = "response"))

##############################################################################
## gs - Tri
##############################################################################
df$gsw[c(86)] <- NA

gsw.tri <- lmer(
 gsw ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

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
  gsw ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

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

##############################################################################
## stomatal limitation - Tri
##############################################################################
df$stom.lim[c(228, 229)] <- NA

stomlim.tri <- lmer(
  log(stom.lim) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri" & stom.lim > 0))

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


##############################################################################
## stomatal limitation - Mai
##############################################################################
stom.lim.mai <- lmer(
  log(stom.lim) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai" & stom.lim > 0))

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
emmeans(stom.lim.mai, pairwise~canopy)
emmeans(stom.lim.mai, pairwise~gm.trt*canopy, type = "response")

##############################################################################
## Vcmax - Tri
##############################################################################
df$vcmax25[183] <- NA

vcmax.tri <- lmer(log(vcmax25) ~ 
                    canopy * (phosphate_ppm_day + n_plantAvail_day) +
                    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))


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
  log(vcmax25) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))


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

jmax.tri <- lmer(log(jmax25) ~ 
                    canopy * (phosphate_ppm_day + n_plantAvail_day) +
                    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))


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
  jmax25 ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))


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
df$jmax.vcmax[184] <- NA

jmax.vcmax.tri <- lmer(log(jmax.vcmax) ~ 
                   canopy * (phosphate_ppm_day + n_plantAvail_day) +
                   gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))


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

jmax.vcmax.mai <- lmer(log(jmax.vcmax) ~ 
                         canopy * (phosphate_ppm_day + n_plantAvail_day) +
                         gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))


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
df$iwue[c(80, 86, 228)] <- NA

iwue.tri <- lmer(log(iwue) ~ 
                         canopy * (phosphate_ppm_day + n_plantAvail_day) +
                         gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))


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

iwue.mai <- lmer(log(iwue) ~ 
                   canopy * (phosphate_ppm_day + n_plantAvail_day) +
                   gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))


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
df$iwue[c(80, 86, 228)] <- NA

vcmax.gs.tri <- lmer(log(vcmax.gs) ~ 
                   canopy * (phosphate_ppm_day + n_plantAvail_day) +
                   gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))


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
vcmax.gs.mai <- lmer(log(vcmax.gs) ~ 
                       canopy * (phosphate_ppm_day + n_plantAvail_day) +
                       gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))


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
