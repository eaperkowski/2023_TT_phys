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
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "invaded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")))
head(df)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, total_subplot, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:n_plantAvail_day),
               .funs = mean)


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

# Prep file for figure
nitrate_results <- cld(emmeans(nitrate, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = c("B", "B", "A", "A"))

ggplot() +
  geom_point(data = df.soil, 
             aes(x = canopy, y = nitrate_ppm_day, fill = gm.trt),
             position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2, shape = 21) +
  geom_errorbar(data = nitrate_results, aes(x = canopy, y = emmean,
                                            ymin = lower.CL, ymax = upper.CL,
                                            group = gm.trt),
                linewidth = 1, position = position_dodge(width = 0.75), width = 0.25) +
  geom_point(data = nitrate_results, aes(x = canopy, y = emmean, fill = gm.trt),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = nitrate_results, aes(x = canopy, y = 1, group = gm.trt,
                                        label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Canopy status", 
       y = expression(bold("[NO"["3"]*"-N] (ppm day"^"-1"*")")),
       fill = "Garlic mustard\n treatment") +
  theme_bw(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


##############################################################################
## Ammonium
##############################################################################
ammonium <- lmer(ammonium_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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

# Prep file for figure
ammonium_results <- cld(emmeans(ammonium, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = c("B", "B", "A", "A"))

ggplot() +
  geom_point(data = df.soil, 
             aes(x = canopy, y = ammonium_ppm_day, fill = gm.trt),
             position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2, shape = 21) +
  geom_errorbar(data = ammonium_results, aes(x = canopy, y = emmean,
                                            ymin = lower.CL, ymax = upper.CL,
                                            group = gm.trt),
                linewidth = 1, position = position_dodge(width = 0.75), width = 0.25) +
  geom_point(data = ammonium_results, aes(x = canopy, y = emmean, fill = gm.trt),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = ammonium_results, aes(x = canopy, y = 0.02, group = gm.trt,
                                        label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.005)) +
  labs(x = "Canopy status", 
       y = expression(bold("[NH"["4"]*"-N] (ppm day"^"-1"*")")),
       fill = "Garlic mustard\n treatment") +
  theme_bw(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

##############################################################################
## Phosphate
##############################################################################
phosphate <- lmer(phosphate_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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
phosphate_results <- cld(emmeans(phosphate, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both"))

phosphate_plot <- ggplot() +
  geom_point(data = df.soil, 
             aes(x = canopy, y = phosphate_ppm_day, fill = gm.trt),
             position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2, shape = 21) +
  geom_errorbar(data = phosphate_results, aes(x = canopy, y = emmean,
                                             ymin = lower.CL, ymax = upper.CL,
                                             group = gm.trt),
                linewidth = 1, position = position_dodge(width = 0.75), width = 0.25) +
  geom_point(data = phosphate_results, aes(x = canopy, y = emmean, fill = gm.trt),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = phosphate_results, aes(x = canopy, y = 0.06, group = gm.trt,
                                         label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 0.06), breaks = seq(0, 0.06, 0.02)) +
  labs(x = "Canopy status", 
       y = expression(bold("Phosphate concentration (ppm day"^"-1"*")")),
       fill = "Garlic mustard\n treatment") +
  theme_bw(base_size = 16) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


##############################################################################
## N availability (nitrate + ammonium)
##############################################################################
plant_availableN <- lmer(n_plantAvail_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

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
  mutate(.group = c("A", "A", "B", "B"))

n_plantavailable_plot <- ggplot() +
  geom_point(data = df.soil, 
             aes(x = canopy, y = n_plantAvail_day, fill = gm.trt),
             position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2, shape = 21) +
  geom_errorbar(data = n_results, aes(x = canopy, y = emmean,
                                              ymin = lower.CL, ymax = upper.CL,
                                              group = gm.trt),
                linewidth = 1, position = position_dodge(width = 0.75), width = 0.25) +
  geom_point(data = n_results, aes(x = canopy, y = emmean, fill = gm.trt),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = n_results, aes(x = canopy, y = 1, group = gm.trt,
                                          label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Canopy status", 
       y = expression(bold("Plant-available N (ppm day"^"-1"*")")),
       fill = "Garlic mustard\n treatment") +
  theme_bw(base_size = 16) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


## Write plot for soil phosphate and plant-available N
png("../figs/TT23_soil_nutrients.png", width = 12, height = 4.5,
    units = "in", res = 600)
ggarrange(n_plantavailable_plot, phosphate_plot, common.legend = TRUE,
          legend = "right", ncol = 2)
dev.off()


##############################################################################
## Anet - Tri
##############################################################################
df$anet[c(35, 80, 86, 116)] <- NA

anet.tri <- lmer(log(anet) ~ gm.trt * canopy * (phosphate_ppm_day + n_plantAvail_day) + 
               (1 | plot), data = subset(df, spp == "Tri"))

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
emmeans(anet.tri, pairwise~canopy)
cld(emmeans(anet.tri, pairwise~canopy*gm.trt, type = "response"))

##############################################################################
## Anet - Mai
##############################################################################
df$anet[c(49, 132)] <- NA

anet.mai <- lmer(log(anet) ~ gm.trt * canopy * (phosphate_ppm_day + n_plantAvail_day) + 
                   (1 | plot), data = subset(df, spp == "Mai"))

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
emmeans(anet.mai, pairwise~canopy)

##############################################################################
## Anet - Tri continuous garlic mustard
##############################################################################
anet.tri.cont <- lmer(log(anet) ~ total_subplot * canopy * (phosphate_ppm_day + n_plantAvail_day) + 
                   (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(anet.tri.cont)
qqnorm(residuals(anet.tri.cont))
qqline(residuals(anet.tri.cont))
densityPlot(residuals(anet.tri.cont))
shapiro.test(residuals(anet.tri.cont))
outlierTest(anet.tri.cont)

# Model output
summary(anet.tri.cont)
Anova(anet.tri.cont)
r.squaredGLMM(anet.tri.cont)

# Pairwise comparisons
test(emtrends(anet.tri.cont, ~n_plantAvail_day*canopy, "total_subplot",
              at = list(n_plantAvail_day = seq(0.1, 0.9, 0.2))))
emmeans(anet.tri.cont, pairwise~canopy)

##############################################################################
## Anet - Mai
##############################################################################
anet.mai.cont <- lmer(log(anet) ~ total_subplot * canopy * (phosphate_ppm_day + n_plantAvail_day) + 
                   (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(anet.mai.cont)
qqnorm(residuals(anet.mai.cont))
qqline(residuals(anet.mai.cont))
densityPlot(residuals(anet.mai.cont))
shapiro.test(residuals(anet.mai.cont))
outlierTest(anet.mai.cont)

# Model output
summary(anet.mai.cont)
Anova(anet.mai.cont)
r.squaredGLMM(anet.mai.cont)

# Pairwise comparisons
test(emtrends(anet.mai.cont, ~1, "total_subplot", type = "response"))
test(emtrends(anet.mai.cont, ~canopy, "total_subplot", type = "response"))


##############################################################################
## Vcmax regressed against GM treatment (cat.), canopy status (cat.),
## phosphate concentration (cont.), and plant-available N concentration (cont.)
##############################################################################
df$vcmax25[231] <- NA

vcmax.tri <- lmer(log(vcmax25) ~ 
                    gm.trt * canopy * (phosphate_ppm_day + n_plantAvail_day) + 
                    (1 | plot), data = subset(df, spp == "Tri"))


# Check model assumptions
plot(vcmax.tri)
qqnorm(residuals(vcmax.tri))
qqline(residuals(vcmax.tri))
densityPlot(residuals(vcmax.tri))
shapiro.test(residuals(vcmax.tri))
outlierTest(vcmax)


# Model output
summary(vcmax.tri)
Anova(vcmax.tri)
r.squaredGLMM(vcmax.tri)

# Pairwise comparisons
test(emtrends(vcmax.tri, pairwise~gm.trt*canopy, "n_plantAvail_day"))
cld(emmeans(vcmax.tri, pairwise~gm.trt*canopy, type = "response"))


##############################################################################
## Vcmax regressed against GM treatment (cat.), canopy status (cat.),
## phosphate concentration (cont.), and plant-available N concentration (cont.)
##############################################################################
vcmax.mai <- lmer(log(vcmax25) ~ 
                    gm.trt * canopy * (phosphate_ppm_day + n_plantAvail_day) + 
                    (1 | plot), data = subset(df, spp == "Mai"))


# Check model assumptions
plot(vcmax.mai)
qqnorm(residuals(vcmax.mai))
qqline(residuals(vcmax.mai))
densityPlot(residuals(vcmax.mai))
shapiro.test(residuals(vcmax.mai))
outlierTest(vcmax)

# Model output
summary(vcmax.mai)
Anova(vcmax.mai)
r.squaredGLMM(vcmax.mai)

# Pairwise comparisons
cld(emmeans(vcmax.mai, pairwise~gm.trt*canopy))
test(emtrends(vcmax.mai, pairwise~gm.trt, "n_plantAvail_day"))
test(emtrends(vcmax.mai, pairwise~canopy, "phosphate_ppm_day"))

# Single factor patterns
emmeans(vcmax.mai, pairwise~canopy)

##############################################################################
## Vcmax regressed against garlic mustard density
##############################################################################
vcmax.cont.tri <- lmer(log(vcmax25) ~ 
                         total_subplot * canopy * (phosphate_ppm_day + n_plantAvail_day) 
              + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(vcmax.cont.tri)
qqnorm(residuals(vcmax.cont.tri))
qqline(residuals(vcmax.cont.tri))
densityPlot(residuals(vcmax.cont.tri))
shapiro.test(residuals(vcmax.cont.tri))
outlierTest(vcmax.cont.tri)

# Model output
summary(vcmax.cont.tri)
Anova(vcmax.cont.tri)
r.squaredGLMM(vcmax.cont.tri)

# Post-hoc comparisons
test(emtrends(vcmax.cont.tri, ~n_plantAvail_day*canopy, "total_subplot", 
              at = list(n_plantAvail_day = seq(0.1, 0.9, 0.2))))

# Individual trends
emmeans(vcmax.cont.tri, pairwise~canopy)

##############################################################################
## Vcmax regressed against garlic mustard density
##############################################################################
vcmax.cont.tri <- lmer(log(vcmax25) ~ total_subplot * canopy * (phosphate_ppm_day + n_plantAvail_day) 
                       + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(vcmax.cont.tri)
qqnorm(residuals(vcmax.cont.tri))
qqline(residuals(vcmax.cont.tri))
densityPlot(residuals(vcmax.cont.tri))
shapiro.test(residuals(vcmax.cont.tri))
outlierTest(vcmax.cont.tri)

# Model output
summary(vcmax.cont.tri)
Anova(vcmax.cont.tri)
r.squaredGLMM(vcmax.cont.tri)

# Post-hoc comparisons
test(emtrends(vcmax.cont.tri, ~n_plantAvail_day*canopy, "total_subplot", 
              at = list(n_plantAvail_day = seq(0.1, 0.9, 0.1))))

# Individual trends
emmeans(vcmax.cont.tri, pairwise~canopy)




