##############################################################################
## Prepare libraries  
##############################################################################
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

## Remove outliers
df$anet[c(35, 49, 80, 86, 132, 145)] <- NA
df$gsw[c(86, 120)] <- NA
df$stom.lim[c(228, 229)] <- NA
df$vcmax25[c(183, 231)] <- NA
df$jmax25[183] <- NA
df$jmax.vcmax[c(94, 184, 225, 231)] <- NA
df$iwue[c(80, 86, 104, 228)] <- NA

## Create models for soil data
nitrate <- lmer(
  nitrate_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)
ammonium <- lmer(
  ammonium_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)
phosphate <- lmer(
  phosphate_ppm_day ~ gm.trt * canopy + (1 | plot), data = df.soil)
plant_availableN <- lmer(
  n_plantAvail_day ~ gm.trt * canopy + (1 | plot), data = df.soil)

## Create models for photosynthesis data
anet.tri <- lmer(
  log(anet) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
anet.mai <- lmer(
  log(anet) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))
gsw.tri <- lmer(
  gsw ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
gsw.mai <- lmer(
  gsw ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))
stomlim.tri <- lmer(
  log(stom.lim) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri" & stom.lim > 0))
stom.lim.mai <- lmer(
  log(stom.lim) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai" & stom.lim > 0))
vcmax.tri <- lmer(
  log(vcmax25) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
vcmax.mai <- lmer(
  log(vcmax25) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))
jmax.tri <- lmer(
  log(jmax25) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
jmax.mai <- lmer(
  jmax25 ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))
jmax.vcmax.tri <- lmer(
  log(jmax.vcmax) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
jmax.vcmax.mai <- lmer(
  log(jmax.vcmax) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))
iwue.tri <- lmer(
  log(iwue) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
iwue.mai <- lmer(
  log(iwue) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))
vcmax.gs.tri <- lmer(
  log(vcmax.gs) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))
vcmax.gs.mai <- lmer(
  log(vcmax.gs) ~ canopy * (phosphate_ppm_day + n_plantAvail_day) +
    gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))


##############################################################################
## Soil N availability 
##############################################################################
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

##############################################################################
## Phosphate figure  
##############################################################################
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
## Net photosynthesis - Trillium
##############################################################################
anet_results <- cld(emmeans(anet.tri, ~gm.trt*canopy, type = "response"), 
                Letters = LETTERS) %>% 
  data.frame() %>%
  mutate(.group = c("B", "C", "A", "A")) %>% data.frame()


ggplot(data = subset(df, spp == "Tri")) +
  geom_point(aes(x = canopy, y = anet, fill = gm.trt),
       position = position_jitterdodge(dodge.width = 0.75,
                                       jitter.width = 0.1),
       alpha = 0.5, size = 2.5, shape = 21) +
  geom_errorbar(data = anet_results, 
                aes(x = canopy, y = response, ymin = lower.CL, ymax = upper.CL,
                    group = gm.trt), linewidth = 1, 
                position = position_dodge(width = 0.75), 
                width = 0.25) +
  geom_point(data = anet_results, aes(x = canopy, y = response, fill = gm.trt),
             size = 5, position = position_dodge(width = 0.75), shape = 21) +
  geom_text(data = anet_results, aes(x = canopy, y = 20, group = gm.trt,
                                          label = .group),
            position = position_dodge(width = 0.75), fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_fill_discrete(labels = c("Ambient", "Weeded")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(fill = expression(paste(bolditalic("Alliaria"), bold(" treatment"))),
       y = expression(bold(italic("A")["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       x = "Canopy") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25))
  




##############################################################################
## Figure 1: Soil nutrients 
##############################################################################
png("../figs/TT23_soil_nutrients.png", width = 12, height = 5,
    units = "in", res = 600)
ggarrange(n_plantavailable_plot, phosphate_plot, common.legend = TRUE,
          legend = "right", ncol = 2, labels = c("(a)", "(b)"),
          font.label = list(size = 18))
dev.off()








