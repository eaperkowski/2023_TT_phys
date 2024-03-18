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
library(gghalves)

## Read compiled data file
df <- read.csv("../data/TT23_compiled_datasheet.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "invaded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")),
         spp = factor(spp, levels = c("Tri", "Mai", "Ari")),
         vcmax.gs = vcmax25 / gsw,
         spad.gs = SPAD / gsw) %>%
  unite(col = "gm.canopy", gm.trt, canopy, sep = "_", remove = FALSE) %>%
  mutate(gm.canopy = factor(gm.canopy, levels = c("weeded_pre_closure",
                                                  "invaded_pre_closure",
                                                  "weeded_post_closure",
                                                  "invaded_post_closure")))
head(df)

# helper fxn to change "NaN" to "NA"
NaN_to_NA <- function(x) ifelse(is.nan(x), NA, x)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:inorg_n_ppm),
               .funs = mean, na.rm = TRUE) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "invaded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure"))) %>%
  mutate(across(nitrate_ppm:inorg_n_ppm, .fns = NaN_to_NA),
         np.ratio = inorg_n_ppm/phosphate_ppm)

## Remove outliers
df.soil$ammonium_ppm[c(40, 56)] <- NA
df.soil$np.ratio[c(55)] <- NA
df$anet[c(35, 71, 80, 86, 116, 120)] <- NA
df$gsw[c(86, 120)] <- NA
df$stom.lim[c(228, 229)] <- NA
df$vcmax25[c(183, 231)] <- NA
df$jmax25[183] <- NA
df$jmax.vcmax[c(94, 184, 224, 225, 231)] <- NA
df$SPAD[c(122)] <- NA
df$Phi2[92] <- NA
df$iwue[c(80, 86, 104)] <- NA
df$spad.gs[c(49, 132)] <- NA

## Create models for soil data
nitrate <- lmer(nitrate_ppm ~ gm.trt * canopy + (1 | plot), 
                data = df.soil)
ammonium <- lmer(ammonium_ppm ~ gm.trt * canopy + (1 | plot), 
                 data = df.soil)
phosphate <- lmer(phosphate_ppm ~ gm.trt * canopy + (1 | plot), 
                  data = df.soil)
plant_availableN <- lmer(inorg_n_ppm ~ gm.trt * canopy + (1 | plot), 
                         data = df.soil)
n_to_p_ratio <- lmer(np.ratio ~ gm.trt * canopy + (1 | plot), 
                     data = df.soil)

## Create models for photosynthesis data
anet.tri <- lmer(anet ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Tri"))
anet.mai <- lmer(anet ~ gm.trt * canopy  + (1 | plot),
                 data = subset(df, spp == "Mai"))
gsw.tri <- lmer(gsw ~ gm.trt * canopy + (1 | plot),
                data = subset(df, spp == "Tri"))
gsw.mai <- lmer(gsw ~ gm.trt * canopy + (1 | plot),
                data = subset(df, spp == "Mai"))
stomlim.tri <- lmer(log(stom.lim) ~ gm.trt * canopy  + (1 | plot),
                    data = subset(df, spp == "Tri" & stom.lim > 0))
stom.lim.mai <- lmer(log(stom.lim) ~ gm.trt * canopy + (1 | plot),
                     data = subset(df, spp == "Mai" & stom.lim > 0))
vcmax.tri <- lmer(log(vcmax25) ~ gm.trt * canopy + (1 | plot),
                  data = subset(df, spp == "Tri"))
vcmax.mai <- lmer(log(vcmax25) ~ gm.trt * canopy  + (1 | plot),
                  data = subset(df, spp == "Mai"))
jmax.tri <- lmer(log(jmax25) ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Tri"))
jmax.mai <- lmer(jmax25 ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Mai"))
jmax.vcmax.tri <- lmer(log(jmax.vcmax) ~ gm.trt * canopy + (1 | plot),
                       data = subset(df, spp == "Tri"))
jmax.vcmax.mai <- lmer(log(jmax.vcmax) ~ gm.trt * canopy + (1 | plot),
                       data = subset(df, spp == "Mai"))
spad.tri <- lmer(SPAD ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Tri"))
spad.mai <- lmer(SPAD ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Mai"))
phips2.tri <- lmer(Phi2 ~ gm.trt * canopy + (1 | plot), 
                   data = subset(df, spp == "Tri"))
phips2.mai <- lmer(Phi2 ~ gm.trt * canopy + (1 | plot), 
                   data = subset(df, spp == "Mai"))
iwue.tri <- lmer(log(iwue) ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Tri"))
iwue.mai <- lmer(log(iwue) ~ gm.trt * canopy  + (1 | plot),
                 data = subset(df, spp == "Mai"))
vcmax.gs.tri <- lmer(log(vcmax.gs) ~ gm.trt * canopy + (1 | plot),
                     data = subset(df, spp == "Tri"))
vcmax.gs.mai <- lmer(log(vcmax.gs) ~ gm.trt * canopy + (1 | plot),
                     data = subset(df, spp == "Mai"))
spad.gs.tri <- lmer(log(spad.gs) ~ gm.trt * canopy + (1 | plot),
                    data = subset(df, spp == "Tri"))
spad.gs.mai <- lmer(log(spad.gs) ~ gm.trt * canopy + (1 | plot),
                    data = subset(df, spp == "Mai"))

## Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

## Color palettes
gm.colors <- c("#F7FCB9", "#D95F0E")
canopy.colors <- c("#EDF8FB", "#3182bd")

##############################################################################
## Soil N availability 
##############################################################################
# Prep file for gm.trt figure
inorgN_gmtrt_results <- cld(emmeans(plant_availableN, pairwise~gm.trt), 
                               Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

nitrogen_gmtrt_plot <- ggplot(data = df.soil,
                              aes(x = gm.trt, y = inorg_n_ppm, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = inorgN_gmtrt_results, 
            aes(x = gm.trt, y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "Soil inorg. N (ppm)") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
nitrogen_gmtrt_plot

# Prep file for canopy figure
inorgN_canopy_results <- cld(emmeans(plant_availableN, pairwise~canopy),
                             Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"))

nitrogen_canopy_plot <- ggplot(data = df.soil,
                               aes(x = canopy, y = inorg_n_ppm, fill = canopy)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = inorgN_canopy_results, 
            aes(x = canopy, y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = canopy.colors) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = "Measurement period",
       y = "Soil inorg. N (ppm)") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
nitrogen_canopy_plot

##############################################################################
## Phosphate figure  
##############################################################################
# Prep file for gm.trt figure
phosphate_gmtrt_results <- cld(emmeans(phosphate, pairwise~gm.trt), 
                         Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

phosphate_gmtrt_plot <- ggplot(data = df.soil,
                         aes(x = gm.trt, y = phosphate_ppm, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = phosphate_gmtrt_results, 
            aes(x = gm.trt, y = 2, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "Soil phosphate (ppm)") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
phosphate_gmtrt_plot

# Prep file for canopy figure
phosphate_canopy_results <- cld(emmeans(phosphate, pairwise~canopy), 
                               Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"))

phosphate_canopy_plot <- ggplot(data = df.soil,
                               aes(x = canopy, y = phosphate_ppm, fill = canopy)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = phosphate_canopy_results, 
            aes(x = canopy, y = 2, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = canopy.colors) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = "Measurement period",
       y = "Soil phosphate (ppm)") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
phosphate_canopy_plot

##############################################################################
## Soil N:P figure  
##############################################################################
Anova(n_to_p_ratio)

# Prep file for gm.trt figure
soil_np_results_gmtrt <- cld(emmeans(n_to_p_ratio, pairwise~gm.trt), 
                         Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both"))

soil_np_gmtrt_plot <- ggplot(data = df.soil,
                       aes(x = gm.trt, y = np.ratio, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = soil_np_results_gmtrt, 
            aes(x = gm.trt, y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "Soil N:P ratio (unitless)") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
soil_np_gmtrt_plot

# Prep file for measurement period figure
soil_np_results_canopy <- cld(emmeans(n_to_p_ratio, pairwise~canopy), 
                             Letters = LETTERS, reversed = TRUE)

soil_np_canopy_plot <- ggplot(data = df.soil,
                             aes(x = canopy, y = np.ratio, fill = canopy)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = soil_np_results_canopy, 
            aes(x = canopy, y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = canopy.colors) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = "Measurement period",
       y = "Soil N:P ratio (unitless)") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
soil_np_canopy_plot

##############################################################################
## Net photosynthesis - Trillium
##############################################################################
Anova(anet.tri)

anet_tri_results <- cld(emmeans(anet.tri, ~gm.trt, type = "response"), 
                Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.125) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_tri_results, 
            aes(x = gm.trt, y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = expression(bold(italic("A")["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
anet_tri_plot

##############################################################################
## Net photosynthesis - Maianthemum
##############################################################################
Anova(anet.mai)

anet_mai_results <- cld(emmeans(anet.mai, ~gm.trt, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                        aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_mai_results, 
            aes(x = gm.trt, y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = expression(bold(italic("A")["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
anet_mai_plot

##############################################################################
## Stomatal conductance - Tri
##############################################################################
Anova(gsw.tri)

gsw_tri_results <- cld(emmeans(gsw.tri, ~gm.trt, type = "response"), 
                        Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

gsw_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                  aes(x = gm.trt, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = gsw_tri_results, 
            aes(x = gm.trt, y = 0.3, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = expression(bold(italic("g")["sw"]*" (mol m"^"-2"*" s"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
gsw_tri_plot

##############################################################################
## Stomatal conductance - Maianthemum
##############################################################################
Anova(gsw.mai)

gsw_mai_results <- cld(emmeans(gsw.mai, ~gm.trt, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

gsw_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                       aes(x = gm.trt, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = gsw_mai_results, 
            aes(x = gm.trt, y = 0.3, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = expression(bold(italic("g")["sw"]*" (mol m"^"-2"*" s"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
gsw_mai_plot

##############################################################################
## Stomatal limitation - Tri
##############################################################################
Anova(stomlim.tri)

stomlim_tri_results <- cld(emmeans(stomlim.tri, ~gm.trt, type = "response"), 
                       Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

stomlim_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                           aes(x = gm.trt, y = stom.lim, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = stomlim_tri_results, 
            aes(x = gm.trt, y = 1, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "Stom. limitation (unitless)") +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
stomlim_tri_plot

##############################################################################
## Stomatal limitation - Mai
##############################################################################
Anova(stom.lim.mai)

stomlim_mai_results <- cld(emmeans(stom.lim.mai, ~gm.trt, type = "response"), 
                       Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

stomlim_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                           aes(x = gm.trt, y = stom.lim, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = stomlim_mai_results, 
            aes(x = gm.trt, y = 1, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "Stom. limitation (unitless)") +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
stomlim_mai_plot

##############################################################################
## Vcmax - Tri
##############################################################################
Anova(vcmax.tri)

vcmax_tri_results <- cld(emmeans(vcmax.tri, pairwise~gm.trt*canopy, type = "response"), 
                           Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

vcmax_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                           aes(x = canopy, y = vcmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = vcmax_tri_results, 
            aes(x = canopy, y = 200, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
vcmax_tri_plot

##############################################################################
## Vcmax - Mai
##############################################################################
Anova(vcmax.mai)

vcmax_mai_results <- cld(emmeans(vcmax.mai, ~gm.trt*canopy, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

vcmax_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                         aes(x = canopy, y = vcmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = vcmax_mai_results, 
            aes(x = canopy, y = 100, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
vcmax_mai_plot

##############################################################################
## Jmax - Tri
##############################################################################
jmax_tri_results <- cld(emmeans(jmax.tri, ~gm.trt*canopy, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jmax_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                         aes(x = canopy, y = jmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jmax_tri_results, 
            aes(x = canopy, y = 300, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 100)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jmax_tri_plot

##############################################################################
## Jmax - Mai
##############################################################################
Anova(jmax.mai)

jmax_mai_results <- cld(emmeans(jmax.mai, ~gm.trt*canopy, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jmax_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                         aes(x = canopy, y = jmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jmax_mai_results, 
            aes(x = canopy, y = 150, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jmax_mai_plot

##############################################################################
## Jmax: Vcmax - Tri
##############################################################################
jvmax_tri_results <- cld(emmeans(jmax.vcmax.tri, ~gm.trt*canopy, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jvmax_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = canopy, y = jmax.vcmax, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jvmax_tri_results, 
            aes(x = canopy, y = 0.7, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0.4, 0.70), breaks = seq(0.4, 0.70, 0.1)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"]*" (unitless)")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jvmax_tri_plot

##############################################################################
## Jmax:Vcmax - Mai
##############################################################################
jvmax_mai_results <- cld(emmeans(jmax.vcmax.mai, ~gm.trt*canopy, type = "response"), 
                        Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jvmax_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                        aes(x = canopy, y = jmax.vcmax, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jvmax_mai_results, 
            aes(x = canopy, y = 0.7, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0.4, 0.7), breaks = seq(0.4, 0.7, 0.1)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"]*" (unitless)")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jvmax_mai_plot

##############################################################################
## SPAD - Tri
##############################################################################
Anova(spad.tri)

## Canopy plot
spad_tri_canopy_results <- cld(emmeans(spad.tri, ~canopy, type = "response"),
                               Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spad_tri_canopy_plot <- ggplot(data = subset(df, spp == "Tri"),
                       aes(x = canopy, y = SPAD, fill = canopy)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = canopy.colors) +
  geom_text(data = spad_tri_canopy_results, 
            aes(x = canopy, y = 60, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(15, 60), breaks = seq(15, 60, 15)) +
  labs(x = "Measurement period",
       y = "SPAD (unitless)") +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
spad_tri_canopy_plot

## GM trt plot
spad_tri_gmtrt_results <- cld(emmeans(spad.tri, ~gm.trt, type = "response"),
                               Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spad_tri_gmtrt_plot <- ggplot(data = subset(df, spp == "Tri"),
                               aes(x = gm.trt, y = SPAD, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = spad_tri_gmtrt_results, 
            aes(x = gm.trt, y = 60, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(15, 60), breaks = seq(15, 60, 15)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "SPAD (unitless)") +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
spad_tri_gmtrt_plot

##############################################################################
## SPAD - Mai
##############################################################################
Anova(spad.mai)

## Canopy plot
spad_mai_canopy_results <- cld(emmeans(spad.mai, ~canopy, type = "response"),
                               Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spad_mai_canopy_plot <- ggplot(data = subset(df, spp == "Mai"),
                               aes(x = canopy, y = SPAD, fill = canopy)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = spad_mai_canopy_results, 
            aes(x = canopy, y = 60, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(15, 60), breaks = seq(15, 60, 15)) +
  labs(x = "Measurement period",
       y = "SPAD (unitless)") +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
spad_mai_canopy_plot

## GM trt plot
spad_mai_gmtrt_results <- cld(emmeans(spad.mai, ~gm.trt, type = "response"),
                              Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spad_mai_gmtrt_plot <- ggplot(data = subset(df, spp == "Mai"),
                              aes(x = gm.trt, y = SPAD, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = gm.colors) +
  geom_text(data = spad_mai_gmtrt_results, 
            aes(x = gm.trt, y = 60, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(15, 60), breaks = seq(15, 60, 15)) +
  labs(x = expression(bolditalic("A. petiolata")*bold(" treatment")),
       y = "SPAD (unitless)") +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
spad_mai_gmtrt_plot


##############################################################################
## iWUE - Tri
##############################################################################
Anova(iwue.tri)

## GM trt plot
iwue_tri_results <- cld(emmeans(iwue.tri, ~gm.trt*canopy, type = "response"),
                              Letters = LETTERS, reverse = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

iwue_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = canopy, y = iwue, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = iwue_tri_results, 
            aes(x = canopy, y = 150, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("i")*"WUE ("*mu*"mol mol"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
iwue_tri_plot

##############################################################################
## iWUE - Mai
##############################################################################
Anova(iwue.mai)

## GM trt plot
iwue_mai_results <- cld(emmeans(iwue.mai, ~gm.trt*canopy, type = "response"),
                        Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

iwue_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                        aes(x = canopy, y = iwue, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = iwue_mai_results, 
            aes(x = canopy, y = 150, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("i")*"WUE ("*mu*"mol mol"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
iwue_mai_plot

##############################################################################
## SPAD:gs - Tri
##############################################################################
Anova(spad.gs.tri)

spadgs_tri_results <- cld(emmeans(spad.gs.tri, ~gm.trt*canopy, type = "response"),
                        Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spadgs_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = canopy, y = spad.gs, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = spadgs_tri_results, 
            aes(x = canopy, y = 600, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, 150)) +
  labs(x = "Measurement period",
       y = expression(bold("SPAD:g"["sw"]*" (1/ mol m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
spadgs_tri_plot

##############################################################################
## SPAD:gs - Mai
##############################################################################
Anova(spad.gs.mai)

spadgs_mai_results <- cld(emmeans(spad.gs.mai, ~gm.trt*canopy, type = "response"),
                          Letters = LETTERS, reverse = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spadgs_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                          aes(x = canopy, y = spad.gs, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = iwue_mai_results, 
            aes(x = canopy, y = 4000, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 4000), breaks = seq(0, 4000, 1000)) +
  labs(x = "Measurement period",
       y = expression(bold("SPAD:g"["sw"]*" (1/ mol m"^"-2"*"s"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
spadgs_mai_plot

##############################################################################
## Vcmax:gs - Tri
##############################################################################
Anova(vcmax.gs.tri)

vcmaxgs_tri_results <- cld(emmeans(vcmax.gs.tri, ~gm.trt*canopy, type = "response"),
                          Letters = LETTERS, reverse = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

vcmaxgs_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                           aes(x = canopy, y = vcmax.gs, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = vcmaxgs_tri_results, 
            aes(x = canopy, y = 1600, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("V")["cmax25"]*":g"["sw"]*" ("*mu*"mol mol"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
vcmaxgs_tri_plot

##############################################################################
## Vcmax:gs - Mai
##############################################################################
Anova(vcmax.gs.mai)

vcmaxgs_mai_results <- cld(emmeans(vcmax.gs.mai, ~gm.trt*canopy, type = "response"),
                           Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

vcmaxgs_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                           aes(x = canopy, y = vcmax.gs, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = vcmaxgs_mai_results, 
            aes(x = canopy, y = 1600, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open canopy", "closed canopy")) +
  scale_y_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  labs(x = "Measurement period",
       y = expression(bold(italic("V")["cmax25"]*":g"["sw"]*" ("*mu*"mol mol"^"-1"*")")),
       fill = expression(bolditalic("A. petiolata")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
vcmaxgs_mai_plot

##############################################################################
## Figure 1: Soil nutrients 
##############################################################################
png("../drafts/figs/TT23_fig1_soilNutrients.png", width = 12, height = 8,
    units = "in", res = 600)
ggarrange(nitrogen_canopy_plot, phosphate_canopy_plot, soil_np_canopy_plot,
          nitrogen_gmtrt_plot, phosphate_gmtrt_plot, soil_np_gmtrt_plot,
          ncol = 3, nrow = 2, hjust = 0,
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 2: Gas exchange
##############################################################################
png("../drafts/figs/TT23_fig2_gasExchange.png", width = 12, height = 8,
    units = "in", res = 600)
ggarrange(anet_tri_plot, gsw_tri_plot, stomlim_tri_plot, anet_mai_plot,  
          gsw_mai_plot, stomlim_mai_plot, common.legend = TRUE, hjust = 0,
          legend = "right", ncol = 3, nrow = 2, align = "hv",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 3: Photosynthetic capacity
##############################################################################
png("../drafts/figs/TT23_fig3_photoCapacity.png", 
    width = 12, height = 8, units = "in", res = 600)
ggarrange(vcmax_tri_plot, jmax_tri_plot, jvmax_tri_plot, vcmax_mai_plot,  
          jmax_mai_plot, jvmax_mai_plot, common.legend = TRUE, 
          legend = "bottom", ncol = 3, nrow = 2, align = "hv", hjust = 0,
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 4: Resource use efficiency
##############################################################################
png("../drafts/figs/TT23_fig4_nitrogenWater_tradeoffs.png", 
    width = 12, height = 8, units = "in", res = 600)
ggarrange(iwue_tri_plot, spadgs_tri_plot, vcmaxgs_tri_plot, 
          iwue_mai_plot, spadgs_mai_plot, vcmaxgs_mai_plot,
          common.legend = TRUE, legend = "bottom", ncol = 3, 
          nrow = 2, hjust = 0, align = "hv", 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure SX: Chlorophyll fluorescence
##############################################################################
png("../drafts/figs/TT23_figSX_chlor_fluor.png", 
    width = 8, height = 8, units = "in", res = 600)
ggarrange(spad_tri_gmtrt_plot, spad_tri_canopy_plot,
          spad_mai_gmtrt_plot, spad_mai_canopy_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 2, 
          align = "hv", font.label = list(size = 18), hjust = 0,
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))
dev.off()

