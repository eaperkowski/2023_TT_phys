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
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "invaded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")),
         spp = factor(spp, levels = c("Tri", "Mai", "Ari")),
         vcmax.gs = vcmax25 / gsw) %>%
  unite(col = "gm.canopy", gm.trt, canopy, sep = "_", remove = FALSE) %>%
  mutate(gm.canopy = factor(gm.canopy, levels = c("weeded_pre_closure",
                                                  "invaded_pre_closure",
                                                  "weeded_post_closure",
                                                  "invaded_post_closure")))
head(df)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, total_subplot, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:n_plantAvail_day),
               .funs = mean) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "invaded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")))

## Remove outliers
df$anet[c(35, 71, 80, 86, 116, 120)] <- NA
df$gsw[c(86, 120)] <- NA
df$stom.lim[c(228, 229)] <- NA
df$vcmax25[c(183, 231)] <- NA
df$jmax25[183] <- NA
df$jmax.vcmax[c(94, 184, 224, 225, 231)] <- NA
df$iwue[c(80, 86, 104)] <- NA

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
iwue.tri <- lmer(log(iwue) ~ gm.trt * canopy + (1 | plot),
                 data = subset(df, spp == "Tri"))
iwue.mai <- lmer(log(iwue) ~ gm.trt * canopy  + (1 | plot),
                 data = subset(df, spp == "Mai"))
vcmax.gs.tri <- lmer(log(vcmax.gs) ~ gm.trt * canopy + (1 | plot),
                     data = subset(df, spp == "Tri"))
vcmax.gs.mai <- lmer(log(vcmax.gs) ~ gm.trt * canopy + (1 | plot),
                     data = subset(df, spp == "Mai"))

## Add code for facet labels
facet.labs <- c("Trillium", "Maianthemum")
names(facet.labs) <- c("Tri", "Mai")

##############################################################################
## Soil N availability 
##############################################################################
# Prep file for figure
n_results <- cld(emmeans(plant_availableN, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = c("B", "B", "A", "A"))

n_plantavailable_plot <- ggplot(data = df.soil,
                                aes(x = canopy, y = n_plantAvail_day, 
                                    fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = n_results, 
            aes(x = canopy, y = 1, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) + 
  #annotate(geom = "text", x = 2.5, y = 0.875, size = 2.5, hjust = 1,
  #         label = expression(italic("Alliaria")*" presence: "*italic("p")*">0.05")) +
  #annotate(geom = "text", x = 2.5, y = 0.825, size = 2.5, hjust = 1,
  #         label = expression(bold("Canopy status: ")*bolditalic("p")*bold("<0.001"))) +
  #annotate(geom = "text", x = 2.5, y = 0.775, size = 2.5, hjust = 1,
  #         label = expression(italic("Alliaria")*"*Canopy: "*italic("p")*">0.05")) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Tree canopy status",
       y = expression(bold("Soil inorg. N (ppm day"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
n_plantavailable_plot

##############################################################################
## Phosphate figure  
##############################################################################
# Prep file for figure
phosphate_results <- cld(emmeans(phosphate, pairwise~canopy*gm.trt), 
                         Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both"))

phosphate_plot <- ggplot(data = df.soil,
                         aes(x = canopy, y = phosphate_ppm_day, 
                             fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = phosphate_results, 
            aes(x = canopy, y = 0.06, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 0.06), breaks = seq(0, 0.06, 0.015)) +
  labs(x = "Tree canopy status",
       y = expression(bold("Soil P (ppm day"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
phosphate_plot

##############################################################################
## Soil nitrate
##############################################################################
# Prep file for figure
no3_results <- cld(emmeans(nitrate, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = c("B", "B", "A", "A"))

no3_plot <- ggplot(data = df.soil,
                                aes(x = canopy, y = nitrate_ppm_day, 
                                    fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = no3_results, 
            aes(x = canopy, y = 1, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Tree canopy status",
       y = expression(bold("Soil NO"["3"]*"-N (ppm day"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
no3_plot

##############################################################################
## Soil ammonium
##############################################################################
# Prep file for figure
nh4_results <- cld(emmeans(ammonium, ~canopy*gm.trt), Letters = LETTERS) %>%
  mutate(.group = c("B", "B", "A", "A"))

nh4_plot <- ggplot(data = df.soil,
                   aes(x = canopy, y = ammonium_ppm_day, 
                       fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = nh4_results, 
            aes(x = canopy, y = 0.02, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.005)) +
  labs(x = "Tree canopy status",
       y = expression(bold("Soil NH"["4"]*"-N (ppm day"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
nh4_plot

##############################################################################
## Net photosynthesis - Trillium
##############################################################################
Anova(anet.tri)

anet_tri_results <- cld(emmeans(anet.tri, ~gm.trt, type = "response"), 
                Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = c("A", "A"))

anet_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.125) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  annotate(geom = "text", x = 1.5, y = 18, size = 4,
           label = expression(italic("Alliaria")*" treatment: "*italic("p")*">0.05")) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
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
                        Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                        aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 18, size = 4,
           label = expression(bolditalic("Alliaria")*bold(" treatment: ")*bolditalic("p")*bold("<0.001"))) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
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

gsw_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                  aes(x = gm.trt, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 0.3, size = 4,
           label = expression(italic("Alliaria")*" treatment: "*italic("p")*">0.05")) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
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

gsw_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                       aes(x = gm.trt, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 0.3, size = 4,
           label = expression(bolditalic("Alliaria")*bold(" treatment: ")*bolditalic("p")*bold("<0.001"))) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
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

stomlim_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                           aes(x = gm.trt, y = stom.lim, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 1, size = 4,
           label = expression(italic("Alliaria")*" treatment: "*italic("p")*">0.05")) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
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

stomlim_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                           aes(x = gm.trt, y = stom.lim, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 1, size = 4,
           label = expression(bolditalic("Alliaria")*bold(" treatment: ")*bolditalic("p")*bold("<0.05"))) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
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

vcmax_tri_results <- cld(emmeans(vcmax.tri, ~gm.trt*canopy, type = "response"), 
                           Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = c("C", "B", "A", "A"), spp = "Tri")

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
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
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
                         Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = c("B", "B", "A", "A"), spp = "Mai")

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
            aes(x = canopy, y = 200, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
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
                         Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = c("C", "B", "A", "A"))

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
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 100)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
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
                         Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = c("B", "B", "A", "A"))

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
            aes(x = canopy, y = 300, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 100)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
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
                        Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = c("B", "AB", "AB", "A"))

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
            aes(x = canopy, y = 0.8, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0.4, 0.8), breaks = seq(0.4, 0.8, 0.1)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"]*" (unitless)")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
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
            aes(x = canopy, y = 0.8, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D"),
                    labels = c(expression(italic("Alliaria")*" weeded"), 
                               expression(italic("Alliaria")*" present"))) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  scale_y_continuous(limits = c(0.4, 0.8), breaks = seq(0.4, 0.8, 0.1)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"]*" (unitless)")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
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
## iWUE - Tri
##############################################################################
Anova(iwue.tri)

iwue_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = gm.trt, y = iwue, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 200, size = 4,
           label = expression(italic("Alliaria")*italic(" treatment: ")*italic("p<0.1"))) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold("Intrinsic WUE ("*mu*"mol mol"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
iwue_tri_plot

##############################################################################
## iWUE - Mai
##############################################################################
Anova(iwue.mai)

iwue_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
         aes(x = gm.trt, y = iwue, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 200, size = 4,
           label = expression(bolditalic("Alliaria")*bold(" treatment: ")*bolditalic("p")*bold("<0.001"))) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold("Intrinsic WUE ("*mu*"mol mol"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
iwue_mai_plot

##############################################################################
## Vcmax:gs - Tri
##############################################################################
Anova(vcmax.gs.tri)

vcmaxgs_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                           aes(x = gm.trt, y = vcmax.gs, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 2000, size = 4,
           label = expression(italic("Alliaria")*" treatment: "*italic("p")*">0.05")) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 2000), breaks = seq(0, 2000, 500)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold(italic("V")["cmax25"]*" : g"["sw"]*" ("*mu*"mol mol"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
vcmaxgs_tri_plot

##############################################################################
## Vcmax:gs - Mai
##############################################################################
Anova(vcmax.gs.mai)

vcmaxgs_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                           aes(x = gm.trt, y = vcmax.gs, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  scale_fill_manual(values = c("#7BAFDE", "#F1932D")) +
  annotate(geom = "text", x = 1.5, y = 2000, size = 4,
           label = expression(bolditalic("Alliaria")*bold(" treatment: ")*bolditalic("p")*bold("<0.001"))) +
  scale_x_discrete(labels = c("weeded", "ambient")) +
  scale_y_continuous(limits = c(0, 2000), breaks = seq(0, 2000, 500)) +
  labs(x = expression(bolditalic("Alliaria")*bold(" treatment")),
       y = expression(bold(italic("V")["cmax25"]*" : g"["sw"]*" ("*mu*"mol mol"^"-1"*")"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
vcmaxgs_mai_plot

##############################################################################
## Figure 1: Soil nutrients 
##############################################################################
png("../drafts/figs/TT23_fig1_soilNutrients.png", width = 10, height = 8,
    units = "in", res = 600)
ggarrange(n_plantavailable_plot, phosphate_plot, no3_plot, nh4_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 2,
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 2: Gas exchange
##############################################################################
png("../drafts/figs/TT23_fig2_gasExchange.png", width = 12, height = 8,
    units = "in", res = 600)
ggarrange(anet_tri_plot, gsw_tri_plot, stomlim_tri_plot, anet_mai_plot, 
          gsw_mai_plot, stomlim_mai_plot, common.legend = TRUE, 
          legend = "right", ncol = 3, nrow = 2, align = "hv",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 3: Photosynthetic capacity
##############################################################################
png("../drafts/figs/TT23_fig3_photoCapacity.png", 
    width = 14, height = 8, units = "in", res = 600)
ggarrange(vcmax_tri_plot, jmax_tri_plot, jvmax_tri_plot, vcmax_mai_plot,  
          jmax_mai_plot, jvmax_mai_plot, common.legend = TRUE, 
          legend = "right", ncol = 3, nrow = 2, align = "hv",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 4: Resource use efficiency
##############################################################################
png("../drafts/figs/TT23_fig4_nitrogenWater_tradeoffs.png", 
    width = 8, height = 8, units = "in", res = 600)
ggarrange(iwue_tri_plot, vcmaxgs_tri_plot, iwue_mai_plot, vcmaxgs_mai_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 2, 
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)"), 
          font.label = list(size = 18))
dev.off()

