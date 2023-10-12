## Load libraries
library(tidyverse)
library(ggpubr)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)
library(glmmTMB)

## Read in compiled curve fit file
df <- read.csv("../data/TT23_phys_data_long.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "invaded")),
         canopy = factor(canopy, levels = c("pre_closure", "post_closure")))
head(df)

##############################################################################
## Vcmax regressed against garlic mustard treatment and canopy openness
##############################################################################
df$vcmax25[225] <- NA

vcmax <- lmer(log(vcmax25) ~ gm.trt * spp * canopy + (1 | plot), 
              data = subset(df, spp != "Ari"))

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

# Post-hoc comparisons
cld(emmeans(vcmax, pairwise~spp*canopy, type = "response"))
cld(emmeans(vcmax, pairwise~canopy, type = "response"))


# Create compact letters for plot
vcmax.letters <- cld(emmeans(vcmax, pairwise~canopy*spp), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both")) %>% data.frame()

## Create plot
vcmax.cat <- ggplot(data = subset(df, spp != "Ari"), 
                    aes(x = canopy, y = vcmax25, fill = spp)) +
  geom_boxplot(alpha = 0.75) +
  geom_point(shape = 21, size = 2, alpha = 0.75, 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75)) +
  geom_text(data = vcmax.letters, aes(y = 200, label = .group),
            position = position_dodge(width = 0.75), size = 5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 200)) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  labs(x = "Canopy status", 
       y = expression(bold("V"["cmax25"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
vcmax.cat

##############################################################################
## Vcmax regressed against subplot GM densitiy and canopy openness
##############################################################################
df$vcmax25[c(20, 225)] <- NA

vcmax.cont <- lmer(log(vcmax25) ~ total_subplot * spp * canopy + (1 | plot), 
              data = subset(df, spp != "Ari"))

# Check model assumptions
plot(vcmax.cont)
qqnorm(residuals(vcmax.cont))
qqline(residuals(vcmax.cont))
densityPlot(residuals(vcmax.cont))
shapiro.test(residuals(vcmax.cont))
outlierTest(vcmax.cont)

# Model output
summary(vcmax.cont)
Anova(vcmax.cont)
r.squaredGLMM(vcmax.cont)

# Post-hoc comparisons
test(emtrends(vcmax.cont, ~spp*canopy, "total_subplot"))



# Create compact letters for plot
vcmax.trend <- emmeans(vcmax.cont, ~total_subplot*spp*canopy, 
                       at = list(total_subplot = seq(0, 145, 1)),
                       type = "response") %>% data.frame() %>%
  mutate(linetype = ifelse(spp == "Mai" & canopy == "pre_closure",
                           "solid", "dashed"))
vcmax.trend

## Create plot
vcmax.cat <- ggplot(data = subset(df, spp != "Ari"), 
                    aes(x = total_subplot, y = vcmax25, fill = spp)) +
  geom_point(aes(shape = canopy), size = 2, alpha = 0.75) +
  geom_smooth(data = vcmax.trend, 
              aes(y = response, linetype = linetype, color = spp)) +
  facet_grid(~canopy)
  geom_boxplot(alpha = 0.75) +
  geom_point(shape = 21, size = 2, alpha = 0.75, 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75)) +
  geom_text(data = vcmax.letters, aes(y = 200, label = .group),
            position = position_dodge(width = 0.75), size = 5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 200)) +
  scale_x_discrete(labels = c("Open", "Closed")) +
  labs(x = "Canopy status", 
       y = expression(bold("V"["cmax25"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
vcmax.cat





##############################################################################
## Jmax
##############################################################################
df$jmax25[c(20, 178)] <- NA

jmax <- lmer(log(jmax25) ~ total_subplot * spp * canopy + (1 | plot), 
             data = subset(df, spp != "Ari"))

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

# Post-hoc comparisons
test(emtrends(jmax, ~spp*canopy, "total_subplot"))


cld(emtrends(jmax, 
             
             
             
             pairwise~total_subplot*spp*canopy))

# Create compact letters for plot
jmax.letters <- cld(emmeans(jmax, pairwise~gm.trt*spp), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both")) %>% data.frame()

## Create plot
jmax.cat <- ggplot(data = df, 
                    aes(x = gm.trt, y = jmax, fill = spp)) +
  geom_boxplot(alpha = 0.75) +
  geom_point(shape = 21, size = 2, alpha = 0.75, 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75)) +
  geom_text(data = jmax.letters, aes(y = 300, label = .group),
            position = position_dodge(width = 0.75), size = 5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 75)) +
  labs(x = "Garlic mustard treatment", 
       y = expression(bold("J"["max25"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
jmax.cat

##############################################################################
## Write plot for Vcmax and Jmax
##############################################################################
png("../figs/TT23_preclosure_photo.png", width = 12, height = 4.5, units = "in",
    res = 600)
ggarrange(vcmax.cat, jmax.cat, ncol = 2, common.legend = TRUE, legend = "right")
dev.off()

##############################################################################
## Jmax:Vcmax
##############################################################################
df$jmax.vcmax[91] <- NA

jvmax <- lmer(log(jmax.vcmax) ~ gm.trt * spp + (1 | plot), data = df)

# Check model assumptions
plot(jvmax)
qqnorm(residuals(jvmax))
qqline(residuals(jvmax))
densityPlot(residuals(jvmax))
shapiro.test(residuals(jvmax))
outlierTest(jvmax)

# Model output
summary(jvmax)
Anova(jvmax)
r.squaredGLMM(jvmax)

# Post-hoc comparisons
emmeans(jvmax, pairwise~gm.trt, type = "response") ## Lower ratio in weeded plots


##############################################################################
## Ci:Ca
##############################################################################
df$ci.ca[c(47, 83)] <- NA

ci.ca <- lmer(ci.ca ~ gm.trt * spp + (1 | plot), data = df)

# Check model assumptions
plot(ci.ca)
qqnorm(residuals(ci.ca))
qqline(residuals(ci.ca))
densityPlot(residuals(ci.ca))
shapiro.test(residuals(ci.ca))
outlierTest(ci.ca)

# Model output
summary(ci.ca)
Anova(ci.ca)
r.squaredGLMM(ci.ca)

# Post-hoc comparisons
cld(emmeans(ci.ca, pairwise~spp*gm.trt), Letters = LETTERS)
## GM invasion insignificantly increases Ci:Ca in Trillium; decreases
## Ci:Ca in Maianthemum

##############################################################################
## gsw
##############################################################################
df$gsw[c(47,91,111,117,128)] <- NA

gsw <- lmer(log(gsw) ~ gm.trt * spp + (1 | plot), data = df)

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

# Post-hoc comparisons
cld(emmeans(gsw, pairwise~spp*gm.trt, type = "response"), Letters = LETTERS)
## GM invasion decreases gsw, but only in Maianthemum

# Create compact letters for plot
gsw.letters <- cld(emmeans(gsw, pairwise~gm.trt*spp), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both")) %>% data.frame()

## Create plot
gsw.cat <- ggplot(data = df, 
                    aes(x = gm.trt, y = gsw, fill = spp)) +
  geom_boxplot(alpha = 0.75) +
  geom_point(shape = 21, size = 2, alpha = 0.75, 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75)) +
  geom_text(data = gsw.letters, aes(y = 0.3, label = .group),
            position = position_dodge(width = 0.75), size = 5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(x = "Garlic mustard treatment", 
       y = expression(bold("g"["sw"]*" (mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
gsw.cat

##############################################################################
## iwue
##############################################################################
df$iwue[c(39, 78, 83, 94, 101)] <- NA

iwue <- lmer(log(iwue) ~ gm.trt * spp + (1 | plot), data = df)

# Check model assumptions
plot(iwue)
qqnorm(residuals(iwue))
qqline(residuals(iwue))
densityPlot(residuals(iwue))
shapiro.test(residuals(iwue))
outlierTest(iwue)

# Model output
summary(iwue)
Anova(iwue)
r.squaredGLMM(iwue)

# Post-hoc comparisons
cld(emmeans(iwue, pairwise~spp*gm.trt, type = "response"), Letters = LETTERS)
## GM invasion increases water use efficiency, but only in Maianthemum

# Create compact letters for plot
iwue.letters <- cld(emmeans(iwue, pairwise~gm.trt*spp), Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both")) %>% data.frame()

## Create plot
iwue.cat <- ggplot(data = df, 
                  aes(x = gm.trt, y = iwue, fill = spp)) +
  geom_boxplot(alpha = 0.75) +
  geom_point(shape = 21, size = 2, alpha = 0.75, 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             dodge.width = 0.75)) +
  geom_text(data = iwue.letters, aes(y = 150, label = .group),
            position = position_dodge(width = 0.75), size = 5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 150)) +
  labs(x = "Garlic mustard treatment", 
       y = expression(bold("iWUE ("*mu*"mol"*" mol"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
iwue.cat

##############################################################################
## Write plot for gsw and iWUE
##############################################################################
png("../figs/TT23_preclosure_waterUse.png", width = 12, height = 4.5, units = "in",
    res = 600)
ggarrange(gsw.cat, iwue.cat, ncol = 2, common.legend = TRUE, legend = "right")
dev.off()

##############################################################################
##############################################################################
## Summary of categorical treatments
##############################################################################
##############################################################################
# Garlic mustard treatments have limited effects on leaf physiology, namely
# Rubisco carboxylation and RuBP regeneration rates. However, treatments seem
# to have a strong impact on water usage, which are particularly apparent in 
# Maianthemum. GM presence decreased Maianthemum stomatal conductance and water
# use efficiency.

# However, I am curious if the lack of garlic mustard presence in plot 7 might
# be skewing these results. To address this, the following code substitutes
# garlic mustard density for garlic mustard treatment. Note that this approach
# creates large zero inflation in models, so a zero-inflated, generalized linear
# mixed effect model is used.


##############################################################################
## Zero-inflated GLMM for Vcmax
##############################################################################
df$vcmax[c(20, 34)] <- NA

vcmax.gmdens <- lmer(sqrt(vcmax) ~ gm.total.dens * spp + (1 | plot), data = df)

# Check model assumptions
plot(vcmax.gmdens)
qqnorm(residuals(vcmax.gmdens))
qqline(residuals(vcmax.gmdens))
densityPlot(residuals(vcmax.gmdens))
shapiro.test(residuals(vcmax.gmdens))
outlierTest(vcmax.gmdens)

# Model output
summary(vcmax.gmdens)
Anova(vcmax.gmdens)
r.squaredGLMM(vcmax.gmdens)

# Post-hoc comparisons
emmeans(vcmax.gmdens, pairwise~spp, type = "response")
test(emtrends(vcmax.gmdens, ~spp, "gm.total.dens", type = "response"))

# Model predictions across range in x-axis values for plotting
vcmax.gmdens.trend <- emmeans(vcmax.gmdens, ~spp, "gm.total.dens", 
                              type = "response", 
                              at = list(gm.total.dens = seq(0,164, 1))) %>%
  data.frame() %>% mutate(linetype = ifelse(spp == "Mai", "solid", "dashed"))

# Plot
vcmax.gmdens.plot <- ggplot(data = df, aes(x = gm.total.dens, y = vcmax, fill = spp)) +
  geom_jitter(shape = 21, size = 3, alpha = 0.75) +
  geom_ribbon(data = vcmax.gmdens.trend, aes(y = response, ymin = lower.CL,
                                             ymax = upper.CL), alpha = 0.25) +
  geom_smooth(data = vcmax.gmdens.trend, aes(y = response, color = spp, 
                                             linetype = linetype), size = 2) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Garlic mustard subplot density", 
       y = expression(bold("V"["cmax25"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
vcmax.gmdens.plot


##############################################################################
## Zero-inflated GLMM for Jmax
##############################################################################
df$jmax[18] <- NA

jmax.gmdens <- lmer(log(jmax) ~ gm.total.dens * spp + (1 | plot), data = df)

# Check model assumptions
plot(jmax.gmdens)
qqnorm(residuals(jmax.gmdens))
qqline(residuals(jmax.gmdens))
densityPlot(residuals(jmax.gmdens))
shapiro.test(residuals(jmax.gmdens))
outlierTest(jmax.gmdens)

# Model output
summary(jmax.gmdens)
Anova(jmax.gmdens)
r.squaredGLMM(jmax.gmdens)

# Post-hoc comparisons
test(emtrends(jmax.gmdens, ~spp, "gm.total.dens"))

# Model predictions across range in x-axis values for plotting
jmax.gmdens.trend <- emmeans(jmax.gmdens, ~spp, "gm.total.dens", 
                              type = "response", 
                              at = list(gm.total.dens = seq(0, 164, 1))) %>%
  data.frame() %>% mutate(linetype = ifelse(spp == "Mai", "solid", "dashed"))

# Plot
jmax.gmdens.plot <- ggplot(data = df, aes(x = gm.total.dens, y = jmax, fill = spp)) +
  geom_jitter(shape = 21, size = 3, alpha = 0.75) +
  geom_ribbon(data = jmax.gmdens.trend, aes(y = response, ymin = lower.CL,
                                             ymax = upper.CL), alpha = 0.25) +
  geom_smooth(data = jmax.gmdens.trend, 
              aes(y = response, color = spp, linetype = linetype), size = 2) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Garlic mustard subplot density", 
       y = expression(bold("J"["max25"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
jmax.gmdens.plot

##############################################################################
## Write plot for Vcmax and Jmax
##############################################################################
png("../figs/TT23_preclosure_photo_continuous.png", 
    width = 12, height = 4.5, units = "in", res = 600)
ggarrange(vcmax.gmdens.plot, jmax.gmdens.plot, ncol = 2, common.legend = TRUE, 
          legend = "right", align = "hv")
dev.off()

##############################################################################
## Zero-inflated GLMM for gsw
##############################################################################
gsw.gmdens <- lmer(sqrt(gsw) ~ gm.total.dens * spp + (1 | plot),  data = df)

# Check model assumptions
plot(gsw.gmdens)
qqnorm(residuals(gsw.gmdens))
qqline(residuals(gsw.gmdens))
densityPlot(residuals(gsw.gmdens))
shapiro.test(residuals(gsw.gmdens))
outlierTest(gsw.gmdens)

# Model output
summary(gsw.gmdens)
Anova(gsw.gmdens)
r.squaredGLMM(gsw.gmdens)

# Post-hoc comparisons
test(emtrends(gsw.gmdens, ~spp, "gm.total.dens"))

# Model predictions across range in x-axis values for plotting
gsw.gmdens.trend <- emmeans(gsw.gmdens, ~spp, "gm.total.dens", 
                             type = "response", 
                             at = list(gm.total.dens = seq(0, 164, 1))) %>%
  data.frame()

# Plot
gsw.gmdens.plot <- ggplot(data = df, aes(x = gm.total.dens, y = gsw, fill = spp)) +
  geom_jitter(shape = 21, size = 3, alpha = 0.75) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Garlic mustard subplot density", 
       y = expression(bold("g"["sw"]*" (mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
gsw.gmdens.plot

##############################################################################
## Zero-inflated GLMM for iWUE
##############################################################################
iwue.gmdens <- lmer(log(iwue) ~ gm.total.dens * spp + (1 | plot), data = df)

# Check model assumptions
plot(iwue.gmdens)
qqnorm(residuals(iwue.gmdens))
qqline(residuals(iwue.gmdens))
densityPlot(residuals(iwue.gmdens))
shapiro.test(residuals(iwue.gmdens))
outlierTest(iwue.gmdens)

# Model output
summary(iwue.gmdens)
Anova(iwue.gmdens)
r.squaredGLMM(iwue.gmdens)

# Post-hoc comparisons
test(emtrends(iwue.gmdens, ~spp, "gm.total.dens", type = "response"))

# Model predictions across range in x-axis values for plotting
iwue.gmdens.trend <- emmeans(iwue.gmdens, ~spp, "gm.total.dens", 
                            type = "response", 
                            at = list(gm.total.dens = seq(0, 164, 1))) %>%
  data.frame() %>% mutate(linetype = ifelse(spp == "Mai", "solid", "dashed"))

# Plot
iwue.gmdens.plot <- ggplot(data = df, aes(x = gm.total.dens, y = iwue, fill = spp)) +
  geom_jitter(shape = 21, size = 3, alpha = 0.75) +
  geom_ribbon(data = iwue.gmdens.trend, aes(y = response, ymin = lower.CL,
                                            ymax = upper.CL), alpha = 0.25) +
  geom_smooth(data = iwue.gmdens.trend, 
              aes(y = response, color = spp, linetype = linetype), size = 2) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 150)) +
  labs(x = "Garlic mustard subplot density", 
       y = expression(bold("iWUE ("*mu*"mol"*" mol"^"-1"*")")),
       fill = "Species", color = "Species") +
  guides(linetype = "none") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
iwue.gmdens.plot

##############################################################################
## Write plot for Vcmax and Jmax
##############################################################################
png("../figs/TT23_preclosure_waterUsage_continuous.png", 
    width = 12, height = 4.5, units = "in", res = 600)
ggarrange(gsw.gmdens.plot, iwue.gmdens.plot, ncol = 2, common.legend = TRUE, 
          legend = "right", align = "hv")
dev.off()




##############################################################################
## Stomatal limitation regressed against GM total density
##############################################################################
stom.lim.gmdens <- lmer(log(stom.lim) ~ gm.total.dens * spp + (1 | plot), data = df)

# Check model assumptions
plot(stom.lim.gmdens)
qqnorm(residuals(stom.lim.gmdens))
qqline(residuals(stom.lim.gmdens))
densityPlot(residuals(stom.lim.gmdens))
shapiro.test(residuals(stom.lim.gmdens))
outlierTest(stom.lim.gmdens)

# Model output
summary(stom.lim.gmdens)
Anova(stom.lim.gmdens)
r.squaredGLMM(stom.lim.gmdens)

# Post-hoc comparisons
emmeans(stom.lim.gmdens, pairwise~spp)
test(emtrends(stom.lim.gmdens, ~spp, "gm.total.dens"))



## Vcmax
ggplot(data = subset(df, gm.trt == "invaded"), aes(x = gm.total.dens, y = vcmax, fill = spp)) +
  geom_point(shape = 21, size = 4) +
  geom_smooth(method = 'lm') +
  labs(x = "Garlic mustard density (rosette + adults)", 
       y = expression(bold("V"["cmax"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))




## Jmax
ggplot(data = df, aes(x = spp, y = jmax, fill = gm.trt)) +
  geom_jitter(shape = 21, position = position_jitterdodge(jitter.width = 0.1, 
                                                          dodge.width = 0.75)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Garlic mustard treatment", 
       y = expression(bold("J"["max"]*" ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")),
       fill = "Species", color = "Species") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

## Jmax:Vcmax
ggplot(data = df, aes(x = gm.trt, y = jmax.vcmax)) +
  geom_boxplot(aes(fill = spp))

## Ci:Ca
ggplot(data = df, aes(x = gm.trt, y = ci.ca)) +
  geom_boxplot(aes(fill = spp))

## iWUE
ggplot(data = df, aes(x = gm.trt, y = iwue)) +
  geom_boxplot(aes(fill = spp))

## Rd:Vcmax
ggplot(data = df, aes(x = gm.trt, y = rd.vcmax)) +
  geom_boxplot(aes(fill = spp))