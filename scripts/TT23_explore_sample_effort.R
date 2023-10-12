library(ggplot2)
library(dplyr)
library(ggpubr)

df.open <- read.csv("../data/TT23_pre_closure_phys_data.csv")
df.closed <- read.csv("../log/TT_tag_metadata_canopy_closed.csv")

open.summary <- df.open %>%
  count(plot, gm.trt, spp) %>%
  select(plot:spp, n_open = n)


closed.summary <- df.closed %>%
  count(plot, gm.trt, spp) %>%
  select(plot:spp, n_closed = n)

summary.total <- open.summary %>%
  full_join(closed.summary)

## Open canopy summary plot
open.plot <- ggplot(data = summary.total) +
  geom_bar(aes(x = spp, y = n_open, fill = gm.trt),
           stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Individuals w/ A/Ci curves",
       fill = "GM treatment", title = "Open canopy") +
  theme_bw(base_size = 18)

## Closed canopy summary plot
closed.plot <- ggplot(data = summary.total) +
  geom_bar(aes(x = spp, y = n_closed, fill = gm.trt),
           stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Species", y = "Individuals w/ A/Ci curves",
       fill = "GM treatment", title = "Closed canopy") +
  theme_bw(base_size = 18)

png("../figs/TT23_sample_effort_plots.png", width = 12, height = 5,
    units = "in", res = 600)
ggarrange(open.plot, closed.plot, nrow = 1, ncol = 2, common.legend = TRUE,
          legend = "right", align = "hv")
dev.off()


