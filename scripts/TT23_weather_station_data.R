## Libraries
library(tidyverse)
library(lubridate)
library(bigleaf)

## Load weather data
wx <- read.csv("../data/TT_weather_data_05222023.csv")

## Convert solar radiation from w m^-2 to PPFD
wx$ppfd <- Rg.to.PPFD(Rg = wx$solar.rad.wm2)
wx$date <- ymd_hms(wx$date, tz = "US/Eastern")

# Visualize data 
head(wx)

# Calculate mean daytime growing conditions
wx %>%
  mutate(month = month(date),
         day = day(date)) %>%
  unite("month.day", month:day) %>%
  group_by(month.day) %>%
  summarize(ppfd.max = max(ppfd, na.rm = TRUE),
            temp.max = max(air.temp.c, na.rm = TRUE),
            temp.min = min(air.temp.c, na.rm = TRUE),
            vpd.max = max(vpd.kpa, na.rm = TRUE),
            vpd.min = min(vpd.kpa, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(ppfd = mean(ppfd.max),
            temp.max = mean(temp.max, na.rm =TRUE),
            temp.min = mean(temp.min, na.rm = TRUE),
            vpd.max = mean(vpd.max, na.rm = TRUE),
            vpd.min = mean(vpd.min, na.rm = TRUE))


precip <- wx %>% mutate(month = month(date),
                         day = day(date),
                         year = year(date)) %>%
  unite("date", month:year, sep = "/") %>%
  mutate(date = mdy(date)) %>%
  group_by(date) %>%
  summarize(precip.total = sum(precip.mm))

# Plot PPFD
ppfd.plot <- ggplot(data = wx, aes(x = date, y = ppfd)) +
  geom_hline(yintercept = 349, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) +
  geom_line(size = 2) +
  scale_y_continuous(limits = c(0,800)) +
  labs(x = "Date", 
       y = expression("PPFD ("*mu*"mol"*" m"^"-2"*" s"^"-1"*")")) +
  theme_bw(base_size = 18)
ppfd.plot

# Plot temperature
airtemp.plot <- ggplot(data = wx, aes(x = date, y = air.temp.c)) +
  geom_hline(yintercept = 17.1, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 6.99, color = "blue", linetype = "dashed", size = 1) +
  geom_line(size = 2) +
  scale_y_continuous(limits = c(-1.5, 30)) +
  labs(x = "Date", 
       y = expression("Temperature ("*degree*"C)")) +
  theme_bw(base_size = 18)
airtemp.plot

# Plot VPD
vpd.plot <- ggplot(data = wx, aes(x = date, y = vpd.kpa)) +
  geom_hline(yintercept = 1.09, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0.07, color = "blue", linetype = "dashed", size = 1) +
  geom_line(size = 2) +
  scale_y_continuous(limits = c(0, 3)) +
  labs(x = "Date", 
       y = "VPD (kPa)") +
  theme_bw(base_size = 18)
vpd.plot

# Plot precipitation
prcp.plot <- ggplot(data = precip, aes(x = date, y = precip.total)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 15)) +
  labs(x = "Date", 
       y = "Daily precipitation (mm)") +
  theme_bw(base_size = 18)
prcp.plot

png("../figs/TT23_wx_as_of_05222023.png",
    width = 12, height = 8, units = "in", res = 600)
ggarrange(ppfd.plot, airtemp.plot, vpd.plot, prcp.plot, nrow = 2, ncol = 2,
          align = "hv")
dev.off()


