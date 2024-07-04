## Libraries
library(tidyverse)
library(lubridate)

## Data frame containing sensor IDs and their plot # / treatment
sensor_meta <- data.frame(sensor_id = seq(94232156, 94232165, 1),
                          plot = c(7, 7, 3 , 3, 4, 
                                   6, 4, 5, 5, 6),
                          position = c("NW", "W", "W", "NW", "NW",
                                       "W", "W", "W", "NW", "NW"))
sensor_meta$gm.trt <- ifelse(sensor_meta$position == "NW",
                             "ambient", "weeded")

## Read in files
files <- list.files("../data/tomst_probe_raw/", 
                    recursive = T,
                    pattern = "\\.csv$", 
                    full.names = T)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))
tomst <- lapply(files, read.csv)
tomst <- purrr::map(names(tomst), ~tomst[[.x]] %>%
                      mutate(sensor_id = .x))

## Add sensor ID 
for(i in seq_along(tomst)) {
  tomst[[i]]$sensor_id <- str_extract(string = tomst[[i]]$sensor_id, 
                                      pattern = "[0-9]{8}")
}

## Merge files into single data.frame, filter to appropriate
## temporal scale
tomst_df <- reshape::merge_all(tomst) %>%
  dplyr::select(sensor_id, date:errFlag)

tomst_df_clean <- tomst_df %>%
  mutate(sensor_id = as.numeric(sensor_id),
         date = with_tz(ymd_hm(date, tz = "UTC"), "US/Eastern"),
         quad_a = 0.000000017, # quad params for silt loam
         quad_b = 0.000118119, # quad params for silt loam
         quad_c = -0.101168511) %>% # quad params for silt loam
  group_by(sensor_id) %>%
  filter(date > "2023-04-26" & date < "2023-07-01") %>%
  full_join(sensor_meta) %>%
  dplyr::select(sensor_id, date, plot, position, gm.trt, 
                t1:soil_moisture_raw, quad_a:quad_c)

## Calculate volumetric soil moisture
sm_volumetric <- tomst_df_clean %>%
  mutate(soil_moisture_volum = (quad_a * (soil_moisture_raw)^2) +
           (quad_b * soil_moisture_raw) + quad_c)

## Daily soil moisture summary by plot, trt
daily_sm <- sm_volumetric %>%
  mutate(day = date(date)) %>%
  group_by(day, sensor_id, plot, position, gm.trt) %>%
  summarize(daily_sm = mean(soil_moisture_volum)) %>%
  mutate(doy = yday(day))

## Daily soil moisture summary by plot, trt only including
## plots from 2023 field season (plots 3, 5, 7)
daily_sm_tt23 <- sm_volumetric %>%
  filter(plot != 4 & plot != 6) %>%
  mutate(day = date(date)) %>%
  group_by(day, sensor_id, plot, position, gm.trt) %>%
  summarize(daily_sm = mean(soil_moisture_volum)) %>%
  mutate(doy = yday(day))

# write.csv(daily_sm_tt23, "../data/TT23_tomst_probe_sm_daily.csv", row.names = F)