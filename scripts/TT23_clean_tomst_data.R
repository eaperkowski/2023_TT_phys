## Libraries
library(tidyverse)
library(lubridate)

## Data frame containing sensor IDs and their plot # / treatment
sensor_meta <- data.frame(sensor_id = seq(94232156, 94232165, 1),
                          plot = c(7, 7, 3, 3, 4, 
                                   6, 4, 5, 5, 6),
                          position = c("NW", "W", "W", "NW", "NW",
                                       "W", "W", "W", "NW", "NW"))
sensor_meta$gm.trt <- ifelse(sensor_meta$position == "NW",
                             "ambient", "weeded")

## Read in files
files <- list.files("../data/tomst_probe_calibrated/", 
                    recursive = T,
                    pattern = "\\.csv$", 
                    full.names = T)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))
tomst <- lapply(files, read.csv)
tomst <- purrr::map(names(tomst), ~tomst[[.x]] %>%
                      mutate(sensor_id = .x))


## Add sensor ID as column to each data frame
for(i in seq_along(tomst)) {
  tomst[[i]]$sensor_id <- str_extract(string = tomst[[i]]$sensor_id, 
                                      pattern = "[0-9]{8}")
}

## Merge files into single data.frame, filter to appropriate
## temporal scale
tomst_df <- reshape::merge_all(tomst) %>%
  dplyr::select(sensor_id, date_time = Date...time, 
         soil_temp = Soil.temperature..6.cm,
         surf_temp = Surface.temperature,
         air_temp = Air.temperature..12.cm,
         vwc = Vol..moisture)

tomst_fully_cleaned <- tomst_df %>%
  mutate(sensor_id = as.numeric(sensor_id),
         date_time = with_tz(ymd_hm(date_time, tz = "UTC"), "US/Eastern")) %>%
  filter(date_time > "2023-04-26" & date_time < "2023-07-01") %>%
  full_join(sensor_meta) %>%
  dplyr::select(sensor_id, date_time, plot, gm.trt, soil_temp:vwc)

## Daily soil moisture summary by plot, trt
daily_sm <- tomst_fully_cleaned %>%
  mutate(day = date(date_time)) %>%
  group_by(day, sensor_id, plot, gm.trt) %>%
  summarize(daily_vwc = mean(vwc)) %>%
  mutate(doy = yday(day))

## Daily soil moisture summary by plot, trt only including
## plots from 2023 field season (plots 3, 5, 7)
daily_sm_tt23 <- tomst_fully_cleaned %>%
  filter(plot != 4 & plot != 6) %>%
  mutate(day = date(date_time)) %>%
  group_by(day, sensor_id, plot, gm.trt) %>%
  summarize(daily_vwc = mean(vwc)) %>%
  mutate(doy = yday(day))

# write.csv(daily_sm_tt23, "../data/TT23_tomst_probe_sm_daily.csv", 
#           row.names = F)
