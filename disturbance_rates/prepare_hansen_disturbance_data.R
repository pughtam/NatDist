# Preparing Hansen disturbance data for non-spatial analysis, i.e., aggregating disturbance/forest pixels per year/landscape
#
# C. Senf
# 09.09.2019

# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)

# Get annual disturbance data ---------------------------------------------

# Studysite names

load("data/studysites_names_boreal.RData")
studysites_names_boreal <- studysites_names
load("data/studysites_names_temperate.RData")
studysites_names_temperate <- studysites_names
rm(studysites_names)

# Disturbance data

load("data/disturbances_boreal.RData")
disturbances_boreal <- disturbances
load("data/disturbances_temperate.RData")
disturbances_temperate <- disturbances
rm(disturbances)

# Forestcover data

load("data/forestcover_10perc_boreal.RData")
forestcover_boreal_10 <- forestcover_10
forestcover_boreal_10 <- forestcover_boreal_10[,-3]
names(forestcover_boreal_10) <- c("site", "forestcover_10")
load("data/forestcover_10perc_temperate.RData")
forestcover_temperate_10 <- forestcover_10
forestcover_temperate_10 <- forestcover_temperate_10[,-3]
names(forestcover_temperate_10) <- c("site", "forestcover_10")
rm(forestcover_10)

load("data/forestcover_50perc_boreal.RData")
forestcover_boreal_50 <- forestcover
names(forestcover_boreal_50) <- c("site", "forestcover_50")
load("data/forestcover_50perc_temperate.RData")
forestcover_temperate_50 <- forestcover
names(forestcover_temperate_50) <- c("site", "forestcover_50")
rm(forestcover)

forestcover_boreal <- forestcover_boreal_50 %>% left_join(forestcover_boreal_10)
forestcover_temperate <- forestcover_temperate_50 %>% left_join(forestcover_temperate_10)

# Caclualte annual disturbance rates --------------------------------------

disturbance_rates_boreal <- disturbances_boreal %>%
  purrr::map(~ freq(.) %>%
               as.data.frame(.) %>%
               right_join(., expand.grid(value = 1:14), by = "value") %>%
               mutate(disturbance = ifelse(is.na(count), 0, count)) %>%
               dplyr::select(-count)) %>%
  set_names(studysites_names_boreal) %>%
  bind_rows(.id = "site") %>%
  left_join(forestcover_boreal, by = "site") %>%
  mutate(year = value + 2000,
         forest10 = forestcover_10 / 900,
         forest50 = forestcover_50 / 900)

disturbance_rates_temperate <- disturbances_temperate %>%
  purrr::map(~ freq(.) %>%
               as.data.frame(.) %>%
               right_join(., expand.grid(value = 1:14), by = "value") %>%
               mutate(disturbance = ifelse(is.na(count), 0, count)) %>%
               dplyr::select(-count)) %>%
  set_names(studysites_names_temperate) %>%
  bind_rows(.id = "site") %>%
  left_join(forestcover_temperate, by = "site") %>%
  mutate(year = value + 2000,
         forest10 = forestcover_10 / 900,
         forest50 = forestcover_50 / 900)

disturbance_rates <- list(boreal = disturbance_rates_boreal,
                          temperate = disturbance_rates_temperate) %>%
  bind_rows(.id = "biome") %>%
  mutate(rate10 = disturbance / forest10,
         rate50 = disturbance / forest50) %>%
  filter(!is.na(disturbance))

# Add cluster information and select only northern hemisphere sites -------

cluster <- read_csv("data/studysites_combined.csv") %>%
  mutate(biome = tolower(biome)) %>%
  dplyr::select(site, cluster, x, y) %>%
  mutate(cluster = case_when(cluster == "Low" ~ 1,
                             cluster == "Moderate" ~ 2,
                             cluster == "High" ~ 3,
                             TRUE ~ 4))

disturbance_rates <- disturbance_rates %>%
  left_join(cluster, by = c("site")) %>%
  filter(y > 0 & cluster != 4)

# Save to disc ------------------------------------------------------------

write_csv(disturbance_rates, "data/disturbance_data_landscapes.csv")
