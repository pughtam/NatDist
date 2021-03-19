# Aggregate forest cover and disturbances from Hansen to definition of open/closed-canopy forests used by Pugh et al. 2019 Nature Geosciences
#
# C. Senf

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(raster)

# Get annual disturbance data ---------------------------------------------

# Studysites 

load("data/studysites_boreal.RData")
studysites_boreal <- studysites
load("data/studysites_temperate.RData")
studysites_temperate <- studysites
rm(studysites)

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

# Forests

load("data/forests_boreal.RData")
forests_boreal <- forest
load("data/forests_temperate.RData")
forests_temperate <- forest
rm(forest)

# Calculate disturbance rates for closed forests  ---------------------------

# Boreal

forests_boreal_aggregated <- forests_boreal %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(is.na(x), 0, x))) %>%
  map2(studysites_boreal, ~ mask(.x, .y)) %>%
  map(~ raster::aggregate(., fact = 40, fun = function(x, na.rm = TRUE) {
    if (mean(is.na(x)) < 0.5) { # Remove pixels with >50% NA-values
      return(mean(x, na.rm = TRUE))
    } else {
      return(NA)
    }
  })) %>%
  map(~ reclassify(., matrix(c(-1, 10, 0, 10, 50, 1, 50, 100, 2), ncol = 3, byrow = TRUE))) 

forestcover_boreal_aggregated <- forests_boreal_aggregated %>%
  map(~ cellStats(., stat = function(x, na.rm = TRUE) sum(x == 2, na.rm = TRUE) * 40*40*0.09)) %>%
  map(~ data.frame(closedforest_ha = .)) %>%
  set_names(., studysites_names_boreal) %>%
  bind_rows(.id = "site")

disturbances_boreal_closed_canopy <- forests_boreal_aggregated %>%
  map2(.y = disturbances_boreal, ~ resample(.x, .y, method = "ngb")) %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(x == 2, 1, NA))) %>%
  map2(.y = disturbances_boreal, ~ mask(.y, .x))

# Temperate

forests_temperate_aggregated <- forests_temperate %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(is.na(x), 0, x))) %>%
  map2(studysites_temperate, ~ mask(.x, .y)) %>%
  map(~ raster::aggregate(., fact = 40, fun = function(x, na.rm = TRUE) {
    if (mean(is.na(x)) < 0.5) { # Remove pixels with >50% NA-values 
      return(mean(x, na.rm = TRUE))
    } else {
      return(NA)
    }
  })) %>%
  map(~ reclassify(., matrix(c(-1, 10, 0, 10, 50, 1, 50, 100, 2), ncol = 3, byrow = TRUE))) 

forestcover_temperate_aggregated <- forests_temperate_aggregated %>%
  map(~ cellStats(., stat = function(x, na.rm = TRUE) sum(x == 2, na.rm = TRUE) * 40*40*0.09)) %>%
  map(~ data.frame(closedforest_ha = .)) %>%
  set_names(., studysites_names_temperate) %>%
  bind_rows(.id = "site")

disturbances_temperate_closed_canopy <- forests_temperate_aggregated %>%
  map2(.y = disturbances_temperate, ~ resample(.x, .y, method = "ngb")) %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(x == 2, 1, NA))) %>%
  map2(.y = disturbances_temperate, ~ mask(.y, .x))
    
### Calculate disturbance rates

disturbance_rates_boreal_closedforest <- disturbances_boreal_closed_canopy %>%
  purrr::map(~ freq(.) %>%
               as.data.frame(.) %>%
               right_join(., expand.grid(value = 1:14), by = "value") %>%
               mutate(disturbance_ha = ifelse(is.na(count), 0, count * 0.09)) %>%
               dplyr::select(-count)) %>%
  set_names(studysites_names_boreal) %>%
  bind_rows(.id = "site") %>%
  left_join(forestcover_boreal_aggregated, by = "site") %>%
  mutate(year = value + 2000)

disturbance_rates_temperate_closedforest <- disturbances_temperate_closed_canopy %>%
  purrr::map(~ freq(.) %>%
               as.data.frame(.) %>%
               right_join(., expand.grid(value = 1:14), by = "value") %>%
               mutate(disturbance_ha = ifelse(is.na(count), 0, count * 0.09)) %>%
               dplyr::select(-count)) %>%
  set_names(studysites_names_temperate) %>%
  bind_rows(.id = "site") %>%
  left_join(forestcover_temperate_aggregated, by = "site") %>%
  mutate(year = value + 2000)

### Add cluster information and select only northern hemisphere sites

disturbance_rates_closedforest <- list(disturbance_rates_boreal_closedforest, disturbance_rates_temperate_closedforest) %>%
  set_names(c("boreal", "temperate")) %>%
  bind_rows(.id = "biome")

cluster <- read_csv("data/studysites_combined.csv") %>%
  mutate(biome = tolower(biome)) %>%
  dplyr::select(site, cluster, x, y) %>%
  mutate(cluster = ifelse(cluster == "Low", 1, ifelse(cluster == "Moderate", 2, 3)))

disturbance_rates_closedforest <- disturbance_rates_closedforest %>%
  left_join(cluster, by = c("site")) %>%
  filter(y > 0)

### Save

write_csv(disturbance_rates_closedforest, "data/modeldat_closedforest.csv")


# Calculate disturbance rates for open forests  ---------------------------

# Boreal

forests_boreal_aggregated <- forests_boreal %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(is.na(x), 0, x))) %>%
  map2(studysites_boreal, ~ mask(.x, .y)) %>%
  map(~ raster::aggregate(., fact = 40, fun = function(x, na.rm = TRUE) {
    if (mean(is.na(x)) < 0.5) { # Remove pixels with >50% NA-values
      return(mean(x, na.rm = TRUE))
    } else {
      return(NA)
    }
  })) %>%
  map(~ reclassify(., matrix(c(-1, 10, 0, 10, 50, 1, 50, 100, 2), ncol = 3, byrow = TRUE))) 

forestcover_boreal_aggregated <- forests_boreal_aggregated %>%
  map(~ cellStats(., stat = function(x, na.rm = TRUE) sum(x == 1, na.rm = TRUE) * 40*40*0.09)) %>%
  map(~ data.frame(openforest_ha = .)) %>%
  set_names(., studysites_names_boreal) %>%
  bind_rows(.id = "site")

disturbances_boreal_open_canopy <- forests_boreal_aggregated %>%
  map2(.y = disturbances_boreal, ~ resample(.x, .y, method = "ngb")) %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(x == 1, 1, NA))) %>%
  map2(.y = disturbances_boreal, ~ mask(.y, .x))

# Temperate

forests_temperate_aggregated <- forests_temperate %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(is.na(x), 0, x))) %>%
  map2(studysites_temperate, ~ mask(.x, .y)) %>%
  map(~ raster::aggregate(., fact = 40, fun = function(x, na.rm = TRUE) {
    if (mean(is.na(x)) < 0.5) { # Remove pixels with >50% NA-values 
      return(mean(x, na.rm = TRUE))
    } else {
      return(NA)
    }
  })) %>%
  map(~ reclassify(., matrix(c(-1, 10, 0, 10, 50, 1, 50, 100, 2), ncol = 3, byrow = TRUE))) 

forestcover_temperate_aggregated <- forests_temperate_aggregated %>%
  map(~ cellStats(., stat = function(x, na.rm = TRUE) sum(x == 1, na.rm = TRUE) * 40*40*0.09)) %>%
  map(~ data.frame(openforest_ha = .)) %>%
  set_names(., studysites_names_temperate) %>%
  bind_rows(.id = "site")

disturbances_temperate_open_canopy <- forests_temperate_aggregated %>%
  map2(.y = disturbances_temperate, ~ resample(.x, .y, method = "ngb")) %>%
  map(~ calc(., fun = function(x, na.rm = TRUE) ifelse(x == 1, 1, NA))) %>%
  map2(.y = disturbances_temperate, ~ mask(.y, .x))

### Calculate disturbance rates

disturbance_rates_boreal_openforest <- disturbances_boreal_open_canopy %>%
  purrr::map(~ freq(.) %>%
               as.data.frame(.) %>%
               right_join(., expand.grid(value = 1:14), by = "value") %>%
               mutate(disturbance_ha = ifelse(is.na(count), 0, count * 0.09)) %>%
               dplyr::select(-count)) %>%
  set_names(studysites_names_boreal) %>%
  bind_rows(.id = "site") %>%
  left_join(forestcover_boreal_aggregated, by = "site") %>%
  mutate(year = value + 2000)

disturbance_rates_temperate_openforest <- disturbances_temperate_open_canopy %>%
  purrr::map(~ freq(.) %>%
               as.data.frame(.) %>%
               right_join(., expand.grid(value = 1:14), by = "value") %>%
               mutate(disturbance_ha = ifelse(is.na(count), 0, count * 0.09)) %>%
               dplyr::select(-count)) %>%
  set_names(studysites_names_temperate) %>%
  bind_rows(.id = "site") %>%
  left_join(forestcover_temperate_aggregated, by = "site") %>%
  mutate(year = value + 2000)

### Add cluster information and select only northern hemisphere sites

disturbance_rates_openforest <- list(disturbance_rates_boreal_openforest, disturbance_rates_temperate_openforest) %>%
  set_names(c("boreal", "temperate")) %>%
  bind_rows(.id = "biome")

cluster <- read_csv("data/studysites_combined.csv") %>%
  mutate(biome = tolower(biome)) %>%
  dplyr::select(site, cluster, x, y) %>%
  mutate(cluster = ifelse(cluster == "Low", 1, ifelse(cluster == "Moderate", 2, 3)))

disturbance_rates_openforest <- disturbance_rates_openforest %>%
  left_join(cluster, by = c("site")) %>%
  filter(y > 0)

### Save

write_csv(disturbance_rates_openforest, "data/modeldat_openforest.csv")

