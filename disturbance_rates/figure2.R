
# Libraries and functions -------------------------------------------------

library(raster)
library(tidyverse)
library(sf)
library(lme4)
library(geometry)
library(patchwork)

# General data ------------------------------------------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Landscapes with return-intervalls ---------------------------------------

modeldat <- read_csv("data/modeldat.csv")
site_locations <- read_csv("data/studysites_combined.csv")
#site_shp <- read_sf("../Boreal Forests/data/studysites/studysites_combined_centers.shp")

# Read LPJ-Output ---------------------------------------------------------

## rotation period

rotation <- "results/lpj_outputs/best_est_adjparam_latosa4000_20patch_10pCanopyCover.txt" %>%
  read_table(., col_names = FALSE) %>%
  as.matrix(.) %>%
  raster(.)

rotation[rotation == -9999] <- NA
extent(rotation) <- c(-180, 180, -90, 90)
projection(rotation) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

plot(rotation)

## Rotation period in closed canopy forests

rotation_closedcanopy <- "results/lpj_outputs/best_est_adjparam_latosa4000_closedcan_20patch_5pClosedCanopyCover_1deg.txt" %>%
  read_table(., col_names = FALSE) %>%
  as.matrix(.) %>%
  raster(.)

rotation_closedcanopy[rotation_closedcanopy == -9999] <- NA
extent(rotation_closedcanopy) <- c(-180, 180, -90, 90)
projection(rotation_closedcanopy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

plot(rotation_closedcanopy)

## Observed rotation periods (Pugh et al. 2019)

rotation_pugh2019 <- "results/pugh2019_natgeo/tauO_standard_forest-area_LUcorrected.nc" %>%
  raster(.)

## Forest cover

forestfraction <- "results/pugh2019_natgeo/hansen_forested_canopy_frac_0p5deg.nc4" %>%
  raster(.)

## Biomes

biomes <- "results/lpj_outputs/simplemodel_best_est_100patch_10pCanopyCover_biomes_v2.txt" %>%
  read_table(., col_names = FALSE) %>%
  as.matrix(.) %>%
  raster(.)

biomes[biomes == -9999] <- NA
extent(biomes) <- c(-180, 180, -90, 90)
projection(biomes) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

plot(biomes)


# Convex hull -------------------------------------------------------------

dat_landscape <- read_csv("data/traits_landscape.csv")

traits_lpj <- "results/lpj_outputs/distprob_2001_2014" %>%
  read_table()

traits_lpj_raster <- rasterFromXYZ(traits_lpj, crs = projection(rotation))

chull <- convhulln(dat_landscape %>% 
                     dplyr::select(wooddensity = wooddensity.species, temp_range) %>%
                     as.matrix())

inhull <- inhulln(ch = chull, 
                  p = traits_lpj %>% 
                    dplyr::select(wooddensity = WDmean, temp_range = Trange) %>%
                    mutate(wooddensity = wooddensity * 1000) %>%
                    as.matrix(.))

inhull <- traits_lpj %>% 
  dplyr::select(Lon, Lat) %>%
  mutate(inhull = inhull)

inhull_sf <- inhull %>%
  st_as_sf(., coords = c("Lon", "Lat"), crs = projection(rotation))

p_chull <- ggplot() +
  geom_sf(data = world, color = gray(0.8), fill = gray(0.8)) +
  geom_point(data = inhull, aes(x = Lon, y = Lat, col = ifelse(inhull, "Inside", "Outside")), alpha = 0.5, size = 0.1, pch = 19) +
  theme_void() +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, 
       col = NULL,
       fill = NULL) +
  xlim(-170, 170) +
  ylim(min(inhull$Lat), max(inhull$Lat))

ggsave("results/convex_hull.pdf", p_chull, width = 7.5, height = 3)

rotation_sf <- rotation %>%
  rasterToPolygons(.) %>%
  st_as_sf(., ) %>%
  set_names("rotation", "geometry")

rotation_inhull <- rotation_sf %>%
  st_join(., inhull_sf)

p_chull <- rotation_inhull %>%
  st_drop_geometry() %>%
  ggplot(data = .,  aes(x = rotation, 
                        col = factor(inhull, labels = c("Extrapolation", "Interpolation")))) +
  stat_density(geom = "line", position = "identity") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1)) +
  labs(x = bquote("Disturbance rotation period"~(tau)~"for all forests"),
       y = "Density",
       col = NULL) +
  xlim(0, 1000)

ggsave("results/rotation_period_chull_cutat1000.pdf", p_chull, width = 4, height = 3)

rotation_closedcanopy_1deg <- aggregate(rotation_closedcanopy, fact = 2, fun = mean)

rotation_closedcanopy_diff <- rotation_closedcanopy_1deg - rotation_pugh2019

#rotation_closedcanopy_diff <- disaggregate(rotation_closedcanopy_diff, fact = 2)

rotation_diff_sf <- rotation_closedcanopy_diff %>%
  rasterToPolygons(.) %>%
  st_as_sf(., ) %>%
  set_names("rotation_diff", "geometry")

rotation_diff_inhull <- rotation_diff_sf %>%
  st_join(., inhull_sf)

p_chull <- rotation_diff_inhull %>%
  st_drop_geometry() %>%
  filter(!is.na(inhull)) %>%
  ggplot(data = .,  aes(x = rotation_diff, col = factor(inhull))) +
  stat_density(geom = "line", position = "identity") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1)) +
  labs(x = bquote("Disturbance rotation period"~(tau)~"for all forests"),
       y = "Density",
       col = NULL)

ggsave("results/rotation_period_chull.pdf", p_chull, width = 3.5, height = 3.5)

# Figure 2 -----------------------------------------------------------------

### Panel A

rotation_sf <- rotation %>%
  rasterToPolygons(.) %>%
  st_as_sf(., ) %>%
  set_names("rotation", "geometry") %>%
  mutate(rotation_capped = ifelse(rotation > 1000, 1000, rotation),
         rotation_log10 = log10(rotation))

forestfraction_sf <- forestfraction %>%
  rasterToPolygons(.) %>%
  st_as_sf(., ) %>%
  set_names("forestcover", "geometry")

rotation_sf <- rotation_sf %>%
  st_join(., forestfraction_sf)

rotation_sf <- rotation_sf %>%
  mutate(rotation_log10_forest = ifelse(forestcover > 10, rotation_log10, NA))

panel_a <- ggplot() +
  geom_sf(data = world, color = gray(0.8), fill = gray(0.8)) +
  geom_sf(data = rotation_sf,
          aes(fill = log10(rotation_capped),
              col = log10(rotation_capped)), size = 0.2) +
  scale_color_viridis_c(direction = -1, limits = c(2, 3), breaks = c(2, 2.25, 2.5, 2.75, 3), labels = round(10^c(2, 2.25, 2.5, 2.75, 3), 0)) +
  scale_fill_viridis_c(direction = -1, limits = c(2, 3), breaks = c(2, 2.25, 2.5, 2.75, 3), labels = round(10^c(2, 2.25, 2.5, 2.75, 3), 0)) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "lightblue", color = NULL),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(0.125, "cm"),
        legend.key.height = unit(0.35, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, 
       fill = bquote(log[10](tau)),
       col = bquote(log[10](tau)),
       # col = NULL,
       # fill = NULL,
       title = bquote(bold(.("A"))~"Disturbance rotation period"~(tau)~"for all forests")) +
  xlim(-170, 170) +
  ylim(17.5, 72.5)

cutter <- 1500

panel_a_nonlog <- ggplot() +
  geom_sf(data = world, color = gray(0.8), fill = gray(0.8)) +
  geom_sf(data = rotation_sf,
          aes(fill = rotation_capped,
              col = rotation_capped), size = 0.2) +
  # geom_point(data = inhull %>% filter(inhull == FALSE), 
  #            aes(x = Lon, y = Lat), alpha = 0.1, size = 0.05, pch = 19, col = "black") +
  scale_color_viridis_c(direction = -1, rescaler = function(x, to = c(0, 1), from = NULL) {
    ifelse(x < cutter, 
           scales::rescale(x,
                           to = to,
                           from = c(min(x, na.rm = TRUE), cutter)),
           1)}) +
  scale_fill_viridis_c(direction = -1, rescaler = function(x, to = c(0, 1), from = NULL) {
    ifelse(x < cutter, 
           scales::rescale(x,
                           to = to,
                           from = c(min(x, na.rm = TRUE), cutter)),
           1)}) +
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(0.125, "cm"),
        legend.key.height = unit(0.35, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, 
       col = bquote(tau),
       fill = bquote(tau)) +
  scale_x_continuous(limits = c(-155, 170)) +
  scale_y_continuous(limits = c(-55, 85))

### Panel B

rotation_closedcanopy[rotation_closedcanopy > 1000] <- 1000

rotation_closedcanopy_diff <- rotation_closedcanopy - rotation_pugh2019

rotation_diff_sf <- rotation_closedcanopy_diff %>%
  rasterToPolygons(.) %>%
  st_as_sf(., ) %>%
  set_names("rotation_diff", "geometry")

panel_b <- ggplot() +
  geom_sf(data = world, color = gray(0.8), fill = gray(0.8)) +
  geom_sf(data = rotation_diff_sf,
          aes(fill = rotation_diff,
              col = rotation_diff), size = 0.2) +
  # geom_point(data = inhull %>% filter(inhull == FALSE), 
  #            aes(x = Lon, y = Lat), alpha = 0.1, size = 0.05, pch = 19, col = "black") +
  scale_fill_gradient2(limits = c(-1000, 1000)) +
  scale_color_gradient2(limits = c(-1000, 1000)) +
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(0.125, "cm"),
        legend.key.height = unit(0.35, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, 
       fill = bquote(Delta*tau),
       col = bquote(Delta*tau)) +
  scale_x_continuous(limits = c(-155, 170)) +
  scale_y_continuous(limits = c(-55, 85))

### Panel C 

rotation_data <- data.frame(rotation = values(rotation)) %>%
  na.omit() %>%
  mutate(rotation = ifelse(rotation > 1000, 1000, rotation)) %>%
  arrange(rotation) %>%
  mutate(p = 1/n(),
         pcum = cumsum(p))

rotation_pugh2019_data <- data.frame(rotation = values(rotation_pugh2019)) %>%
  na.omit() %>%
  arrange(rotation) %>%
  mutate(p = 1/n(),
         pcum = cumsum(p))

rotation_closedcanopy_data <- data.frame(rotation = values(rotation_closedcanopy)) %>%
  na.omit() %>%
  mutate(rotation = ifelse(rotation > 1000, 1000, rotation)) %>%
  arrange(rotation) %>%
  mutate(p = 1/n(),
         pcum = cumsum(p))

panel_c_1.1 <- ggplot() +
  geom_step(data = rotation_data, aes(x = rotation, y = pcum, col = "All forests")) +
  labs(x = bquote(tau), 
       y = "Empirical PDF",
       col = NULL,
       title = bquote(bold(.("B")))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = -0.4),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        legend.key.size = unit(0.125, "cm"),
        legend.text = element_text(size = 7)) +
  scale_color_manual(values = "black")

panel_c_1.2 <- ggplot() +
  geom_histogram(data = rotation_data, aes(x = rotation)) +
  labs(x = bquote(tau), 
       y = "Counts",
       col = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = -0.4),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.size = unit(0.125, "cm"),
        legend.text = element_text(size = 7)) +
  scale_color_manual(values = "black")

panel_c_2.1 <- ggplot() +
  geom_step(data = rotation_pugh2019_data, aes(x = rotation, y = pcum, col = "Current")) +
  geom_step(data = rotation_closedcanopy_data, aes(x = rotation, y = pcum, col = "Natural")) +
  labs(x = bquote(tau), 
       y = "Empirical PDF",
       col = "Closed canopy",
       title = bquote(bold(.("D")))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = -0.4),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        legend.key.size = unit(0.125, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  scale_color_brewer(palette = "Set1")

panel_c_2.2 <- ggplot() +
  geom_histogram(data = rotation_closedcanopy_data, aes(x = rotation, fill = "Natural"), alpha = 0.5) +
  geom_histogram(data = rotation_pugh2019_data, aes(x = rotation, fill = "Current"), alpha = 0.5) +
  labs(x = bquote(tau), 
       y = "Counts",
       col = NULL,
       fill = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = -0.4),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.size = unit(0.125, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  scale_fill_brewer(palette = "Set1")


### Combine all graphs

p <- panel_a_nonlog + 
  plot_spacer() + 
  (plot_spacer() + panel_c_1.2 + plot_spacer() + plot_layout(ncol = 1, heights = c(0.01, 0.98, 0.01))) + 
  panel_b + 
  plot_spacer() + 
  (plot_spacer() + panel_c_2.2 + plot_spacer() + plot_layout(ncol = 1, heights = c(0.01, 0.98, 0.01))) + 
  plot_layout(ncol = 3, widths = c(0.8, 0.01, 0.19), height = c(0.5, 0.5)) +
  plot_annotation(tag_levels = "A")
  
ggsave("results/figure2.pdf", p, width = 7.5, height = 5)

# Combine with Curtis -----------------------------------------------------

curtis <- raster("data/Curtis_Drivers_Deforest_wgs84.tif")
rotation_grid <- rasterToPolygons(rotation)
rotation_grid$layer <- 1:ncell(rotation_grid)

rotation_closedcanopy_diff_curtis <- resample(rotation_closedcanopy_diff, curtis, method = "ngb")

rotation_closedcanopy_diff_curtis <- stack(rotation_closedcanopy_diff_curtis,
                                           curtis) %>%
  as.data.frame() %>%
  set_names(., c("rotation_diff", "agent")) %>%
  na.omit() %>%
  mutate(agent = factor(agent, levels = 0:5 ,
                        labels = c("Minor loss",
                                   "Deforestation", 
                                   "Shifting agriculture",
                                   "Forestry",
                                   "Wildfire",
                                   "Urbanization"))) %>%
  mutate(rotation_diff_cap = ifelse(rotation_diff > 1000, 1000, rotation_diff))

p <- ggplot(rotation_closedcanopy_diff_curtis, 
       aes(x = reorder(agent, rotation_diff_cap), 
           y = rotation_diff_cap,
           fill = agent)) +
  geom_boxplot(outlier.colour = NA) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Set1")[c(1, 5, 6, 3, 4, 2)]) +
  labs(x = NULL, 
       y = bquote("Difference to observed"~tau)) +
  theme(legend.position = "none")

ggsave("results/curtis_difference_rotation.pdf", p, width = 3.5, height = 3)

# Validation ----------------------------------------------------

sites_xy <- modeldat %>% filter(y > 0 & cluster != "Ephemeral")
sites_xy <- unique(as.matrix(sites_xy[, c("x", "y")]))

biomes_df <- biomes %>%
  rasterToPolygons(.)
biomes_df <- cbind(coordinates(biomes_df), biomes_df$layer)
colnames(biomes_df) <- c("x", "y", "biomes")

closest <- RANN::nn2(data = biomes_df[, c("x", "y")], 
                     query = sites_xy, 
                     k = 1)

biomes_extr <- biomes_df[closest$nn.idx, "biomes"]

observed_rotations <- modeldat %>%
  filter(., y > 0 & cluster != "Ephemeral") %>%
  group_by(site) %>%
  summarise(distrate_observed = mean(rate10),
            rotation_observed = rotperiod(distrate_observed)) %>%
  ungroup() %>%
  mutate(biomes = biomes_extr) %>%
  mutate(rotation_observed_capped = ifelse(rotation_observed > 1000, 1000, rotation_observed)) %>%
  left_join(modeldat %>% group_by(site) %>% summarise(cluster = factor(unique(cluster))), by = "site")

ggplot(observed_rotations, aes(x = reorder(site, rotation_observed), y = rotation_observed)) +
  geom_point(aes(col = cluster)) +
  coord_flip() +
  scale_y_log10()

ggplot(observed_rotations, aes(x = cluster, y = rotation_observed)) +
  geom_boxplot() +
  scale_y_log10()

ggplot(observed_rotations, aes(x = factor(biomes), y = rotation_observed)) +
  geom_boxplot() +
  scale_y_log10()

dat_obs_summary <- dat_obs %>%
  group_by(biomes) %>%
  summarize(obs_mn = mean(obs),
            obs_capped_mn = mean(obs_capped)) %>%
  ungroup()

ggplot(dat_obs, aes(x = factor(biomes), y = obs)) +
  geom_boxplot()

dat_mod <- data.frame(biomes = values(biomes),
                      mod = values(rotation),
                      mod_capped = ifelse(values(rotation) > 1000, 1000, values(rotation))) 

dat_mod_summary <- dat_mod %>%
  group_by(biomes) %>%
  summarize(mod_mn = mean(mod, na.rm = TRUE),
            mod_capped_mn = mean(mod_capped, na.rm = TRUE)) %>%
  ungroup()

ggplot(dat_mod, aes(x = factor(biomes), y = mod)) +
  geom_boxplot()

dat_obs_summary %>%
  left_join(dat_mod_summary) %>%
  ggplot(.) +
  geom_point(aes(x = mod_mn, y = obs_mn)) +
  # geom_errorbarh(aes(y = obs_mn, xmin = mod_mn - mod_sd, xmax = mod_mn + mod_sd)) +
  # geom_errorbar(aes(x = mod_mn, ymin = obs_mn - obs_sd, ymax = obs_mn + obs_sd)) +
  #xlim(0, 1250) + ylim(0, 1250) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(aes(x = mod_capped_mn, y = obs_capped_mn), method = "lm", se = FALSE)

dat_obs %>%
  right_join(dat_mod %>% na.omit(.)) %>%
  gather(key = key, value = value, -site, -cluster, -biomes) %>%
  ggplot(., aes(x = factor(biomes), y = value, fill = key)) +
  geom_boxplot() +
  scale_y_log10()

ggplot(dat_obs, aes(x = obs_capped, y = mod_capped)) +
  geom_point() +
  facet_wrap(~cluster) +
  geom_abline(intercept = 0, slope = 1)
