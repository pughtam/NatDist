# Code for producing Figure 2 based on LPJ-Output
#
# C. Senf
# 01.02.2021

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

modeldat <- read_csv("../disturbance_rates/data/modeldat_openforest.csv")
site_locations <- read_csv("../disturbance_rates/data/studysites_combined.csv")

# Read LPJ-Output ---------------------------------------------------------

## rotation period

rotation <- "best_est_adjparam_100patch_10pCanopyCover.txt" %>%
  read_table(., col_names = FALSE) %>%
  as.matrix(.) %>%
  raster(.)

rotation[rotation == -9999] <- NA
extent(rotation) <- c(-180, 180, -90, 90)
projection(rotation) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

plot(rotation)

## Rotation period in closed canopy forests

rotation_closedcanopy <- "best_est_adjparam_closedcan_20patch_5pClosedCanopyCover_1deg.txt" %>%
  read_table(., col_names = FALSE) %>%
  as.matrix(.) %>%
  raster(.)

rotation_closedcanopy[rotation_closedcanopy == -9999] <- NA
extent(rotation_closedcanopy) <- c(-180, 180, -90, 90)
projection(rotation_closedcanopy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

plot(rotation_closedcanopy)

## Observed rotation periods (Pugh et al. 2019)

rotation_pugh2019 <- "/Users/pughtam/Documents/GAP_and_other_work/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc" %>%
  raster(.)

## Forest cover

forestfraction <- "/Users/pughtam/data/hansen/hansen_canopy_cover_frac.nc4" %>%
  raster(.)


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
       title = bquote(bold(.("A"))~"Disturbance rotation period"~(tau)~"for all forests")) +
  xlim(-170, 170) +
  ylim(17.5, 72.5)

cutter <- 1500

panel_a_nonlog <- ggplot() +
  geom_sf(data = world, color = gray(0.8), fill = gray(0.8)) +
  geom_sf(data = rotation_sf,
          aes(fill = rotation_capped,
              col = rotation_capped), size = 0.2) +
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
  
ggsave("figure2.pdf", p, width = 7.5, height = 5)

