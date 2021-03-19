# Code for producing Figure 3
#
# C. Senf
# 01.02.2021

# Libraries and functions -------------------------------------------------

library(tidyverse)
library(patchwork)
library(raster)
library(sf)

# Get data ----------------------------------------------------------------

regions <- list.files("results/age_structure", ".csv") %>% 
  gsub("age_dist_adjparam_latosa4000_region_", "", .) %>% 
  gsub(".csv", "", .)

dat <- list.files("results/age_structure", ".csv", full.names = TRUE) %>%
  map(read_csv, skip = 1) %>%
  set_names(regions) %>%
  bind_rows(.id = "region") %>%
  gather(key = age, value = value, -region, -Simulation) %>%
  separate("Simulation", c("simulation", "stat"), " ") %>%
  mutate(stat = gsub("\\.", "", stat)) %>%
  mutate(simulation = gsub("\\.", "", simulation)) %>%
  spread(key = stat, value = value) %>%
  mutate(age = ifelse(age == "OG", ">140", age)) %>%
  mutate(age = factor(age, levels = c(paste(seq(1, 131, 10), seq(10, 140, 10), sep = "-"), ">140")))

dat <- dat %>%
  mutate(simulation = ifelse(simulation == "LUH2", "Current", "Natural"))

# Plot --------------------------------------------------------------------

p <- ggplot() +
  geom_ribbon(data = dat %>% 
                mutate(Min = ifelse(age == ">140", NA, Min),
                       Max = ifelse(age == ">140", NA, Max)), 
              aes(x = age, 
                  ymin = Min, 
                  ymax = Max, 
                  group = interaction(simulation, region), 
                  fill = simulation), alpha = 0.3) +
  geom_line(data = dat %>% 
              mutate(Mean = ifelse(age == ">140", NA, Mean)), 
            aes(x = age, 
                y = Mean, 
                group = interaction(region, simulation), 
                col = simulation)) +
  geom_point(data = dat %>% filter(age == ">140"), 
             aes(x = age, 
                 y = Mean / 10, 
                 group = interaction(region, simulation), 
                 col = simulation, 
                 shape = simulation)) +
  geom_errorbar(data = dat %>% filter(age == ">140"), 
                aes(x = age, 
                    ymin = Min / 10, 
                    ymax = Max / 10, 
                    group = interaction(region, simulation), 
                    col = simulation),
                width = 0) +
  facet_wrap(~region, scales = "free", ncol = 4) +
  theme_linedraw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 5, color = "grey30"),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 5, color = "grey30"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(x = "Age class (years)",
       y = bquote("Young forest area (M"~km^2*")"),
       col = NULL, shape = NULL, fill = NULL) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(sec.axis = sec_axis( trans = ~.*10, name = bquote("Old-growth forest area (M"~km^2*")")))

ggsave("results/age_structure/age_structure.pdf", p, width = 7.5, height = 4.5)

regions_map <- "results/age_structure/gfad_region_map.txt" %>%
  read_table(., col_names = FALSE) %>%
  as.matrix(.) %>%
  raster(.)

regions_map[regions_map == -9999] <- NA
extent(regions_map) <- c(-180, 180, -90, 90)
projection(regions_map) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

regions_map_sf <- regions_map %>%
  rasterToPolygons(.) %>%
  st_as_sf(., ) %>%
  set_names("region", "geometry")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

for (i in unique(regions_map_sf$region)) {
  
  p <- ggplot() +
    geom_sf(data = world, color = gray(0.8), fill = gray(0.8), size = 0.01) +
    geom_sf(data = regions_map_sf %>% filter(region == i), col = NA, stroke = 0, fill = "black", size = 0.1) +
    theme_void() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    coord_sf(expand = FALSE) +
    scale_x_continuous(limits = c(-155, 170)) +
    scale_y_continuous(limits = c(-55, 85))
  
  ggsave(paste0("results/age_structure/region_map_", i, ".pdf"), p, width = 0.5, height = 0.25)
  
}


