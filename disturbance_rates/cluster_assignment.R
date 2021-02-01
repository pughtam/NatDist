
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(nnet)
library(patchwork)
library(rpart)

# Load data and wrangel into form -----------------------------------------

species <- readxl::read_excel("data/species.xlsx")
traits_pft <- readxl::read_excel("data/tree_species_traits_PFTmapping_tempbor.xlsx")
traits_pft <- traits_pft %>% dplyr::rename(., species = "Tree species", pdf = "PFT(TeBE=1,TeBS=2,IBS=3,TeNE/BNE=4,BINE=5,BNS=6)")
wooddensity <- data.frame(pdf_name = c("TeBE", "TeBS", "IBS", "TeNE/BNE", "BINE", "BNS"),
                          pdf = c(1, 2, 3, 4, 5, 6),
                          wooddensity_mean = c(298.995, 351.514, 262.939, 229.419, 224.408, 263.705),
                          wooddensity_sd = c(79.292, 41.564, 38.547, 20.772, 16.646, 13.959))
maximumheight <- data.frame(pdf_name = c("TeBE", "TeBS", "IBS", "TeNE/BNE", "BINE", "BNS"),
                            pdf = c(1, 2, 3, 4, 5, 6),
                            maximumheight_mean = c(38.906, 20.781, 25.871, 46.356, 45.367, 54.605),
                            maximumheight_sd = c(11.777, 14.445,  7.966, 12.202, 11.961,  5.030))
traits_pft <- traits_pft %>% left_join(wooddensity) %>% left_join(maximumheight)

traits_species <- readxl::read_excel("data/traits.xlsx")

cluster <- read_csv("data/studysites_combined.csv") %>%
  mutate(biome = tolower(biome))

load("data/climate_variability.RData")

climate <- climate_variability %>%
  filter(year %in% 1980:2010) %>%
  group_by(site, biome) %>%
  summarize(temp_mean_sd = sd(temp_mean, na.rm = TRUE),
            prec_mean_sd = sd(prec_mean, na.rm = TRUE),
            temp_mean = mean(temp_mean, na.rm = TRUE),
            temp_range = mean(temp_range, na.rm = TRUE),
            prec_mean = mean(prec_mean, na.rm = TRUE),
            prec_range = mean(prec_range, na.rm = TRUE))

dat_landscape <- species %>%
  group_by(site, biome, tree_num) %>% 
  summarise(species = unique(tree),
            tree_perc = mean(tree_perc)) %>%
  mutate(genus = stringr::str_split(species, " ") %>% map(~.[[1]]) %>% unlist(.)) %>%
  left_join(traits_pft %>% dplyr::select(species, 
                                         broadleaved.pft = broadleaved, 
                                         wooddensity.pft = wooddensity_mean, 
                                         maximumheight.pft = maximumheight_mean), 
            by = "species") %>%
  left_join(traits_species %>% dplyr::select(species, 
                                            coniferous.species = conifer, 
                                            wooddensity.species = wd, 
                                            maximumheight.species = hmax90), by = "species") %>%
  split(.$genus) %>% # Fill missing species trait values with genus averages
  map(., ~ mutate(., coniferous.species = ifelse(is.na(coniferous.species), raster::modal(.$coniferous.species, na.rm = TRUE), coniferous.species),
                  wooddensity.species = ifelse(is.na(wooddensity.species), mean(.$wooddensity.species, na.rm = TRUE), wooddensity.species),
                  maximumheight.species = ifelse(is.na(maximumheight.species), mean(.$maximumheight.species, na.rm = TRUE), maximumheight.species))) %>% 
  bind_rows(.) %>%
  na.omit(.) %>%
  group_by(site, biome) %>%
  mutate(weight = tree_perc / sum(tree_perc)) %>%
  summarize(height.pft = sum(maximumheight.pft * weight, na.rm = TRUE),
            wooddensity.pft = sum(wooddensity.pft * weight, na.rm = TRUE),
            conifershare.pft = mean(broadleaved.pft == "no", na.rm = TRUE),
            height.species = sum(maximumheight.species * weight, na.rm = TRUE),
            wooddensity.species = sum(wooddensity.species * 1000 * weight, na.rm = TRUE),
            conifershare.species = mean(coniferous.species == 1, na.rm = TRUE)) %>%
  left_join(climate, by = c("site", "biome")) %>%
  left_join(cluster, by = c("site", "biome")) %>%
  na.omit(.)

dat_landscape <- dat_landscape %>%
  filter(cluster != "Ephemeral") %>%
  filter(y > 0)

write_csv(dat_landscape, "data/traits_landscape.csv")

# Model -------------------------------------------------------------------

dat_landscape <- read_csv("data/traits_landscape.csv")

dat_landscape %>%
  dplyr::select(-x, -y, -site_id) %>%
  gather(key = trait, value = value, -site, -biome, -cluster) %>%
  separate("trait", c("trait", "source"), "\\.") %>%
  ggplot(., aes(x = cluster, y = value, fill = source)) +
  geom_boxplot() +
  facet_wrap(~trait, ncol = 5, scales = "free")

# Model using PFT trait values

multinomModel_null_pft <- multinom(cluster ~ 1, data = dat_landscape)
multinomModel1_pft <- multinom(cluster ~ height.pft + wooddensity.pft + conifershare.pft + temp_mean + temp_range + prec_mean + prec_range, data = dat_landscape)
multinomModel2_pft <- multinom(cluster ~ wooddensity.pft + temp_range, data = dat_landscape)
multinomModel3_pft <- multinom(cluster ~ wooddensity.pft + temp_range + height.pft, data = dat_landscape)
multinomModel4_pft <- multinom(cluster ~ wooddensity.pft + temp_range + prec_range + height.pft, data = dat_landscape)
multinomModel5_pft <- multinom(cluster ~ wooddensity.pft + temp_range + conifershare.pft, data = dat_landscape)
multinomModel6_pft <- multinom(cluster ~ wooddensity.pft + temp_mean, data = dat_landscape)
multinomModel7_pft <- multinom(cluster ~ wooddensity.pft + temp_mean + height.pft, data = dat_landscape)
multinomModel8_pft <- multinom(cluster ~ temp_mean + height.pft, data = dat_landscape)
multinomModel9_pft <- multinom(cluster ~ temp_mean, data = dat_landscape)

mod_list <- list(multinomModel1_pft,
                 multinomModel2_pft,
                 multinomModel3_pft,
                 multinomModel4_pft,
                 multinomModel5_pft,
                 multinomModel6_pft,
                 multinomModel7_pft,
                 multinomModel8_pft,
                 multinomModel9_pft)

aics <- mod_list %>% 
  map(AIC) %>%
  unlist() %>%
  data.frame(aic = .)

preds <- mod_list %>%
  map(~ .$coefnames %>% paste(., collapse = ",")) %>%
  unlist()

aics$model <- preds  

aics <- aics %>% arrange(., aic)

write_csv(aics, "results/cluster_assignement/cluster_model_comparison.csv")

final_model_pft <- multinomModel2_pft

summary(final_model_pft)

anova(multinomModel_null_pft, final_model_pft)

conf <- table(factor(predict(final_model_pft, dat_landscape), levels = c("Low", "Moderate", "High")), 
              factor(dat_landscape$cluster, levels = c("Low", "Moderate", "High")))
sum(diag(conf)) / sum(conf)
1 - diag(conf) / rowSums(conf)
1 - diag(conf) / colSums(conf)

coefs_pft <- coef(final_model_pft) %>%
  as.data.frame() %>%
  rownames_to_column(var = "cluster")

write_csv(as.data.frame(coefs_pft), "results/cluster_assignement/multi_mod_coefficients_pft-traits.csv")

# Model using species trait values

multinomModel_null_species <- multinom(cluster ~ 1, data = na.omit(dat_landscape))
multinomModel1_species <- multinom(cluster ~ height.species + wooddensity.species + conifershare.species + temp_mean + temp_range + prec_mean + prec_range, data = na.omit(dat_landscape))
multinomModel2_species <- multinom(cluster ~ wooddensity.species + temp_range, data = na.omit(dat_landscape))
multinomModel3_species <- multinom(cluster ~ wooddensity.species + temp_range + height.species, data = na.omit(dat_landscape))
multinomModel4_species <- multinom(cluster ~ wooddensity.species + temp_range + prec_range + height.species, data = na.omit(dat_landscape))
multinomModel5_species <- multinom(cluster ~ wooddensity.species + temp_range + conifershare.species, data = na.omit(dat_landscape))
multinomModel6_species <- multinom(cluster ~ wooddensity.species + temp_mean, data = na.omit(dat_landscape))
multinomModel7_species <- multinom(cluster ~ wooddensity.species + temp_mean + height.species, data = na.omit(dat_landscape))
multinomModel8_species <- multinom(cluster ~ temp_mean + height.species, data = na.omit(dat_landscape))

mod_list <- list(multinomModel1_species,
                 multinomModel2_species,
                 multinomModel3_species,
                 multinomModel4_species,
                 multinomModel5_species,
                 multinomModel6_species,
                 multinomModel7_species,
                 multinomModel8_species)

aics <- mod_list %>% 
  map(AIC) %>%
  unlist() %>%
  data.frame(aic = .)

preds <- mod_list %>%
  map(~ .$coefnames %>% paste(., collapse = ",")) %>%
  unlist()

aics$model <- preds  

aics <- aics %>% arrange(., aic)

write_csv(aics, "results/cluster_assignement/cluster_model_comparison.csv")

final_model_species <- multinomModel2_species

summary(final_model_species)

anova(multinomModel_null_species, final_model_species)

conf <- table(factor(predict(final_model_species), levels = c("Low", "Moderate", "High")), 
              factor(na.omit(dat_landscape)$cluster, levels = c("Low", "Moderate", "High")))
sum(diag(conf)) / sum(conf)
1 - diag(conf) / rowSums(conf) # Commission
1 - diag(conf) / colSums(conf) # Omission

coefs_species <- coef(final_model_species) %>%
  as.data.frame() %>%
  rownames_to_column(var = "cluster")

write_csv(as.data.frame(coefs_species), "results/cluster_assignement/multi_mod_coefficients_species-traits.csv")

# Create plots ------------------------------------------------------------

# Using PFT trait values

disturbance_rates <- read_csv("results/disturbance_rates_summary.csv")

newdata_pft <- expand.grid(wooddensity.pft = seq(quantile(dat_landscape$wooddensity.pft, 0, na.rm = TRUE),
                                                 quantile(dat_landscape$wooddensity.pft, 1, na.rm = TRUE), length.out = 100),
                           height.pft = c(quantile(dat_landscape$height.pft, 0, na.rm = TRUE),
                                          quantile(dat_landscape$height.pft, 0.5, na.rm = TRUE),
                                          quantile(dat_landscape$height.pft, 1, na.rm = TRUE)),
                           temp_range = seq(quantile(dat_landscape$temp_range, 0, na.rm = TRUE),
                                            quantile(dat_landscape$temp_range, 1, na.rm = TRUE), length.out = 100))

newdata_pft <- cbind(newdata_pft, predict(final_model_pft, newdata = newdata_pft, type = "prob"))

newdata_pft <- cbind(newdata_pft, 
                     mean_dist_rate = as.matrix(newdata_pft[, c("Low", "Moderate", "High")]) %*% disturbance_rates$mean)

#1 - cluster assignment

p <- ggplot(newdata_pft) +
  geom_raster(aes(x = wooddensity.pft, y = temp_range, 
                fill = rgb(High, Low, Moderate, maxColorValue = 1))) +
  scale_fill_identity() +
  labs(x = bquote("Wood density (kg m"^-3~")"), 
       y = "Temperature range (°C)",
       title = "Height (m)") +
  theme_linedraw() +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(color = "grey20", size = 8),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  facet_wrap(~height.pft)

ggsave("results/cluster_assignement/cluster_attribution_pft-traits.pdf", p, width = 7.5, height = 2.5)

#2 - mean disturbance rates

p <- ggplot(newdata_pft) +
  geom_raster(aes(x = wooddensity.pft, y = temp_range, 
                  fill = mean_dist_rate * 100)) +
  labs(x = bquote("Wood density (kg m"^-3~")"), 
       y = "Temperature range (°C)",
       fill = "Average disturbance rate (% per yr.)",
       title = "Height (m)") +
  theme_linedraw() +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(color = "grey20", size = 8),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(title.position = "top")) +
  facet_wrap(~height.pft)

ggsave("results/cluster_assignement/cluster_attribution_disturbance_rates_pft-traits.pdf", p, width = 7.5, height = 3)

# Using species trait values

disturbance_rates <- read_csv("results/disturbance_rates_summary.csv")

newdata_species <- expand.grid(wooddensity.species = seq(quantile(dat_landscape$wooddensity.species, 0, na.rm = TRUE),
                                                 quantile(dat_landscape$wooddensity.species, 1, na.rm = TRUE), length.out = 100),
                           temp_range = seq(quantile(dat_landscape$temp_range, 0, na.rm = TRUE),
                                            quantile(dat_landscape$temp_range, 1, na.rm = TRUE), length.out = 100))

newdata_species <- cbind(newdata_species, predict(final_model_species, newdata = newdata_species, type = "prob"))

newdata_species <- cbind(newdata_species, 
                     mean_dist_rate = as.matrix(newdata_species[, c("Low", "Moderate", "High")]) %*% disturbance_rates$mean)

#1 - cluster assignment

p <- ggplot(newdata_species) +
  geom_raster(aes(x = wooddensity.species, y = temp_range, 
                  fill = rgb(High, Low, Moderate, maxColorValue = 1))) +
  scale_fill_identity() +
  labs(x = bquote("Wood density (kg m"^-3~")"), 
       y = "Temperature range (°C)") +
  theme_linedraw() +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(color = "grey20", size = 8),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001))

ggsave("results/cluster_assignement/cluster_attribution_species-traits.pdf", p, 
       width = 2, height = 2)

#2 - mean disturbance rates

p <- ggplot(newdata_species) +
  geom_raster(aes(x = wooddensity.species, 
                  y = temp_range, 
                  fill = log10(mean_dist_rate * 100))) +
  labs(x = bquote("Wood density (kg m"^-3~")"), 
       y = "Temperature range (°C)",
       fill = bquote("%"~yr^-1)) +
  theme_linedraw() +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.15, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(color = "grey20", size = 8),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_viridis_c(limits = c(log10(0.03), log10(1.25)), 
                       breaks = c(log10(0.05), log10(0.1), log10(0.2), log10(0.4), log10(0.8)), labels = c(0.05, 0.1, 0.2, 0.4, 0.8)) +
  guides(fill = guide_colorbar(title.position = "top"))

ggsave("results/cluster_assignement/cluster_attribution_disturbance_rates_species-traits.pdf", p, width = 2.55, height = 2)

### Map of landscapes

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

sites <- read_csv("data/studysites_combined.csv")
sites <- sites %>% 
  filter(site %in% unique(dat_landscape$site)) %>%
  st_as_sf(., coords = c("x", "y"))
st_crs(sites) <- st_crs(countries)

p <- ggplot() +
  geom_sf(data = countries, color = gray(0.3), fill = gray(0.8)) +
  geom_sf(data = sites, aes(col = factor(cluster, levels = c("Low", "Moderate", "High")),
                            shape = factor(cluster, levels = c("Low", "Moderate", "High")))) +
  theme_void() +
  scale_x_continuous(limits = c(-155, 170)) +
  scale_y_continuous(limits = c(-50, 78)) +
  theme(legend.position = c(0, 0),
        legend.justification = c(-0.3, -0.3),
        legend.background = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  scale_color_manual(values = c("green", "blue", "red")) +
  labs(col = "Cluster", shape = "Cluster")

ggsave("results/cluster_assignement/sites_cluster_map.pdf", p, width = 5, height = 3)

