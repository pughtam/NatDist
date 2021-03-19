# Calibrate disturbance rate model for closed canopy forests
#
# C. Senf
# 19.02.2021

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(raster)
library(lme4)
library(patchwork)

### Functions

inv_logit <- boot::inv.logit
logit <- boot::logit

# Load data and wrangel into form -----------------------------------------

modeldat_closedforests <- read_csv("data/modeldat_closedforest.csv")
modeldat_closedforests$rate <- modeldat_closedforests$disturbance_ha / modeldat_closedforests$closedforest_ha
modeldat_closedforests$rate_model <- round(modeldat_closedforests$disturbance_ha, 0) / modeldat_closedforests$closedforest_ha
modeldat_closedforests$cluster <- factor(modeldat_closedforests$cluster)
modeldat_closedforests$forest <- as.integer(modeldat_closedforests$closedforest_ha)

modeldat <- read_csv("data/modeldat.csv") %>%
  dplyr::select(site, year, rate10, rate50, disturbance, forestcover_10, forestcover_50)

no_closed_forest <- modeldat_closedforests %>%
  group_by(site) %>%
  summarise(ncf = sum(closedforest_ha) == 0) %>%
  filter(ncf == TRUE)

modeldat_closedforests <- modeldat_closedforests %>%
  filter(!(site %in% no_closed_forest$site))

# Compare rates -----------------------------------------------------------

modeldat_compare <- modeldat_closedforests %>%
  left_join(modeldat)

ggplot(modeldat_compare, aes(x = rate, y = rate10, col = cluster)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Closed forest rate", y = "All forest rate") +
  facet_wrap(~biome)

modeldat_compare %>%
  mutate(diff = rate - rate10) %>%
  ggplot(.) +
  geom_histogram(aes(x = diff)) +
  facet_wrap(~biome)

modeldat_compare %>%
  mutate(diff = rate - rate10) %>%
  group_by(biome) %>%
  summarize(mean = mean(diff, na.rm = TRUE),
            median = median(diff, na.rm = TRUE),
            min = min(diff, na.rm = TRUE),
            max = max(diff, na.rm = TRUE))

modeldat_compare %>%
  group_by(site) %>%
  summarize(mn_allforest = mean(rate10, na.rm = TRUE),
            mn_closedforest = mean(rate, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(mn_allforest = mean(mn_allforest, na.rm = TRUE),
            mn_closedforest = mean(mn_closedforest, na.rm = TRUE))

# Fit model ---------------------------------------------------------------

fit_null_null <- glm(rate_model ~ 1,
                     data = modeldat_closedforests,
                     weight = forest,
                     family = binomial(link = "logit"))

fit_null <- glmer(rate_model ~ 1 + (1 | site) + (1 | year),
                  data = modeldat_closedforests,
                  weight = forest,
                  family = binomial(link = "logit"),
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

fit <- glmer(rate_model ~ cluster + (1 | site) + (1 | year),
             data = modeldat_closedforests,
             weight = forest,
             family = binomial(link = "logit"),
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


# Check model -------------------------------------------------------------

# AIC
AIC(fit_null_null, fit_null, fit)

# Pseudo-R2
MuMIn::r.squaredGLMM(fit, null = fit_null)

# Residual checks
plot(DHARMa::simulateResiduals(fittedModel = fit, n = 250))

# Get model parameters ----------------------------------------------------

get_parameter <- function (mod) {
  
  # Fixed effects
  mu <- fixef(mod)
  
  Sigma <- vcov(mod)
  
  # Random intercept among years
  varcorr <- VarCorr(mod)
  sigma_year <- attr(varcorr$year, "stddev")
  
  # Random intercept + slope among sites
  sigma_sites <- attr(varcorr$site, "stddev")
  
  # Return all parameters
  return(list(mu = mu,
              Sigma = Sigma, 
              sigma_year = as.double(sigma_year), 
              sigma_sites = sigma_sites))
}

parameters <- get_parameter(fit)

write_csv(parameters$mu %>% data.frame(mu = .) %>% 
            rownames_to_column(var = "predictor"), paste0("results/closedforests/disturbance_rate_model_model_parameters_mu_closedforests.csv"))
write.csv(parameters$Sigma %>% as.matrix(.) %>% 
            as.data.frame(), paste0("results/closedforests/disturbance_rate_model_parameters_Sigma_closedforests.csv"))
write_csv(parameters$sigma_year %>% 
            data.frame(sigma_year = .), paste0("results/closedforests/disturbance_rate_model_parameters_sigma_year_closedforests.csv"))
write.csv(parameters$sigma_sites %>% 
            data.frame(sigma_site = .), paste0("results/closedforests/disturbance_rate_model_parameters_sigma_site_closedforests.csv"))

# Create model predictions

predict_from_model <- function(params = parameters, mod = fit, cluster) {
  
  names <- names(mod@frame)
  draw_fixef_cluster <- MASS::mvrnorm(n = 1, params$mu, params$Sigma)
  draw_ranef_site <- MASS::mvrnorm(n = 1, 0, params$sigma_sites)
  draw_ranef_year <- MASS::mvrnorm(n = 1, 0, params$sigma_year)
  
  preddat <- data.frame(cluster = factor(cluster, levels = c(1, 2, 3)))
  
  pred <- model.matrix(~ 1 + cluster, data = preddat) %*% fixef(mod) + draw_ranef_site + draw_ranef_year
  
  return(inv_logit(pred))
  
}

nsim <- 10000

sim <- list(data.frame(pred = replicate(nsim, predict_from_model(cluster = 1))),
            data.frame(pred = replicate(nsim, predict_from_model(cluster = 2))),
            data.frame(pred = replicate(nsim, predict_from_model(cluster = 3)))) %>%
  set_names(1:3) %>%
  bind_rows(.id = "cluster")

pred_summary_closedforest <- sim %>%
  group_by(cluster) %>%
  summarize(mean = mean(pred),
            sd = sd(pred),
            se = sd(pred) / sqrt(nsim),
            min = min(pred),
            max = max(pred),
            q0.010 = quantile(pred, 0.01),
            q0.025 = quantile(pred, 0.025),
            q0.250 = quantile(pred, 0.25),
            q0.500 = quantile(pred, 0.5),
            q0.750 = quantile(pred, 0.75),
            q0.975 = quantile(pred, 0.975),
            q0.990 = quantile(pred, 0.99))

write_csv(pred_summary_closedforest, paste0("results/closedforests/disturbance_rates_summary_clopsedforest.csv"))

p <- ggplot(sim) + 
  geom_histogram(aes(x = pred, y = ..density..), bins = 100, alpha = 0.85, fill = "grey") +
  labs(x = "Disturbance probability\n(closed forests)", y = "Count") +
  theme_minimal() +
  facet_wrap(~cluster, scales = "free") +
  geom_vline(data = pred_summary_closedforest, aes(xintercept = mean), linetype = "dashed") +
  geom_text(data = pred_summary_closedforest, aes(x = mean * 1, 
                                     y = c(3000, 1000, 100), 
                                     label = paste0(round(mean * 100, 2), " % per yr.")),
            hjust = -0.15)

ggsave("results/closedforests/disturbance_rates_closedforests.pdf", p, width = 7.5, height = 2.5)

