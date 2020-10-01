## Survival and covariate gams
# Quick mixed effects model analysis of correlations between seal density 
# and Chinook survival to age 2 
# For now focused only on Salish Sea individuals
# Oct. 1 2020

library(mgcv)
library(tidyverse)

cov <- read.csv(here("data/salmonData/survCovariateAnom.csv"), 
                stringsAsFactors = F)
surv <- read.csv(here::here("data/salmonData/cwt_indicator_surv_clean.csv"), 
                 stringsAsFactors = FALSE) %>% 
  filter(agg_reg == "SS",
         !is.na(M)) %>% 
  select(year, stock, smolt, survival, M) 

# raw plot 
plot_dat <- surv %>% 
  left_join(., cov, by = "year")

plot_dat %>% 
  filter(smolt == "streamtype") %>% 
  ggplot(., aes(x = anomaly, y = M)) +
  geom_point() +
  geom_smooth(method = "gam") + 
  facet_wrap(~metric) 
  

# Test fit ---------------------------------------------------------------------

# Fit various GAM structures
temp <- plot_dat %>% 
  filter(metric == "seal_anom") %>% 
  mutate(smolt = as.factor(smolt),
         stock = as.factor(stock))
m1 <- gam(M ~ s(anomaly), data = temp)
m2 <- gam(M ~ s(anomaly, by = smolt), data = temp)
m3 <- gam(M ~ s(anomaly, m = 2) + s(anomaly, by = smolt, m = 1), data = temp)
m4 <- gam(M ~ s(anomaly, m = 2) + s(anomaly, by = smolt, m = 1) + 
            s(stock, bs = "re"), 
          data = temp)
AIC(m1, m2, m3, m4)

# Compare predictions
seal_seq <- seq(-1.5, 1, length.out = 200)
new_dat <- expand.grid(anomaly = seal_seq, smolt = unique(temp$smolt))
pred_dat <- predict(m3, new_dat, se.fit = TRUE)
new_dat$fit <- as.numeric(pred_dat$fit)
new_dat$up <- as.numeric(new_dat$fit + (1.96 * pred_dat$se.fit))
new_dat$low <- as.numeric(new_dat$fit - (1.96 * pred_dat$se.fit))

ggplot(new_dat, aes(x = anomaly)) +
  geom_line(aes(y = fit, color = smolt)) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = smolt, alpha = 0.3)) 

# merge covariate data
dat <- cov %>% 
  select(year, metric, env_anomaly = anomaly) %>% 
  nest(data = c(year, env_anomaly)) %>% 
  mutate(
    #for each covariate combine with survival data
    data = map(data, function (x) {
      surv %>% 
        left_join(., x, by = "year")
    }),
    raw_plot = map2(metric, data, function (name, data) {
      ggplot(data) +
        geom_point(aes(x = M, y = env_anomaly))
    })
  )

# plot data
dat$raw_data <- map2(dat)