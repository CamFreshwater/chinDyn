## Survival and seal density models
# Quick mixed effects model analysis of correlations between seal density 
# and Chinook survival to age 2 
# July 10, 2020

# ------------------------------------------------------------------------------

library(lme4)
library(brms)
library(rethinking)
library(tidyr)
library(modelr)
library(tidybayes)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

seals <- read.csv(here::here("data/salmonData/sealPopEst.csv"), 
                  stringsAsFactors = FALSE) %>% 
  rename(year = "Estimate", mean = "Mean", low = "X2.5th", up = "X97.5th", 
         reg = "Region") %>% 
  mutate(mean = mean / 1000) %>% 
  filter(reg == "SOG")

# eyDat <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
#                   stringsAsFactors = FALSE) 
surv <- read.csv(here::here("data/salmonData/cwt_indicator_surv_clean.csv"), 
                 stringsAsFactors = FALSE)


# CLEAN DATA -------------------------------------------------------------------

ss <- surv %>% 
  filter(agg_reg == "SS") %>% 
  left_join(., 
            seals %>% 
              select(year, mean_seal = mean),
            by = "year")

# plots 
avg <- ss %>% 
  group_by(year) %>% 
  summarize(mean_M = mean(M, na.rm = T),
            mean_seal = mean(mean_seal)) 

ggplot(avg) +
  geom_point(aes(x = mean_seal, y = mean_M))

ggplot(ss) +
  geom_point(aes(x = mean_seal, y = M)) + 
  facet_wrap(~stock)


# FIT MODELS -------------------------------------------------------------------

# models with varying intercepts or slopes
ss1 <- brm(M ~ mean_seal + (1 | stock), data = ss,
           warmup = 500, iter = 1500, 
           cores = 4, chains = 4, 
           seed = 123)
ss2 <- brm(M ~ mean_seal + (1 + mean_seal | stock), data = ss,
           warmup = 500, iter = 1500, 
           cores = 4, chains = 4, 
           seed = 123)

mod_list <- list(ss1, ss2)

# diagnostics
map(mod_list, function(x) {
  x %>%
    ggmcmc::ggs() %>% 
    filter(Parameter %in% c("b_Intercept", "b_mean_seal")) %>% 
    ggplot(.,
           aes(x   = Iteration,
               y   = value, 
               col = as.factor(Chain)))+
    geom_line() +
    geom_vline(xintercept = 500)+
    facet_grid(Parameter ~ . ,
               scale  = 'free_y',
               switch = 'y')
})


# sample from posterior
# stock-specific predictions
post_pred <- ss %>% 
  group_by(stock) %>% 
  data_grid(mean_seal = seq_range(mean_seal, n = 101)) %>% 
  add_fitted_draws(ss2, n = 100) %>% 
  left_join(., ss %>% select(stock, smolt) %>% distinct(),
            by = "stock") 

stock_spec <- ggplot(data = post_pred, aes(x = mean_seal, y = M, color = smolt)) +
  geom_line(aes(y = .value, group = paste(stock, .draw)), alpha = 0.1) +
  geom_point(data = ss) + 
  labs(x = "Seal Abundance") +
  facet_wrap(~stock)

# average fixed effects predictions
pred_fe <- tibble(mean_seal = seq(min(ss$mean_seal, na.rm = T), 
                                  max(ss$mean_seal, na.rm = T), 
                                  length.out = 100)) %>% 
  add_fitted_draws(ss2,
                   re_formula = NA,
                   scale = "response", n = 1e3)
pred_mean_fe <- pred_fe %>% 
  group_by(mean_seal) %>% 
  summarize(.value = mean(.value))

fixed_pred_plot <- ggplot(pred_fe, aes(x = mean_seal, y = .value)) +
  stat_interval(alpha = 0.5) +
  geom_line(data = pred_mean_fe, color = "black", lwd = 2) +
  labs(x = "Seal Abundance", y = "Predicted M",
       title = "Fixed Only Posterior Predictions")
fixed_pred_plot

# average fixed effects predictions
pred_re <- crossing(mean_seal = seq(min(ss$mean_seal, na.rm = T), 
                                    max(ss$mean_seal, na.rm = T), 
                                    length.out = 100),
                    stock = unique(ss$stock)) %>% 
  add_fitted_draws(ss2,
                   scale = "response", n = 1e3)
pred_mean_re <- pred_re %>% 
  group_by(mean_seal) %>% 
  summarize(.value = mean(.value))

mixed_pred_plot <- ggplot(pred_re, aes(x = mean_seal, y = .value)) +
  stat_interval(alpha = 0.5) +
  geom_line(data = pred_mean_re, color = "black", lwd = 2) +
  labs(x = "Seal Abundance", y = "Predicted M", 
       title = "Fixed and Random Posterior Predictions")
mixed_pred_plot

pdf(here::here("figs", "seal_mixed_model.pdf"))
stock_spec
fixed_pred_plot
mixed_pred_plot
dev.off()