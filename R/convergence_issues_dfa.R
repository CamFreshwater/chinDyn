## Non-convergence Example
# DFA model that throws non-convergence errors when data are centered/scaled, but
# not when centered.
# Input is time series of average age at maturity for stocks of Pacific salmon 
# returning to Puget Sound. Similar datasets (i.e. mean age of stocks returning
# to other regions) converge well.
# NOTE: I'm now taking your advice and fitting a model of maximum complexity
# that estimates nu and AR/MA parameters. No convergence issues with centered
# data there either.

library(MARSS)
library(tidyverse)
library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


dat <- readRDS(here::here("data", "salmonData", 
                              "cwt_indicator_surv_clean.RDS")) %>%  
  filter(j_group3b == "puget_oceantype") %>% 
  droplevels()

mat1 <- dat %>%
  select(year, stock, M) %>%
  spread(key = stock, value = M) %>%
  as.matrix() 
mat_in <- t(mat1[, 2:ncol(mat1)])
colnames(mat_in) <- mat1[, "year"]


fit_scale <- fit_dfa(
  y = mat_in, num_trends = 2, zscore = TRUE,
  iter = 3000, chains = 4, thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)
fit_center <- fit_dfa(
  y = mat_in, num_trends = 2, zscore = FALSE,
  iter = 3000, chains = 4, thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)
