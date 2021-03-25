## Non-convergence Example
# DFA model that throws non-convergence errors when unique observation variance 
# is estimated for each time series. Otherwise ok. 
# Input is time series of average age at maturity for stocks of Pacific salmon 
# returning to Puget Sound. Some other similar datasets (i.e. mean age of stocks
# returning to other regions) converge with this model, others don't.

library(tidyverse)
library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


dat <- readRDS(here::here("data", "salmonData", 
                              "cwt_indicator_surv_clean.RDS")) %>%  
  filter(j_group3b == "sog_oceantype", 
         !is.na(gen_length)) %>% 
  droplevels()

mat1 <- dat %>%
  select(year, stock, gen_length) %>%
  spread(key = stock, value = gen_length) %>%
  as.matrix() 
mat_in <- t(mat1[, 2:ncol(mat1)])
colnames(mat_in) <- mat1[, "year"]

var_seq <- seq(1, nrow(mat_in), by = 1)

# fit model with diagonal and unequal R matrix
fit_unequal <- fit_dfa(
  y = mat_in, num_trends = 2, zscore = FALSE,
  varIndx = var_seq,
  iter = 3000, chains = 4, thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

# fit model with diagonal and equal R matrix
fit_diag <- fit_dfa(
  y = mat_in, num_trends = 2, zscore = FALSE,
  iter = 3000, chains = 4, thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)
