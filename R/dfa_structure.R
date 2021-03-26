## BayesDFA Model Structure Evaluation
# Evaluate different structures for Bayesian DFA. Specifically whether nu, phi
# and psi parameters can be estimated if a) observation errors are assumed to 
# vary (i.e. R is diagonal and unequal) and b) observation errors covary (i.e.
# R is unconstrained)
# March 23, 2021

library(tidyverse)
library(bayesdfa)
library(rstan)
library(bayesplot)

ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# import tbls made in m_dfa and gen_dfa
surv_tbl <- readRDS(here::here("data", "mortality_fits", "surv_tbl.RDS"))
gen_tbl <- readRDS(here::here("data", "generation_fits", "gen_tbl.RDS"))


# focus on one stock group and add additional parameters
tbl_in <- rbind(surv_tbl %>% mutate(var = "surv") %>% rename(mat_in = m_mat),
                gen_tbl %>% mutate(var = "age") %>% rename(mat_in = gen_mat)) %>% 
  filter(group == "sog_oceantype")

tbl_in2 <- expand.grid(var = c("age", "surv"),
            saturated = c(TRUE, FALSE),
            r_mat = c("equal", "unequal")
            )  %>% 
   left_join(tbl_in, ., by = "var") %>% 
  filter(var == "age", 
         r_mat == "unequal")

# fit model
furrr::future_pmap(
  list(tbl_in2$mat_in, tbl_in2$var, tbl_in2$r_mat, tbl_in2$saturated),
  .f = function (y, var, r_mat, saturated) {
    varIndx <- if (r_mat %in% c("unequal", "unconstrained")) {
       seq(1, nrow(y), by = 1)
    } else  NULL
    
    sat_name <- ifelse(saturated == TRUE, "sat", "unsat")
    f_name <- paste(var, sat_name, r_mat, "bayesdfa.RDS", sep = "_")
    
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE, 
      varIndx = varIndx, 
      #estimate_nu = saturated, 
      estimate_trend_ar = saturated, estimate_trend_ma = saturated,
      iter = 3000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    saveRDS(fit, here::here("data", "model_structure_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)


files <- list.files(path = here::here("data", "model_structure_fits"), 
                    pattern = "\\.RDS$", full.names = TRUE)
fit_list <- lapply(files, readRDS)


## diagnostics check
# check number of samples and Rhat
map(fit_list, function (x) {
  as.data.frame(summary(x$model)$summary) %>% 
    filter(n_eff < 250 | Rhat > 1.05)
})

fit_tbl <- tbl_in2 %>% 
  arrange(var, -saturated) %>% 
  select(var, saturated, r_mat)
fit_tbl$fits <- fit_list

# check for convergence
fit_tbl$converged <- map(fit_tbl$fits, is_converged) %>% unlist()

# sum div transitions
div_trans <- map(fit_tbl$fits, function (y) {
  par_list <- get_sampler_params(y$model, inc_warmup = FALSE)
  map(par_list, function (x) {
    sum(x[ , "divergent__"]) 
  }) %>% 
    unlist() %>% 
    sum()
})
fit_tbl$divergent <- div_trans %>% unlist()

# plot divergent transitions focusing on candidate parameters
all_pars <- dimnames(fit_tbl$post[[1]])[3] %>% unlist() %>% as.character()
to_match <- c("theta", "phi", "nu", "sigma")
pars <- all_pars[grepl(paste(to_match, collapse = "|"), all_pars)]

fit_tbl$post <- map(fit_tbl$fits, function (x) {
  as.array(x$samples)
})
fit_tbl$np <- map(fit_tbl$fits, function (x) {
  nuts_params(x$model)
})

mcmc_trace(fit_tbl$post[[1]], pars = pars, np = fit_tbl$np[[1]])

age_mat <- gen_tbl$gen_mat[[4]]
r_index <- seq(1, nrow(age_mat), by = 1)

fit1 <- fit_dfa(
  y = age_mat, num_trends = 2, zscore = FALSE, 
  varIndx = r_index, 
  estimate_trend_ar = TRUE, estimate_trend_ma = TRUE,
  iter = 3000, chains = 4, thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)
fit2 <- fit_dfa(
  y = age_mat, num_trends = 2, zscore = FALSE, 
  varIndx = r_index, 
  estimate_trend_ar = FALSE, estimate_trend_ma = FALSE,
  iter = 3000, chains = 4, thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)