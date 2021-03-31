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
surv_tbl <- readRDS(here::here("data", "mortality_fits", "surv_tbl.RDS")) %>% 
  select(group:years)
gen_tbl <- readRDS(here::here("data", "generation_fits", "gen_tbl.RDS")) %>% 
  select(group:years)


# focus on one stock group and add additional parameters
tbl_in <- rbind(surv_tbl %>% mutate(var = "surv") %>% rename(mat_in = m_mat),
                gen_tbl %>% mutate(var = "age") %>% rename(mat_in = gen_mat)) %>% 
  filter(group == "sog_oceantype")

tbl_in2 <- expand.grid(var = c("age", "surv"),
            r_mat = c("equal", "unequal"),
            nu = c(0, 1),
            ar = c(0, 1),
            ma = c(0, 1)
            )  %>% 
   left_join(tbl_in, ., by = "var") %>% 
  filter(var == "surv", 
         r_mat == "equal",
         !(nu == 1 & ar == 1 & ma == 1))

# fit model
furrr::future_pmap(
  list(tbl_in2$mat_in, tbl_in2$var, tbl_in2$r_mat, tbl_in2$nu,
       tbl_in2$ar, tbl_in2$ma),
  .f = function (y, var_in, r_mat, nu, ar, ma) {
    varIndx <- if (r_mat %in% c("unequal", "unconstrained")) {
       seq(1, nrow(y), by = 1)
    } else  NULL
    
    sat_name1 <- if (nu) "nu_" 
    sat_name2 <- if (ar) "ar_" 
    sat_name3 <- if (ma) "ma_" 
    sat_name <- paste(sat_name1, sat_name2, sat_name3, sep = "")
    
    f_name <- paste(var_in, sat_name, r_mat, "bayesdfa.RDS", sep = "_")

    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE, 
      varIndx = varIndx, 
      estimate_nu = nu, estimate_trend_ar = ar, estimate_trend_ma = ma,
      iter = 2000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    saveRDS(fit, here::here("data", "model_structure_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)


files <- list.files(path = here::here("data", "model_structure_fits"), 
                    pattern = "\\.RDS$", full.names = TRUE) 
subfiles <- files[grepl("surv", files)]
fit_list <- lapply(subfiles, readRDS)


## diagnostics check
# check number of samples and Rhat
map(fit_list, function (x) {
  as.data.frame(summary(x$model)$summary) %>% 
    filter(n_eff < 250 | Rhat > 1.05)
})

# fit_tbl <- tbl_in2 %>%
#   arrange(var, nu, ma) %>%
#   select(var, saturated, r_mat)
# fit_tbl$fits <- fit_list

# check for convergence
fit_tbl$converged <- map(fit_list, is_converged) %>% unlist()

# sum div transitions
div_trans <- map(fit_list, function (y) {
  par_list <- get_sampler_params(y$model, inc_warmup = FALSE)
  map(par_list, function (x) {
    sum(x[ , "divergent__"]) 
  }) %>% 
    unlist() %>% 
    sum()
})

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

# loo_tbl <- tibble(group = rep(gen_tbl$group, 2),
#                   m = rep(c(2, 1), each = length(gen_tbl$group)),
#                   fits = c(dfa_fits, dfa_fits1))
loo <- map(fit_list, bayesdfa::loo)
looic <- map(loo, function(x) x$estimates["looic", "Estimate"]) %>%
  unlist()


loo_tbl_out <- loo_tbl %>%
  group_by(group) %>%
  mutate(min_looic = min(looic),
            delta_loo = looic - min_looic)
