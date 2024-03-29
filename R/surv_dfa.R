## Fit MARSS state space and Bayesian DFA models to logit survival data
# MARSS component is used to identify ecological groupings of Chinook salmon 
# with common demographic trends, while DFA used to estimate group-specific 
# trends, generate predictions, and evaluate evidence for regime shifts using 
# hidden Markov models
# All data originates with PST Chinook technical committee and was provided by
# C. Parken

library(MARSS)
library(tidyverse)

source(here::here("R", "functions", "data_cleaning_functions.R"))

ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)

surv <- readRDS(here::here("data/salmon_data/cwt_indicator_surv_clean.RDS")) %>%
  mutate(logit_surv = boot::logit(survival)) %>% 
  #remove stocks that are aggregates of others on CP's advice
  # TST combines STI/TAK and AKS combines SSA and NSA
  filter(!stock %in% c("TST", "AKS")) %>%
  group_by(stock) %>% 
  mutate(logit_surv_z = as.numeric(scale(logit_surv)),
         logit_surv_cent = as.numeric(scale(logit_surv, center = TRUE, 
                                            scale = FALSE))) %>%
  ungroup() %>% 
  droplevels()


# dataframe of only stocks and juvenile groupings
stk_tbl <- surv %>% 
  ungroup() %>% 
  select(stock, stock_name, smolt, run, 
         j_group1:j_group4, a_group1:a_group4, j_group1b, j_group2b, j_group3b,
         j_group4b) %>% 
  distinct() 


## EXPLORATORY -----------------------------------------------------------------

#plot transformed survival data
surv %>%
  filter(!is.na(logit_surv)) %>%
  ggplot(.) +
  geom_point(aes(x = year, y = logit_surv, fill = j_group3), shape = 21) +
  facet_wrap(~ fct_reorder(stock, as.numeric(j_group3))) +
  theme(legend.position = "top") +
  labs(y = "logit(survival rate)")

surv %>%
  filter(!is.na(M)) %>%
  ggplot(.) +
  geom_histogram(aes(M, fill = j_group3)) +
  facet_wrap(~ fct_reorder(stock, as.numeric(j_group3))) +
  theme(legend.position = "top") +
  labs(y = "logit(survival rate)")


## MARSS MODEL RUNS ------------------------------------------------------------

# make matrix of logit survival rates
surv_mat <- surv %>% 
  select(year, stock, logit_surv) %>% 
  pivot_wider(names_from = stock, values_from = logit_surv) %>% 
  arrange(year) %>% 
  as.matrix() %>% 
  t()
surv_mat <- surv_mat[2:nrow(surv_mat), ]
colnames(surv_mat) <- seq(min(surv$year), max(surv$year), by = 1)

n_ts <- nrow(surv_mat)
tt <- ncol(surv_mat)


# specify the z models based on different groupings (smolt, run, a dist, j dist)
z_model_inputs <- colnames(surv)[which(colnames(surv) %in% c("smolt", "run") |
                                         str_detect(colnames(surv), "group"))]
z_models <- map(z_model_inputs, function (x) {
  stk_tbl %>% 
    pull(.data[[x]]) %>% 
    factor()
})
names(z_models) <- z_model_inputs


U <- "unequal"
R <- "diagonal and equal"
A <- "scaling" 
B <- "identity"
x0 <- "unequal"
V0 <- "zero"
Q <- "unconstrained"
model_constants <- list(U = U, B = B, x0 = x0, V0 = V0, Q = Q)


# function to fit models
fit_marss <- function(z_name, z_in) {
  fit_model <- c(list(Z = z_in), model_constants)
  fit <- MARSS(surv_mat, model = fit_model,
               silent = FALSE, 
               control = list(minit = 100, maxit = 500, safe = TRUE))
  if (fit$convergence != 0){
    fit <- MARSS(surv_mat, 
                 model = fit_model, 
                 control = list(maxit = 4000, trace = 1),
                 inits = as.matrix(coef(fit)[1]),
                 method = "BFGS")
  }
  out <- data.frame(
    H = z_name, U = U,
    logLik = fit$logLik, AICc = fit$AICc, num.param = fit$num.params,
    m = length(unique(z_in)),
    num.iter = fit$numIter, 
    converged = !fit$convergence,
    stringsAsFactors = FALSE)
  list(fit = fit, out = out)
}


# tibble containing model combinations
mod_names = expand.grid(z = names(z_models)) %>% 
  mutate(name = paste(z, sep = "-"))
mod_tbl <- tibble(
  mod_name = mod_names$name,
  z_name = mod_names$z,
  z_models = z_models
)


# fit generic MARSS models
marss_list <- furrr::future_pmap(list(z_name = mod_tbl$z_name,
                                      z_in = mod_tbl$z_models),
                                 .f = fit_marss,
                                 .progress = TRUE)

# save AIC table
marss_aic_tab <- purrr::map(marss_list, "out") %>% 
  bind_rows() %>% 
  arrange(AICc) %>% 
  mutate(deltaAICc = AICc - min(AICc),
         rel_like = exp(-1 * deltaAICc / 2),
         aic_weight = rel_like / sum(rel_like))


saveRDS(marss_aic_tab, here::here("data", "survival_fits",
                                  "marss_aic_tab_diag_equal.RDS"))
print(readRDS(here::here("data", "survival_fits", 
                         "marss_aic_tab_diag_equal.RDS")))


## Bayesian DFA ----------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# number of stocks per group
kept_grps <- stk_tbl %>% 
  group_by(j_group3b) %>% 
  tally() %>% 
  filter(n > 2) %>% 
  droplevels()

#generate tbl by group with up to two trends
surv_tbl <- tibble(group = levels(surv$j_group3b)) %>% 
  mutate(
    surv_mat = surv %>% 
      filter(!is.na(M)) %>% 
      group_split(j_group3b) %>% 
      map2(., "logit_surv", make_mat)
  ) %>% 
  filter(group %in% kept_grps$j_group3b) %>% 
  mutate(group = fct_relevel(as.factor(group),"sog_oceantype", after = 1)) %>% 
  arrange(group)
surv_tbl$names <- map(surv_tbl$surv_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., surv %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
surv_tbl$years <- map(surv_tbl$surv_mat, function (x) {
  as.numeric(colnames(x))
})


# fit model
furrr::future_map2(
  surv_tbl$surv_mat,
  surv_tbl$group,
  .f = function (y, group) {
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE,
      estimate_trend_ar = TRUE,
      iter = 4000, warmup = 2000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    f_name <- paste(group, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_")
    saveRDS(fit, here::here("data", "survival_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)

# read outputs
dfa_fits <- map(surv_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_")
  readRDS(here::here("data", "survival_fits", f_name))
})


##  CHECK DFA DIAGNOSTICS ------------------------------------------------------


map2(dfa_fits, surv_tbl$group, function (x, y) {
  data.frame(neff_ratio = bayesplot::neff_ratio(x$model),
             group = y)
}) %>% 
  bind_rows()  %>% 
  ggplot(.) +
  geom_histogram(aes(x = neff_ratio)) +
  facet_wrap(~group)

purrr::map(dfa_fits, function (x) {
  as.data.frame(summary(x$model)$summary) %>% 
    filter(n_eff < 200 | Rhat > 1.03)
})



# chain plots for key pars
post <- map(dfa_fits, function (x) {
  as.array(x$samples)
})
np <- map(dfa_fits, function (x) {
  bayesplot::nuts_params(x$model)
})

trace_list <- pmap(
  list(post, np, surv_tbl$group),
  .f = function(x, y, z) {
    all_pars <- dimnames(x)[3] %>% unlist() %>% as.character()
    to_match <- c(#"theta", 
      "Z", "phi", "xstar", "sigma")
    pars <- all_pars[grepl(paste(to_match, collapse = "|"), all_pars)]
    bayesplot::mcmc_trace(x, pars = pars, np = y) +
      labs(title = z)
  }
)

pdf(
  here::here("figs", "supp_figs", "surv_trace_plots_ar.pdf")
  )
trace_list
dev.off()


##  SAVE DFA OUTPUTS -----------------------------------------------------------

# rotate trends and add to surv_tbl (keep DFA separate because they're huge)
surv_tbl$rot_surv <- map(dfa_fits, rotate_trends)

# test for evidence of regimes 
hmm_list_m <- furrr::future_map(surv_tbl$rot_surv, regime_f, 
                                .options = furrr::furrr_options(seed = TRUE))
map(hmm_list_m, function(x) {
  map(x, ~.$table)
} )

map(hmm_list_m, function(x) {
  map(x, function (y) plot_regime_model(y$best_model))
} )

surv_tbl$regime_trend1 <- map(hmm_list_m, function(x) x[[1]]$best_model)
surv_tbl$regime_trend2 <- map(hmm_list_m, function(x) x[[2]]$best_model)


#export
saveRDS(surv_tbl, here::here("data", "survival_fits", "surv_tbl.RDS"))

