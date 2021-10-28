## Fit MARSS state space and Bayesian DFA models to mean age-at-maturity data
# MARSS component is used to identify ecological groupings of Chinook salmon 
# with common demographic trends, while DFA used to estimate group-specific 
# trends, generate predictions, and evaluate evidence for regime shifts using 
# hidden Markov models
# All data originates with PST Chinook technical committee and was provided by
# C. Parken


library(MARSS)
library(tidyverse)

source(here::here("R", "functions", "data_cleaning_functions.R"))

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)


gen_raw <- readRDS(here::here("data", "salmon_data", 
                              "cwt_indicator_surv_clean.RDS")) 

gen <- gen_raw %>% 
  filter(!is.na(gen_length)) %>% 
  group_by(stock) %>% 
  mutate(gen_z = as.numeric(scale(gen_length)),
         gen_cent = as.numeric(scale(gen_length, center = TRUE, 
                                      scale = FALSE))) %>% 
  ungroup() %>% 
  droplevels()

map(colnames(gen)[which(colnames(gen) %in% c("smolt", "run") |
                          str_detect(colnames(gen), "group"))], 
    function(x) {
      gen %>%
        select(stock_name, .data[[x]]) %>%
        distinct() %>%
        group_by(.data[[x]]) %>%
        tally()
    })

# dataframe of only stocks and adult groupings
stk_tbl <- gen %>% 
  group_by(stock) %>% 
  mutate(max_age = ceiling(max(gen_length)),
         max_ocean_age = ifelse(smolt == "oceantype", max_age - 1, 
                                max_age - 2)) %>% 
  ungroup() %>% 
  select(stock, stock_name, max_age, max_ocean_age, smolt, run, 
         j_group1:j_group4, a_group1:a_group4, j_group1b, j_group2b, j_group3b,
         j_group4b) %>% 
  distinct() 


## EXPLORATORY -----------------------------------------------------------------
 
#plot raw generation length data
raw_gen <- gen %>%
  ggplot(.) +
  geom_point(aes(x = year, y = gen_length, fill = a_group2), shape = 21) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group2))) +
  theme(legend.position = "top") +
  labs(y = "Mean Generation Length") +
  ggsidekick::theme_sleek()

pdf(here::here("figs", "raw_gen_trends.pdf"), width = 10, height = 7)
raw_gen
dev.off()

#distribution of generation length
gen %>%
  group_by(stock) %>%
  ggplot(.) +
  geom_histogram(aes(x = gen_length, fill = a_group4)) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group4))) +
  theme(legend.position = "top") +
  ggsidekick::theme_sleek()

gen %>%
  ggplot(.) +
  geom_point(aes(x = year, y = gen_z, fill = a_group2), shape = 21) +
  geom_point(aes(x = year, y = gen_cent, fill = a_group2), shape = 24) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group2))) +
  theme(legend.position = "top") +
  labs(y = "Mean Generation Length") +
  ggsidekick::theme_sleek()


## MARSS MODEL RUNS ------------------------------------------------------------

# make matrix
gen_mat1 <- gen %>% 
  select(year, stock, gen_length) %>% 
  pivot_wider(names_from = stock, values_from = gen_length) %>% 
  arrange(year) %>% 
  as.matrix() %>% 
  t()
gen_mat <- gen_mat1[2:nrow(gen_mat1), ]
colnames(gen_mat) <- seq(min(gen$year), max(gen$year), by = 1)

n_ts <- nrow(gen_mat)
tt <- ncol(gen_mat)


## Generic MARSS approach (model selection)
# specify the z models based on different groupings (smolt, run, a dist, j dist)
z_model_inputs <- colnames(gen)[which(colnames(gen) %in% c("smolt", "run") |
                                        str_detect(colnames(gen), "group"))]
z_models <- map(z_model_inputs, function (x) {
  stk_tbl %>% 
    pull(.data[[x]]) %>% 
    factor()
})
names(z_models) <- z_model_inputs


# additional MARSS model parameters are fixed
U <- "unequal"
R <- "diagonal and equal" 
A <- "scaling"
B <- "identity"
x0 <- "unequal"
V0 <- "zero"
Q <- 'unconstrained'
model_constants <- list(U = U, B = B, x0 = x0, A = A, V0 = V0, Q = Q)


# function to fit models
fit_marss <- function(z_name, z_in) {
  fit_model <- c(list(Z = z_in), model_constants)
  fit <- MARSS(gen_mat, model = fit_model,
               silent = FALSE, control = list(minit = 100, maxit = 500))
  #use BFGS except when equalvarcov (can't fit)
  if (fit$convergence != 0){
    fit <- MARSS(gen_mat, 
                 model = fit_model, 
                 control = list(maxit=4000, trace=1),
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

age_marss_aic_tab <- purrr::map(marss_list, "out") %>% 
  bind_rows() %>% 
  arrange(AICc) %>% 
  mutate(deltaAICc = AICc - min(AICc),
         rel_like = exp(-1 * deltaAICc / 2),
         aic_weight = rel_like / sum(rel_like))


saveRDS(age_marss_aic_tab, here::here("data", "generation_fits",
                                  "marss_aic_tab_diag_equal.RDS"))
print(readRDS(here::here("data", "generation_fits", 
                         "marss_aic_tab_diag_equal.RDS")))



## FIT BAYESIAN DFA ------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# only retain stock groups that contain more than two stocks
kept_grps <- stk_tbl %>%
  group_by(j_group3b) %>%
  tally() %>%
  filter(n > 2)

#generate tbl by group
gen_tbl <- tibble(group = levels(gen$j_group3b)) %>% 
  mutate(
    gen_mat = gen %>% 
      filter(!is.na(gen_length)) %>% 
      group_split(j_group3b) %>% 
      map(., make_mat, resp = "gen_length")
  ) %>% 
  filter(group %in% kept_grps$j_group3b) %>% 
  mutate(group = fct_relevel(as.factor(group),"sog_oceantype", after = 1)) %>% 
  arrange(group)
gen_tbl$names <- map(gen_tbl$gen_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., gen %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
gen_tbl$years <- map(gen_tbl$gen_mat, function (x) {
  as.numeric(colnames(x))
})



#fit 
furrr::future_map2(
  gen_tbl$gen_mat[c(2,3)],
  gen_tbl$group[2:3],
  .f = function (y, group) {
    
    # convergence issues when estimating phi for this stock group
    if (group == "puget_oceantype") {
      est_ar <- FALSE
    } else {
      est_ar <- TRUE
    }
    
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE,
      estimate_trend_ar = est_ar, 
      iter = 6000, warmup = 2000, chains = 1, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    
    f_name <- paste(group, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_")
    saveRDS(fit, here::here("data", "generation_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)
 
# fit2 <- fit_dfa(
#   y = gen_tbl$gen_mat[[3]], num_trends = 2, zscore = FALSE,
#   # estimate_trend_ar = TRUE, 
#   iter = 4000, chains = 4, thin = 1,
#   control = list(adapt_delta = 0.99, max_treedepth = 20)
# )
# 
# f_name <- paste(gen_tbl$group[3], "two-trend", "ar", "bayesdfa_c.RDS", sep = "_")
# saveRDS(fit2, here::here("data", "generation_fits", f_name))

# read outputs
dfa_fits <- map(gen_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "generation_fits", f_name))
})


##  CHECK DFA DIAGNOSTICS ------------------------------------------------------

# histograms of effective sample size ratio
map2(dfa_fits, gen_tbl$group, function (x, y) {
  data.frame(neff_ratio = bayesplot::neff_ratio(x$model),
             group = y)
}) %>% 
  bind_rows()  %>% 
  ggplot(.) +
  geom_histogram(aes(x = neff_ratio)) +
  facet_wrap(~group)

# check for parameters outside n_eff and rhat threshold
map(dfa_fits, function (x) {
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
  list(post, np, gen_tbl$group),
  .f = function(x, y, z) {
    all_pars <- dimnames(x)[3] %>% unlist() %>% as.character()
    to_match <- c("Z", "phi", "xstar", "sigma")
    pars <- all_pars[grepl(paste(to_match, collapse = "|"), all_pars)]
    bayesplot::mcmc_trace(x, pars = pars, np = y) +
      labs(title = z)
  }
)

# pdf(here::here("figs", "hidden", "dfa", "mean_age", "trace_plots_ar.pdf"))
# trace_list
# dev.off()


##  SAVE DFA OUTPUTS -----------------------------------------------------------

# rotate trends and add to gen_tbl (keep DFA separate because the posterior 
# estimates require a large amount of memory)
gen_tbl$rot_gen <- map(dfa_fits, rotate_trends)

# apply hidden markov model for regimes to rotated trends
hmm_list_g <- furrr::future_map(gen_tbl$rot_gen, regime_f)
map(hmm_list_g, function(x) {
  map(x, ~.$table)
} )

map(hmm_list_g, function(x) {
  map(x, function (y) plot_regime_model(y$best_model))
} )

gen_tbl$regime_trend1 <- map(hmm_list_g, function(x) x[[1]]$best_model)
gen_tbl$regime_trend2 <- map(hmm_list_g, function(x) x[[2]]$best_model)

saveRDS(gen_tbl, here::here("data", "generation_fits", "gen_tbl.RDS"))

