## Fit MARSS then Bayes DFA models to mean generation time data
# Jan 8 2021
# Update on surv_bayesDFA; instead of deciding groupings a priori, first perform
# model selection for groupings with MARSS (faster to converge than bayesDFA), 
# then fit group specific models with bayesdfa
# Updated Feb 23 with groupings edited by AVE

library(MARSS)
library(tidyverse)

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)


gen_raw <- readRDS(here::here("data", "salmonData", 
                              "cwt_indicator_surv_clean.RDS")) 

#remove stocks with no gen_length data and rare life history types (n_stk < 2)
gen <- gen_raw %>% 
  filter(!is.na(gen_length),
         #!j_group3 %in% c("col_streamtype", "north_oceantype", 
                          # "sog_streamtype"),
         #!a_group3 == "north_streamtype"
         ) %>% 
  mutate(gen_z = as.numeric(scale(gen_length))) %>% 
  droplevels()

# map(colnames(gen)[c(13:20, 22:33)], function(x) {
#   gen %>% 
#     select(stock_name, .data[[x]]) %>% 
#     distinct() %>% 
#     group_by(.data[[x]]) %>% 
#     tally()
# })

# dataframe of only stocks and adult groupings
stk_tbl <- gen %>% 
  group_by(stock) %>% 
  mutate(max_age = ceiling(max(gen_length)),
         max_ocean_age = ifelse(smolt == "oceantype", max_age - 1, 
                                max_age - 2)) %>% 
  ungroup() %>% 
  select(stock, stock_name, max_age, max_ocean_age, smolt, run, 
         j_group1:j_group4, a_group1:a_group4, j_group4b:a_group1c) %>% 
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
  geom_histogram(aes(x = log(gen_length), fill = a_group)) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group))) +
  theme(legend.position = "top") +
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

q_models <- c(#"diagonal and equal", 
              "diagonal and unequal", 
              #"equalvarcov",
              "unconstrained")

U <- "unequal"
R <- "diagonal and equal"
A <- "scaling"
B <- "identity"
x0 <- "unequal"
V0 <- "zero"
model_constants <- list(U = U, R = R, A = A, B = B, x0 = x0, V0 = V0)

# function to fit models
fit_marss <- function(z_name, z_in, q_in) {
  fit_model <- c(list(Z = z_in, Q = q_in), model_constants)
  fit <- MARSS(gen_mat, model = fit_model,
               silent = FALSE, control = list(minit = 100, maxit = 500))
  #use BFGS except when equalvarcov (can't fit)
  if (fit$convergence != 0 & q_in != "equalvarcov"){
    fit <- MARSS(gen_mat, 
                 model = fit_model, 
                 control = list(maxit=4000, trace=1),
                 inits = as.matrix(coef(fit)[1]),
                 method = "BFGS")
  }
  out <- data.frame(
    H = z_name, Q = q_in, U = U,
    logLik = fit$logLik, AICc = fit$AICc, num.param = fit$num.params,
    m = length(unique(z_in)),
    num.iter = fit$numIter, 
    converged = !fit$convergence,
    stringsAsFactors = FALSE)
  list(fit = fit, out = out)
}

# tibble containing model combinations
mod_names = expand.grid(q = q_models, 
                        z = names(z_models)) %>% 
  mutate(name = paste(z, q, sep = "-"))
mod_tbl <- tibble(
  mod_name = mod_names$name,
  z_name = mod_names$z,
  q_name = mod_names$q,
  z_models = rep(z_models, each = length(unique(q_models))),
  q_models = rep(q_models, length.out = length(mod_names$name))
)

# fit generic MARSS models
marss_list <- furrr::future_pmap(list(z_name = mod_tbl$z_name,
                                      z_in = mod_tbl$z_models,
                                      q_in = mod_tbl$q_models),
                                 .f = fit_marss,
                                 .progress = TRUE)

marss_aic_tab <- purrr::map(marss_list, "out") %>% 
  bind_rows() %>% 
  arrange(AICc) %>% 
  mutate(deltaAICc = AICc - min(AICc),
         rel_like = exp(-1 * deltaAICc / 2),
         aic_weight = rel_like / sum(rel_like))

saveRDS(marss_aic_tab, here::here("data", "generation_fits",
                                  "marss_aic_tab.RDS"))

marss_aic_tab <- readRDS(here::here("data", "generation_fits",
                                  "marss_aic_tab.RDS"))



## Bayesian DFA ----------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

a_palette <- disco::disco("muted", n = length(unique(gen$j_region)))
names(a_palette) <- unique(gen$j_region)

#helper function to spread and label input matrices for bayesdfa
make_mat <- function(x) {
  mat1 <- x %>%
    select(year, stock, gen_length) %>%
    spread(key = stock, value = gen_length) %>%
    as.matrix() 
  out_mat <- t(mat1[, 2:ncol(mat1)])
  colnames(out_mat) <- mat1[, "year"]
  return(out_mat)
}

# number of stocks per group
# kept_grps <- stk_tbl %>% 
#   group_by(a_group2) %>% 
#   tally() %>% 
#   filter(n > 2)

#generate tbl by group
gen_tbl <- tibble(group = levels(gen$a_group2)) %>% 
  mutate(
    gen_mat = gen %>% 
      filter(!is.na(gen_length)) %>% 
      group_split(a_group2) %>% 
      map(., make_mat)
  ) 
gen_tbl$names <- map(gen_tbl$gen_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., gen %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
gen_tbl$years <- map(gen_tbl$gen_mat, function (x) {
  as.numeric(colnames(x))
})

# specify one trend if there are less than 4 time series, otherwise 2
n_trend_list <- ifelse(unlist(map(gen_tbl$gen_mat, nrow)) < 4, 1, 2) 

dfa_fits <- furrr::future_map2(
  gen_tbl$gen_mat,
  n_trend_list,
  .f = function (y, n_trend) {
    fit_dfa(y = y, num_trends = n_trend, zscore = TRUE,
            iter = 4000, chains = 4, thin = 1,
            control = list(adapt_delta = 0.9, max_treedepth = 15))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)

# # save outputs
map2(dfa_fits, gen_tbl$group, function(x, y) {
  f_name <- paste(y, "bayesdfa.RDS", sep = "_")
  saveRDS(x, here::here("data", "generation_fits", f_name))
})

# read outputs
dfa_fits2 <- map(gen_tbl$group, function(y) {
  f_name <- paste(y, "bayesdfa.RDS", sep = "_") 
  readRDS(here::here("data", "generation_fits", f_name))
})
rot_list <- map(dfa_fits2, rotate_trends)

# make plots 
source(here::here("R", "functions", "plot_fitted_bayes.R"))
fit_list <- pmap(list(dfa_fits2, gen_tbl$names, gen_tbl$years), 
                 function(x, names, obs_years) {
                   plot_fitted_bayes(x, names = names$stock_name, 
                                     years = obs_years) +
                     scale_x_continuous(breaks = c(1970, 1990, 2010), 
                                        limits = c(min(gen$year), 
                                                   max(gen$year)))
                 })
trend_list <- pmap(list(rot_list, gen_tbl$years), 
                   function(x, obs_years) {
                     plot_trends(x, years = obs_years) +
                       scale_x_continuous(breaks = c(1970, 1990, 2010), 
                                          limits = c(min(gen$year), 
                                                     max(gen$year)))
                   }
)
loadings_list <- pmap(list(rot_list, gen_tbl$names), 
                   function(x, names) {
                     plot_loadings(x, names = names$stock_name) +
                       lims(y = c(-2, 2))
                   }
) 


fig_path <- paste("figs", "dfa", "bayes", "generation_length", sep = "/")

pdf(here::here(fig_path, "fits.pdf"))
fit_list
dev.off()

pdf(here::here(fig_path, "trends.pdf"))
trend_list
dev.off()

pdf(here::here(fig_path, "loadings.pdf"))
loadings_list
dev.off()


# Make 4 panel plot of dominant trends by adult grouping
rot_dat <- pmap(list(rot_list, gen_tbl$group, gen_tbl$years), 
                function(x, y, z) {
                  data.frame(
                    mean = x$trends_mean[1, ],
                    lo = x$trends_lower[1, ],
                    hi = x$trends_upper[1, ],
                    group = y,
                    years = z
                  )
                }) %>% 
  bind_rows() %>% 
  glimpse()

trends <- rot_dat %>% 
  ggplot(., aes(x = years, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  facet_wrap(~group) +
  ggsidekick::theme_sleek() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Brood Year", y = "Shared Trend") +
  theme(legend.position = "top")

png(here::here("figs", "dfa", "bayes", "generation_length", "trend_all_groups.png"), 
    height = 5.5, width = 7.5, res = 300, units = "in")
trends
dev.off()
