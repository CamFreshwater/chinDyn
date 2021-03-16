## Fit MARSS then Bayes DFA models to instantaneous mortality data
# Jan 6 2021
# Update on surv_bayesDFA; instead of deciding groupings a priori, first perform
# model selection for groupings with MARSS (faster to converge than bayesDFA), 
# then fit group specific models with bayesdfa
# Updated Feb 23 2020 with Columbia split from other juveniles

library(MARSS)
library(tidyverse)

ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)

surv <- readRDS(here::here("data/salmonData/cwt_indicator_surv_clean.RDS")) %>% 
  #remove stocks that are aggregates of others on CP's advice
  # TST combines STI/TAK and AKS combines SSA and NSA
  filter(!stock %in% c("TST", "AKS"),
         year < 2017) %>%
  group_by(stock) %>% 
  mutate(M_z = as.numeric(scale(M)),
         M_cent = as.numeric(scale(M, center = TRUE, scale = FALSE)),
         gen_cent = as.numeric(scale(gen_length, center = TRUE, 
                                     scale = FALSE))) %>%
  ungroup() %>% 
  droplevels()
  

# dataframe of only stocks and juvenile groupings
stk_tbl <- surv %>% 
  ungroup() %>% 
  select(stock, stock_name, smolt, run, 
         j_group1:j_group4, a_group1:a_group4, j_group4b:j_group1b) %>% 
  distinct() 

# subset of stocks for test run
# keep_stks <- surv %>% 
#   filter(year > 1980,
#          year < 2010,
#          !is.na(M)) %>% 
#   group_by(stock) %>% 
#   tally() %>% 
#   filter(n == 29)
# 
# surv_sub <- surv %>% 
#   filter(#year > 1980,
#          #year < 2010,
#          stock %in% keep_stks$stock)
# stk_tbl_sub <- stk_tbl %>% 
#   filter(stock %in% keep_stks$stock)


## EXPLORATORY -----------------------------------------------------------------

#plot raw survival data
surv %>%
  filter(!is.na(M)) %>%
  ggplot(.) +
  geom_point(aes(x = year, y = M, fill = j_group3), shape = 21) +
  facet_wrap(~ fct_reorder(stock, as.numeric(j_group3))) +
  theme(legend.position = "top") +
  labs(y = "M")

library(mgcv)
s_mod <- gam(gen_cent ~ s(M, m = 2, bs = "tp", k = 4) + 
               s(M, by = stock, m = 1, bs = "tp", k = 4) + 
               s(stock, bs = "re"), 
             data = surv, method = "REML")

new_dat <- expand.grid(
  M = seq(min(surv$M, na.rm = T), max(surv$M, na.rm = T), n = 100))

# fixed effects predictions
preds <- predict(s_mod, new_dat, se.fit = TRUE, 
                 exclude = excl_pars[grepl("stock", excl_pars)],
                 newdata.guaranteed = TRUE)
new_dat2 <- new_dat %>% 
  mutate(link_fit = as.numeric(preds$fit),
         link_se = as.numeric(preds$se.fit),
         pred_surv = link_fit,
         pred_surv_lo = link_fit + (qnorm(0.025) * link_se),
         pred_surv_up = link_fit + (qnorm(0.975) * link_se)
  )



## MARSS MODEL RUNS ------------------------------------------------------------

# make matrix of natural mortality rates
m_mat <- surv %>% 
  select(year, stock, M) %>% 
  pivot_wider(names_from = stock, values_from = M) %>% 
  arrange(year) %>% 
  as.matrix() %>% 
  t()
m_mat <- m_mat[2:nrow(m_mat), ]
colnames(m_mat) <- seq(min(surv$year), max(surv$year), by = 1)

n_ts <- nrow(m_mat)
tt <- ncol(m_mat)


## Generic MARSS approach

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
R <- "diagonal and unequal"
A <- "scaling" 
B <- "identity"
x0 <- "unequal"
V0 <- "zero"
Q <- "unconstrained"
model_constants <- list(U = U, B = B, x0 = x0, V0 = V0, Q = Q)

# function to fit models
fit_marss <- function(z_name, z_in) {
  fit_model <- c(list(Z = z_in), model_constants)
  fit <- MARSS(m_mat, model = fit_model,
               silent = FALSE, 
               control = list(minit = 100, maxit = 500, safe = TRUE))
  if (fit$convergence != 0){
    fit <- MARSS(m_mat, 
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

marss_aic_tab <- purrr::map(marss_list, "out") %>% 
  bind_rows() %>% 
  arrange(AICc) %>% 
  mutate(deltaAICc = AICc - min(AICc),
         rel_like = exp(-1 * deltaAICc / 2),
         aic_weight = rel_like / sum(rel_like))

saveRDS(marss_aic_tab, here::here("data", "mortality_fits",
                               "marss_aic_tab.RDS"))
marss_aic_tab <- readRDS(here::here("data", "mortality_fits",
                                  "marss_aic_tab.RDS"))


## Bayesian DFA ----------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#helper function to spread and label input matrices for bayesdfa
make_mat <- function(x) {
  mat1 <- x %>%
    select(year, stock, M) %>%
    spread(key = stock, value = M) %>%
    as.matrix() 
  out_mat <- t(mat1[, 2:ncol(mat1)])
  colnames(out_mat) <- mat1[, "year"]
  return(out_mat)
}

# number of stocks per group
kept_grps <- stk_tbl %>% 
  group_by(j_group3b) %>% 
  tally() %>% 
  filter(n > 2) %>% 
  droplevels()

#generate tbl by group with up to two trends
surv_tbl <- tibble(group = levels(surv$j_group3b)) %>% 
  mutate(
    m_mat = surv %>% 
      filter(!is.na(M)) %>% 
      group_split(j_group3b) %>% 
      map(., make_mat)
  ) %>% 
  filter(group %in% kept_grps$j_group3b)
surv_tbl$names <- map(surv_tbl$m_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., surv %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
surv_tbl$years <- map(surv_tbl$m_mat, function (x) {
  as.numeric(colnames(x))
})

#export
# saveRDS(surv_tbl, here::here("data", "mortality_fits", "surv_tbl.RDS"))


furrr::future_map2(
  surv_tbl$m_mat,
  surv_tbl$group,
  .f = function (y, group) {
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE, 
      estimate_nu = TRUE, estimate_trend_ar = TRUE, estimate_trend_ma = FALSE,
      iter = 4000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    f_name <- paste(group, "two-trend", "bayesdfa_c.RDS", sep = "_")
    saveRDS(fit, here::here("data", "mortality_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)


# fit <- fit_dfa(
#   y = surv_tbl$m_mat[[2]], num_trends = 2, zscore = FALSE, 
#   estimate_nu = TRUE, estimate_trend_ar = TRUE, estimate_trend_ma = FALSE,
#   iter = 8000, chains = 4, thin = 1, warmup = 2000,
#   control = list(adapt_delta = 0.99, max_treedepth = 20)
# )
# f_name <- paste(surv_tbl[2, ]$group, "two-trend", "bayesdfa_c.RDS", sep = "_")
# saveRDS(fit, here::here("data", "mortality_fits", f_name))


# read outputs
dfa_fits <- map(surv_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "mortality_fits", f_name))
})


# check diagnostics
map2(dfa_fits, surv_tbl$group, function (x, y) {
  data.frame(neff_ratio = bayesplot::neff_ratio(x$model),
             group = y)
}) %>% 
  bind_rows()  %>% 
  ggplot(.) +
  geom_histogram(aes(x = neff_ratio)) +
  facet_wrap(~group)

map(dfa_fits, function (x) {
  as.data.frame(summary(x$model)$summary) %>% 
    filter(n_eff < 500 | Rhat > 1.05)
})


## MAKE PLOTS ------------------------------------------------------------------

source(here::here("R", "functions", "plot_fitted_bayes.R"))


# predicted fits
pred_list <- pmap(list(dfa_fits, surv_tbl$names, surv_tbl$years), 
                  fitted_preds,
                  descend_order = TRUE)
# scale colors based on observed range over entire dataset
col_ramp_in <- pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range()

fit_list <- map(pred_list, plot_fitted_pred, col_ramp = col_ramp_in, facet_col = 5)

mort_fit_panel <- cowplot::plot_grid(
  fit_list[[1]], fit_list[[2]], fit_list[[3]], fit_list[[4]], fit_list[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1 
)


rot_list <- map(dfa_fits, rotate_trends)

# make plots 
trend_list <- pmap(list(rot_list, surv_tbl$years), 
                   function(x, obs_years) {
                     plot_trends(x, years = obs_years) +
                       scale_x_continuous(breaks = c(1970, 1990, 2010), 
                                          limits = c(min(surv$year), 
                                                     max(surv$year)))
                   }
)


loadings_list <- pmap(list(rot_list, surv_tbl$names), 
                      function(x, names) {
                        plot_loadings(x, names = names$stock_name) +
                          lims(y = c(-2, 2))
                      }
) 

fig_path <- paste("figs", "dfa", "bayes", "mortality", sep = "/")

pdf(here::here(fig_path, "fits.pdf"))
fit_list
dev.off()

pdf(here::here(fig_path, "trends.pdf"))
trend_list
dev.off()

pdf(here::here(fig_path, "loadings.pdf"))
loadings_list
dev.off()


# Make 4 panel plot of dominant trends by juvenile grouping
rot_dat <- pmap(list(rot_list, kept_grps$j_group3, surv_tbl$years), 
     function(x, y, z) {
  data.frame(
    mean = c(x$trends_mean[1, ], x$trends_mean[2, ]),
    lo = c(x$trends_lower[1, ], x$trends_lower[2, ]),
    hi = c(x$trends_upper[1, ], x$trends_upper[2, ]),
    trend = rep(c("trend1", "trend2"), each = length(z)),
    j_group = y,
    years = z
  )
}) %>% 
  bind_rows() %>% 
  left_join(., surv %>% select(smolt, j_group = j_group3) %>% distinct, 
            by = "j_group") %>% 
  glimpse()

trends <- rot_dat %>% 
  filter(trend == "trend1") %>% 
  ggplot(., aes(x = years, y = mean, colour = smolt)) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = smolt), alpha = 0.3) +
  facet_wrap(~j_group) +
  ggsidekick::theme_sleek() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Brood Year", y = "Shared Trend") +
  theme(legend.position = "top")

png(here::here("figs", "dfa", "bayes", "mortality", "trend_all_groups.png"), 
    height = 5.5, width = 7.5, res = 300, units = "in")
trends
dev.off()


## DFA model selection approach (NOT USED) -------------------------------------
# 
# # specify z matrix structure
# z_mm <- model.matrix(~ 0 + smolt, 
#                    data = stk_tbl %>% 
#                      filter(stock %in% keep_stks$stock))
# zmat_list <- vector(mode = "list", length = length(z_mm))
# counter <- 0
# for (j in seq_len(nrow(z_mm))) {
#   for (i in seq_len(ncol(z_mm))) {
#     counter <- counter + 1
#     if (z_mm[j, i] == 1) {
#       zmat_list[[counter]] <- paste("z", i, sep = "")
#     } else {
#       zmat_list[[counter]] <- 0
#     }
#   }
# }
# zmat <- matrix(zmat_list, nrow = n_ts, ncol = ncol(z_mm), byrow = TRUE)
# 
# # specify other MARSS parameters
# U <- x0 <- A <- "zero"
# Q <- B <- "identity"
# R <- "diagonal and unequal"
# V0 <- diag(5, ncol(zmat))
# 
# # z-score input matrix
# m_mat_z <- zscore(m_mat)
# 
# model.list <- list(A=A, U=U, x0=x0, Q=Q, B=B,
#                    R=R, Z=zmat, V0=V0)
# kemTestMod <- MARSS(m_mat_z, model = model.list, 
#                     control=list(minit=1000,maxit=5000, trace=1,
#                                  conv.test.slope.tol=0.1))
# while(kemTestMod$convergence!=0){
#   fit <- MARSS(m_mat_z, 
#                model=fit.model, 
#                control=list(maxit=4000, trace=1),
#                inits=as.matrix(coef(fit)[1]),
#                method="BFGS")
# }



## MODEL SELCTION FOR DFA STRUCTURE

m_mat_tibble <- tibble(group = levels(surv$j_group3b)) %>%
  mutate(
    m_mat = surv %>%
      filter(!is.na(M)) %>%
      group_split(j_group3b) %>%
      map(., make_mat)
  )
surv_tbl <- expand.grid(group = levels(surv$j_group3b),
                        m = c(1, 2),
                        est_ar = c(0, 1),
                        est_ma = c(0, 1),
                        est_nu = c(0, 1)) %>%
  as_tibble() %>%
  left_join(., m_mat_tibble, by = "group") %>%
  filter(group %in% kept_grps$j_group3b,
         group == "puget_oceantype",
         m == "2")

test_fits <- furrr::future_pmap(
    list(surv_tbl$m_mat, surv_tbl$m, surv_tbl$group,
         surv_tbl$est_ar,
         surv_tbl$est_ma, surv_tbl$est_nu),
    .f = function(y, m, z,
                  est_ar, est_ma, est_nu) {
      dfa_fits = fit_dfa(y = y, num_trends = m, estimate_nu = est_nu,
                         estimate_trend_ar = est_ar, estimate_trend_ma = est_ma,
                         zscore = FALSE, #only center don't scale
                         iter = 4000, chains = 4, thin = 1,
                         control = list(adapt_delta = 0.99, max_treedepth = 20))
    },
    .progress = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
  )

loo_list <- map(test_fits, bayesdfa::loo)
surv_tbl$looic <- map(loo_list, function(x) x$estimates["looic", "Estimate"]) %>%
  unlist()
