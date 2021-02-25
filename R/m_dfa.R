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
         year < 2017#,
         # !j_group3 %in% c("col_streamtype", "north_oceantype", 
         #                  "sog_streamtype"),
         # !a_group3 == "north_streamtype"
         ) %>% 
  droplevels()
  

# dataframe of only stocks and juvenile groupings
stk_tbl <- surv %>% 
  ungroup() %>% 
  select(stock, stock_name, smolt, run, 
         j_group1:j_group4, a_group1:a_group4, j_group4b:a_group1c) %>% 
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
  fit <- MARSS(m_mat, model = fit_model,
               silent = FALSE, control = list(minit = 100, maxit = 500))
  if (fit$convergence != 0){
    fit <- MARSS(m_mat, 
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

saveRDS(marss_list, here::here("data", "mortality_fits",
                               "marss_selection_fits.RDS"))
saveRDS(marss_aic_tab, here::here("data", "mortality_fits",
                               "marss_aic_tab.RDS"))


## Bayesian DFA ----------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

j_palette <- disco::disco("muted", n = length(unique(surv$j_group3b)))
names(j_palette) <- unique(surv$j_group3b)

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

#generate tbl by group
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

# fit bayesdfa
dfa_fits <- furrr::future_map(surv_tbl$m_mat, .f = fit_dfa,
                              num_trends = 2, zscore = TRUE,
                              iter = 2500, chains = 4, thin = 1,
                              control = list(adapt_delta = 0.92,
                                             max_treedepth = 20),
                              .progress = TRUE,
                              seed = TRUE)

map2(dfa_fits, surv_tbl$group, function(x, y) {
  f_name <- paste(y, "bayesdfa.RDS", sep = "_")
  saveRDS(x, here::here("data", "mortality_fits", f_name))
})

# read outputs
dfa_fits2 <- map(surv_tbl$group, function(y) {
  f_name <- paste(y, "bayesdfa.RDS", sep = "_") 
  readRDS(here::here("data", "mortality_fits", f_name))
})
rot_list <- map(dfa_fits2, rotate_trends)

# make plots 
source(here::here("R", "functions", "plot_fitted_bayes.R"))
fit_list <- pmap(list(dfa_fits2, surv_tbl$names, surv_tbl$years), 
                 function(x, names, obs_years) {
                   plot_fitted_bayes(x, names = names$stock_name, 
                                     years = obs_years) +
                     scale_x_continuous(breaks = c(1970, 1990, 2010), 
                                        limits = c(min(surv$year), 
                                                   max(surv$year)))
                 })
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

