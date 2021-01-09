## Fit MARSS then Bayes DFA models to instantaneous mortality data
# Jan 6 2021
# Update on surv_bayesDFA; instead of deciding groupings a priori, first perform
# model selection for groupings with MARSS (faster to converge than bayesDFA), 
# then fit group specific models with bayesdfa

library(MARSS)
library(tidyverse)

surv <- readRDS(here::here("data/salmonData/cwt_indicator_surv_clean.RDS")) %>% 
  #remove stocks that are aggregates of others on CP's advice
  # TST combines STI/TAK and AKS combines SSA and NSA
  filter(!stock %in% c("TST", "AKS"),
         year < 2017)
  

# dataframe of only stocks and juvenile groupings
stk_tbl <- surv %>% 
  select(stock, stock_name, smolt, j_group:j_group3) %>% 
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

# specify the z models based on different groups
z1 <- factor(stk_tbl$smolt)
z2 <- factor(stk_tbl$j_group)
z3 <- factor(stk_tbl$j_group2)
z4 <- factor(stk_tbl$j_group3)
z_models <- list(z1, z2, z3, z4)
names(z_models) <- c("smolt", "region", "region2", "region2-smolt")

q_models <- c("diagonal and equal", 
              "diagonal and unequal", 
              "equalvarcov",
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
mod_names = expand.grid(q = q_models, z = names(z_models)) %>% 
  mutate(name = paste(z, q, sep = "-"))
mod_tbl <- tibble(
  mod_name = mod_names$name,
  z_name = mod_names$z,
  q_name = mod_names$q,
  z_models = rep(z_models, each = 4),
  q_models = rep(q_models, times = 4)
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

top_model <- marss_list[[14]]$fit


## Bayesian DFA ----------------------------------------------------------------

library(bayesdfa)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

j_palette <- disco::disco("muted", n = length(unique(surv$j_group3)))
names(j_palette) <- unique(surv$j_group3)

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
  group_by(j_group3) %>% 
  tally() %>% 
  filter(n > 2)

#generate tbl by group
surv_tbl <- tibble(group = levels(surv$j_group3)) %>% 
  mutate(
    m_mat = surv %>% 
      filter(!is.na(M)) %>% 
      group_split(j_group3) %>% 
      map(., make_mat)
  ) %>% 
  # remove groups with less than three time series
  filter(
    group %in% kept_grps$j_group3
  )
surv_tbl$names <- map(surv_tbl$m_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., surv %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})

# fit bayesdfa
dfa_fits <- furrr::future_map(surv_tbl$m_mat, .f = fit_dfa, 
                              num_trends = 2, zscore = TRUE, 
                              iter = 1500, chains = 4, thin = 1, 
                              control = list(adapt_delta = 0.95, 
                                             max_treedepth = 20),
                              .progress = TRUE,
                              seed = TRUE)
saveRDS(dfa_fits, here::here("data", "mortality_fits", "bayesdfa_by_group.RDS"))

# temp <- surv_tbl$m_mat[[3]]
# temp_fit <- find_dfa_trends(temp, iter = 2000, kmin = 1, kmax = 5, chains = 4, 
#                             variance =  "equal",
#                             control = list(adapt_delta = 0.95, 
#                                            max_treedepth = 20))

tt <- temp_fit$best_model
rotate_trends(dfa_fits[[1]]) %>% 
  plot_trends()
plot_fitted(tt)
rotate_trends(tt) %>% 
  plot_loadings()



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

