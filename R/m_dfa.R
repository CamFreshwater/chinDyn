## Fit MARSS then DFA models to instantaneous mortality data
# Jan 6 2021
# Update on surv_bayesDFA; instead of deciding groupings a priori, first perform
# model selection for groupings with MARSS (faster to converge than bayesDFA), 
# then fit group specific models with bayesdfa

library(MARSS)
library(tidyverse)

surv <- read.csv(here::here("data/salmonData/cwt_indicator_surv_clean.csv")) %>% 
  mutate_at(vars(stock), list(~ factor(., levels = unique(.)))) %>% 
  mutate(year = ifelse(smolt == "streamtype", brood_year + 2, 
                       brood_year + 1),
         smolt = as.factor(smolt),
         j_group = as.factor(j_group),
         j_group2 = as.factor(j_group2),
         j_group3 = as.factor(j_group3),
         stock_f = fct_reorder(stock, as.numeric(j_group))) %>% 
  filter(!year > 2017)

# dataframe of only stocks and juvenile groupings
stk_tbl <- surv %>% 
  select(stock, stock_name, smolt, j_group:j_group3) %>% 
  distinct()

#test surv dataset
keep_stks <- surv %>% 
  filter(year > 1980,
         year < 2010,
         !is.na(M)) %>% 
  group_by(stock) %>% 
  tally() %>% 
  filter(n == 29)

# subset of stocks for test run
# surv_sub <- surv %>% 
#   filter(#year > 1980,
#          #year < 2010,
#          stock %in% keep_stks$stock)
# stk_tbl_sub <- stk_tbl %>% 
#   filter(stock %in% keep_stks$stock)


## MARSS MODEL RUNS -----------------------------------------------------------

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

#generate tbl by group
surv_tbl <- tibble(group = levels(surv$j_group3)) %>% 
  mutate(
    m_mat = surv %>% 
      filter(!is.na(M)) %>% 
      group_split(j_group3) %>% 
      map(., make_mat)
  )
surv_tbl$names <- map(surv_tbl$m_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., surv %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})

#function to fit 2 trends unless n_groups < 3
pick_and_fit <- function(x) {
  nn <- ifelse(nrow(x) < 3, 1, 2)
  fit_dfa(y = x, num_trends = nn, zscore = TRUE, iter = 2000, chains = 4, 
          thin = 1, control = list(adapt_delta = 0.95, max_treedepth = 20))
}

temp <- surv_tbl$m_mat[[3]]
temp_fit <- find_dfa_trends(temp, iter = 2000, kmin = 1, kmax = 5, chains = 4, 
                            variance =  "equal",
                            control = list(adapt_delta = 0.95, max_treedepth = 20))


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

