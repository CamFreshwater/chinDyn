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

gen <- gen_raw %>% 
  filter(!is.na(gen_length),
         #!j_group3 %in% c("col_streamtype", "north_oceantype", 
         # "sog_streamtype"),
         #!a_group3 == "north_streamtype"
  ) %>% 
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
         j_group1:j_group4, a_group1:a_group4, j_group4b:j_group1b) %>% 
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
  select(year, stock, gen_cent) %>% 
  pivot_wider(names_from = stock, values_from = gen_cent) %>% 
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

# q_models <- c("diagonal and unequal", "unconstrained")
# a_models <- c("scaling", "zero")

U <- "unequal"
R <- "diagonal and unequal"
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

marss_aic_tab <- purrr::map(marss_list, "out") %>% 
  bind_rows() %>% 
  arrange(AICc) %>% 
  mutate(deltaAICc = AICc - min(AICc),
         rel_like = exp(-1 * deltaAICc / 2),
         aic_weight = rel_like / sum(rel_like))

saveRDS(marss_aic_tab, here::here("data", "generation_fits",
                                  "marss_aic_tab_scalingA_centered.RDS"))

# marss_aic_tab1 <- readRDS(here::here("data", "generation_fits",
#                                     "marss_aic_tab_scalingA_centered.RDS"))
# marss_aic_tab2 <- readRDS(here::here("data", "generation_fits", 
#                                      "marss_aic_tab_zeroA_centered.RDS"))
# marss_aic_tab3 <- readRDS(here::here("data", "generation_fits", 
#                                      "marss_aic_tab_scalingA_raw.RDS"))


## FIT BAYESIAN DFA ------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
      map(., make_mat)
  ) %>% 
  filter(group %in% kept_grps$j_group3b)
gen_tbl$names <- map(gen_tbl$gen_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., gen %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
gen_tbl$years <- map(gen_tbl$gen_mat, function (x) {
  as.numeric(colnames(x))
})

# export
saveRDS(gen_tbl, here::here("data", "generation_fits", "gen_tbl.RDS"))


furrr::future_map2(
  gen_tbl$gen_mat,
  gen_tbl$group,
  .f = function (y, group) {
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE, 
      estimate_nu = TRUE, estimate_trend_ar = TRUE, estimate_trend_ma = FALSE,
      # zscore = TRUE,
      iter = 3000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    f_name <- paste(group, "one-trend", "bayesdfa_c.RDS", sep = "_")
    saveRDS(fit, here::here("data", "generation_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)
 

# read outputs
dfa_fits <- map(gen_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "generation_fits", f_name))
})


# model comparison
# loo_tbl <- tibble(group = rep(gen_tbl$group, 2),
#                   m = rep(c(2, 1), each = length(gen_tbl$group)),
#                   fits = c(dfa_fits, dfa_fits1))
# loo_tbl$loo <- map(loo_tbl$fits, bayesdfa::loo)
# loo_tbl$looic <- map(loo_tbl$loo, function(x) x$estimates["looic", "Estimate"]) %>% 
#   unlist()
# loo_tbl_out <- loo_tbl %>%
#   group_by(group) %>% 
#   mutate(min_looic = min(looic),
#             delta_loo = looic - min_looic) 
# saveRDS(loo_tbl_out, 
#         here::here("data", "generation_fits", "gen_bayes_dfa_loo_tbl.RDS"))
# two trend model heavily supported for all groups

loo_tbl_out <- readRDS(here::here("data", "generation_fits", 
                                  "gen_bayes_dfa_loo_tbl.RDS"))

# check diagnostics
map2(dfa_fits, gen_tbl$group, function (x, y) {
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



## FIT ML DFA ------------------------------------------------------------------

# # specify other MARSS parameters
# m <- 2
# x0 <- A <- "zero"
# Q <- B <- "identity"
# R <- "diagonal and unequal"
# V0 <- diag(5, m)
# 
# # z-score input matrix
# gen_mat_in <- zscore(gen_tbl$gen_mat[[2]])
# 
# model.list <- list(A=A, x0=x0, Q=Q, B=B,
#                    R=R, m = m, V0=V0)
# fit <- MARSS(gen_mat_in, model = model.list, form = "dfa", z.score = TRUE,
#              control=list(minit=1000, maxit=3000, trace=1,
#                           conv.test.slope.tol=0.1))
# while(fit$convergence != 0){
#   fit <- MARSS(gen_mat_in,
#                model=fit.model,
#                form = "dfa", z.score = TRUE,
#                control=list(maxit=2000, trace=1),
#                inits=as.matrix(coef(fit)[1]),
#                method="BFGS")
# }
# 
# 
# estZ <- coef(fit, type = "matrix")$Z
# #retrieve rotated matrix
# invH <- if (ncol(estZ) > 1) {
#   varimax(estZ)$rotmat
# } else if (ncol(estZ) == 1) {
#   1
# }
# 
# # Loadings
# rotZ <- (estZ %*% invH) %>% 
#   as.data.frame() %>% 
#   mutate(stock = row.names(gen_mat_in)) %>% 
#   left_join(., gen_tbl$names[[2]], by = "stock") %>% 
#   gather(key = "trend", value = "loading", -stock, -stock_name) %>% 
#   mutate(trend = as.numeric(as.factor(trend)),
#          stock = factor(stock, unique(stock))) %>% 
#   distinct() %>% 
#   mutate(stock = factor(stock, unique(stock))#,
#          # #invert T1 and T2 to make more intuitive
#          # loadingT = case_when(
#          #   trend %in% c("1", "2", "4") ~ (loading * -1), 
#          #   TRUE ~ loading
#          # )
#   )
# 
# ggplot(rotZ, aes(x = stock, y = loading)) +
#   geom_col() +
#   ggsidekick::theme_sleek() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   facet_wrap(~trend)
# 
# rotTrends <- (solve(invH) %*% fit$states) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   mutate(year = gen_tbl$years[[2]]) %>% 
#   gather(key = "trend", value = "est", -year) %>% 
#   mutate(trend = as.numeric(as.factor(trend))) #%>% 
# # mutate(estT = case_when(
# #   trend %in% c("1", "2", "4") ~ (est * -1), 
# #   TRUE ~ est
# # ))
# 
# ggplot(rotTrends, aes(x = year, y = est)) +
#   geom_line() +
#   ggsidekick::theme_sleek() +
#   geom_hline(yintercept = 0, colour = "red") +
#   facet_wrap(~trend)
# 
# 
# ## empty list for results
# fits <- list()
# ## extra stuff for var() calcs
# Ey <- MARSS:::MARSShatyt(fit)
# ## model params
# mod_par <- coef(fit, type="matrix")
# ZZ <- mod_par$Z
# ## number of obs ts
# nn <- dim(Ey$ytT)[1]
# ## number of time steps
# TT <- dim(Ey$ytT)[2]
# ## get the inverse of the rotation matrix
# H_inv <- varimax(ZZ)$rotmat
# ## model expectation
# fits$ex <- ZZ %*% H_inv %*% fit$states + matrix(mod_par$A, nn, TT)
# ## Var in model fits
# VtT <- MARSSkfss(fit)$VtT
# VV <- NULL
# for(tt in 1:TT) {
#   RZVZ <- mod_par$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
#   SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(fit$states[,tt,drop=FALSE])
#   VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))    
# }
# SE <- sqrt(VV)
# ## upper (1-alpha)% CI
# fits$up <- qnorm(1 - 0.05 / 2)*SE + fits$ex    
# ## lower (1-alpha)% CI
# fits$lo <- qnorm(0.05 / 2)*SE + fits$ex
# 
# temp_fits <- map2(fits, names(fits), function (x, y) {
#   rownames(x) <- gen_tbl$names[[2]]$stock
#   tt <- x %>% 
#     t() %>% 
#     as.data.frame %>% 
#     mutate(year =  gen_tbl$years[[2]], 
#            estimate = y) %>% 
#     gather(key = "stock", value = "value", -year, -estimate)
# }) %>% 
#   bind_rows() %>% 
#   left_join(., gen %>% select(stock, year, obs = gen_z), 
#             by = c("stock", "year")) %>% 
#   pivot_wider(., names_from = "estimate", values_from = "value") %>% 
#   glimpse()
# 
# ggplot(temp_fits, aes(x = year)) +
#   geom_line(aes(y = ex)) +
#   geom_point(aes(y = obs), colour = "red") +
#   geom_ribbon(aes(ymin = lo, ymax = up), linetype = 2, 
#               alpha = 0.2) +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~stock)

