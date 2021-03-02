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
  group_by(stock) %>% 
  mutate(gen_z = as.numeric(scale(gen_length)),
         gen_scale = as.numeric(scale(gen_length, center = TRUE, 
                                      scale = FALSE))) %>% 
  ungroup() %>% 
  droplevels()

map(colnames(gen)[c(13:20, 22:33)], function(x) {
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

gen %>% 
  ggplot(.) +
  geom_point(aes(x = year, y = gen_z, fill = a_group2), shape = 21) +
  geom_point(aes(x = year, y = gen_scale, fill = a_group2), shape = 24) +
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


## FIT ML DFA ------------------------------------------------------------------

# specify other MARSS parameters
m <- 2
x0 <- A <- "zero"
Q <- B <- "identity"
R <- "diagonal and unequal"
V0 <- diag(5, m)

# z-score input matrix
gen_mat_in <- zscore(gen_tbl$gen_mat[[2]])

model.list <- list(A=A, x0=x0, Q=Q, B=B,
                   R=R, m = m, V0=V0)
fit <- MARSS(gen_mat_in, model = model.list, form = "dfa", z.score = TRUE,
                    control=list(minit=1000, maxit=3000, trace=1,
                                 conv.test.slope.tol=0.1))
while(fit$convergence != 0){
  fit <- MARSS(gen_mat_in,
               model=fit.model,
               form = "dfa", z.score = TRUE,
               control=list(maxit=2000, trace=1),
               inits=as.matrix(coef(fit)[1]),
               method="BFGS")
}


estZ <- coef(fit, type = "matrix")$Z
#retrieve rotated matrix
invH <- if (ncol(estZ) > 1) {
  varimax(estZ)$rotmat
} else if (ncol(estZ) == 1) {
  1
}

# Loadings
rotZ <- (estZ %*% invH) %>% 
  as.data.frame() %>% 
  mutate(stock = row.names(gen_mat_in)) %>% 
  left_join(., gen_tbl$names[[2]], by = "stock") %>% 
  gather(key = "trend", value = "loading", -stock, -stock_name) %>% 
  mutate(trend = as.numeric(as.factor(trend)),
         stock = factor(stock, unique(stock))) %>% 
  distinct() %>% 
  mutate(stock = factor(stock, unique(stock))#,
         # #invert T1 and T2 to make more intuitive
         # loadingT = case_when(
         #   trend %in% c("1", "2", "4") ~ (loading * -1), 
         #   TRUE ~ loading
         # )
         )

ggplot(rotZ, aes(x = stock, y = loading)) +
  geom_col() +
  ggsidekick::theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~trend)

rotTrends <- (solve(invH) %*% fit$states) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(year = gen_tbl$years[[2]]) %>% 
  gather(key = "trend", value = "est", -year) %>% 
  mutate(trend = as.numeric(as.factor(trend))) #%>% 
  # mutate(estT = case_when(
  #   trend %in% c("1", "2", "4") ~ (est * -1), 
  #   TRUE ~ est
  # ))

ggplot(rotTrends, aes(x = year, y = est)) +
  geom_line() +
  ggsidekick::theme_sleek() +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~trend)


## empty list for results
fits <- list()
## extra stuff for var() calcs
Ey <- MARSS:::MARSShatyt(fit)
## model params
mod_par <- coef(fit, type="matrix")
ZZ <- mod_par$Z
## number of obs ts
nn <- dim(Ey$ytT)[1]
## number of time steps
TT <- dim(Ey$ytT)[2]
## get the inverse of the rotation matrix
H_inv <- varimax(ZZ)$rotmat
## model expectation
fits$ex <- ZZ %*% H_inv %*% fit$states + matrix(mod_par$A, nn, TT)
## Var in model fits
VtT <- MARSSkfss(fit)$VtT
VV <- NULL
for(tt in 1:TT) {
  RZVZ <- mod_par$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
  SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(fit$states[,tt,drop=FALSE])
  VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))    
}
SE <- sqrt(VV)
## upper (1-alpha)% CI
fits$up <- qnorm(1 - 0.05 / 2)*SE + fits$ex    
## lower (1-alpha)% CI
fits$lo <- qnorm(0.05 / 2)*SE + fits$ex

temp_fits <- map2(fits, names(fits), function (x, y) {
  rownames(x) <- gen_tbl$names[[2]]$stock
  tt <- x %>% 
    t() %>% 
    as.data.frame %>% 
    mutate(year =  gen_tbl$years[[2]], 
           estimate = y) %>% 
    gather(key = "stock", value = "value", -year, -estimate)
}) %>% 
  bind_rows() %>% 
  left_join(., gen %>% select(stock, year, obs = gen_z), 
            by = c("stock", "year")) %>% 
  pivot_wider(., names_from = "estimate", values_from = "value") %>% 
  glimpse()

ggplot(temp_fits, aes(x = year)) +
  geom_line(aes(y = ex)) +
  geom_point(aes(y = obs), colour = "red") +
  geom_ribbon(aes(ymin = lo, ymax = up), linetype = 2, 
              alpha = 0.2) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)



## FIT BAYESIAN DFA ------------------------------------------------------------

library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

a_palette <- disco::disco("muted", n = length(unique(gen$j_group3b)))
names(a_palette) <- unique(gen$j_group3b)

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

# specify one trend if there are less than 4 time series, otherwise 2
# n_trend_list <- ifelse(unlist(map(gen_tbl$gen_mat, nrow)) < 4, 1, 2) 
# 
# dfa_fits <- furrr::future_map2(
#   gen_tbl$gen_mat,
#   n_trend_list,
#   .f = function (y, n_trend) {
#     fit_dfa(y = y, num_trends = n_trend, zscore = TRUE,
#             iter = 4250, chains = 4, thin = 1,
#             control = list(adapt_delta = 0.99, max_treedepth = 20))
#   },
#   .progress = TRUE,
#   .options = furrr::furrr_options(seed = TRUE)
# )
# 
# # save outputs
# map2(dfa_fits, gen_tbl$group, function(x, y) {
#   f_name <- paste(y, "bayesdfa.RDS", sep = "_")
#   saveRDS(x, here::here("data", "generation_fits", f_name))
# })

# read outputs
dfa_fits2 <- map(gen_tbl$group, function(y) {
  f_name <- paste(y, "bayesdfa.RDS", sep = "_") 
  readRDS(here::here("data", "generation_fits", f_name))
})

# check diagnostics
library(bayesplot)
posterior_m1 <- as.array(dfa_fits[[1]]$samples)
lp_m1 <- log_posterior(dfa_fits[[1]]$model)
np_m1 <- nuts_params(dfa_fits[[1]]$model)
color_scheme_set("darkgray")
mcmc_parcoord(posterior_m1, np = np_m1)

p_m1 <- posterior_samples(dfa_fits[[1]]$model, add_chain = T) %>% 
  select(-lp__, -iter, -contains("b_"), -contains("sd_"), 
         -contains("Intercept_"))



map(dfa_fits2, function (x) {
  
})


m1 <- dfa_fits[[2]]$model

tt <- find_dfa_trends(y = gen_tbl$gen_mat[[2]], kmin = 1, kmax = 2, 
                      zscore = TRUE, iter = 4250, chains = 4, thin = 1, 
                      convergence_threshold = 1.05, 
                      variance = c("equal", "unequal"),
                      control = list(adapt_delta = 0.99, max_treedepth = 20))
tt1 <- fit_dfa(y = gen_tbl$gen_mat[[2]], num_trends = 1, 
               est_correlation = FALSE,
               zscore = FALSE, iter = 4000, chains = 4, thin = 1, 
               control = list(adapt_delta = 0.99, max_treedepth = 20))
tt2 <- fit_dfa(y = gen_tbl$gen_mat[[2]], num_trends = 2, 
               est_correlation = FALSE,
               zscore = FALSE, iter = 4250, chains = 4, thin = 1, 
               control = list(adapt_delta = 0.99, max_treedepth = 20))


neff_vals <- bayesplot::neff_ratio(tt2$model)
bayesplot::mcmc_neff_data(neff_vals) %>% 
  print(n = Inf)

rhat_vals <- bayesplot::rhat(tt2$model)
bayesplot::mcmc_rhat_data(rhat_vals) %>% 
  print(n = Inf)

n_eff_list <- map(dfa_fits2, function (x) {
  neff_vals <- neff_ratio(x$model)
  mcmc_neff_data(neff_vals) %>% 
    ggplot(.) +
    geom_histogram(aes(x = value))
})

plot_fitted_bayes(tt1, gen_tbl$names[[2]]$stock_name, gen_tbl$years[[2]])

as.data.frame(summary(tt2$model)$summary) %>% 
  filter(n_eff < 500)




# PLOT BAYESIAN DFA ------------------------------------------------------------

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
