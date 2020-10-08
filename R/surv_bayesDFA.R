## Fit Bayesian DFA models to instantaneous mortality data
# Sep 29, 2020
# Data fed to chinBayesDFA_new.Rmd for distribution to co-authors; updated
# version of fitsBayesDFA_old.R

library(tidyverse)
library(bayesdfa)
library(rstan)

#helper functions to fit and post-process DFA
source(here::here("R/functions/dfaFunctions.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

surv <- read.csv(here::here("data/salmonData/cwt_indicator_surv_clean.csv"), 
                 stringsAsFactors = FALSE) %>% 
  mutate_at(vars(stock), list(~ factor(., levels = unique(.)))) %>% 
  mutate(group = as.factor(group),
         stock = fct_reorder(stock, as.numeric(group)))
# saveRDS(surv, here::here("data", "salmonData", "clean_dfa_inputs.rds"))

# plots 
# raw_surv <- surv %>% 
#   filter(!is.na(M)) %>% 
#   ggplot(.) +
#   geom_point(aes(x = year, y = M, fill = group), shape = 21) + 
#   facet_wrap(~ fct_reorder(stock, as.numeric(group))) +
#   theme(legend.position = "top") +
#   labs(y = "Scaled M")
# pdf(here::here("figs", "raw_M_trends.pdf"))
# raw_surv
# dev.off()
# 
# ggplot(surv %>% filter(!is.na(M))) +
#   geom_histogram(aes(x = M, fill = group)) + 
#   facet_wrap(~stock) +
#   theme(legend.position = "top")


make_mat <- function(x) {
  x %>%
    select(year, stock, M) %>%
    spread(key = stock, value = M) %>%
    select(-year) %>% 
    as.matrix() %>% 
    t()
}


## FIT GLOBAL MODEL ------------------------------------------------------------

m_mat <- make_mat(surv)

# n trends = to n groups - 1
global_fit <- fit_dfa(y = m_mat, num_trends = length(unique(surv$group)) - 1, 
                       zscore = TRUE, iter = 2000, chains = 4, thin = 1,
                       control = list(adapt_delta = 0.95, max_treedepth = 20))
r_global <- rotate_trends(global_fit)

pdf(here::here("figs", "dfa", "bayes", "global_6.pdf"))
plot_trends(r_global)
plot_fitted(global_fit6)
plot_loadings(r_global) +
  ylim(-3, 3)
dev.off()


## GROUP MODEL SELECTION -------------------------------------------------------

surv_trim2 <- surv %>%
  filter(!group %in% c("north_oceantype", "south_streamtype")) %>% 
  droplevels()
surv_tbl <- tibble(group = levels(surv_trim2$group)) %>% 
  mutate(
    m_mat = surv_trim2 %>% 
      filter(!is.na(M)) %>% 
      group_split(group) %>% 
      map(., make_mat)
  )
surv_tbl$names <- map(surv_tbl$m_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., surv %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})

# # evaluate up to 3 trends per group
# surv_tbl <- readRDS(here::here("data", "dfaBayesFits", 
#                                "group_top_models.rds")) 
# trend_select_list <- lapply(surv_tbl$m_mat, function(x) {
#   find_dfa_trends(x, iter = 2500, kmin = 1, kmax = 2, chains = 4,
#                   variance =  "equal",
#                   control = list(adapt_delta = 0.95, max_treedepth = 20))
# })
# 

# map(trend_select_list, ~.$summary)
# surv_tbl$top_dfa <- map(trend_select_list, ~.$best_model)
# saveRDS(surv_tbl, here::here("data", "dfaBayesFits", "group_top_models.rds"))

map2(surv_tbl$group, surv_tbl$top_dfa, function (group_name, x) {
  file_name <- paste(group_name, "2trend_fit.pdf", sep = "_")
  
  rot_dum <- rotate_trends(x)
  
  pdf(here::here("figs", "dfa", "bayes", file_name))
  print(plot_fitted(x))
  print(plot_trends(rot_dum))
  print(plot_loadings(rot_dum) +
    ylim(-2, 2))
  dev.off()
})


## FIT TWO TREND MODELS --------------------------------------------------------

# fit models to each region, then export to Rmd
dfa_list <- map(surv_tbl$m_mat, function(x) {
  fit_dfa(y = x, num_trends = 2, zscore = TRUE, iter = 8000, chains = 1, 
          thin = 1, control = list(adapt_delta = 0.95, max_treedepth = 20))
})
surv_tbl$dfa_two_trend <- dfa_list

saveRDS(surv_tbl, here::here("data", "dfaBayesFits", "two-trend-fits.rds"))

# save figs
pmap(list(surv_tbl$group, surv_tbl$dfa_two_trend),
     function (group_name, x) {
       file_name <- paste(group_name, "2trend_fit.pdf", sep = "_")

       rot_dum <- rotate_trends(x)
       pdf(here::here("figs", "dfa", "bayes", file_name))
       print(plot_fitted(x))
       print(plot_trends(rot_dum))
       print(plot_loadings(rot_dum) +
               ylim(-2, 2))
       dev.off()
      }
)


## FIT SALISH SEA MODEL --------------------------------------------------------

# add seal covariates to Salish Sea model 
seals <- read.csv(here::here("data/salmonData/sealPopEst.csv"), 
                  stringsAsFactors = FALSE) %>% 
  rename(year = "Estimate", mean = "Mean", low = "X2.5th", up = "X97.5th", 
         reg = "Region") %>% 
  filter(reg == "SOG",
         !year < 1971) %>% 
  mutate(z_mean = scale(mean)[, 1])

# make new survival tbl with trimmed years and seal covariate
surv_trim3 <- surv %>%
  filter(agg_reg == "SS",
         !year > 2014) %>% 
  droplevels()
# ss_surv_tbl <- tibble(group = levels(surv_trim3$group)) %>% 
#   mutate(
#     m_mat = surv_trim3 %>% 
#       group_split(group) %>% 
#       map(., make_mat),
#     # reformat seal covariate to match bayes dfa input
#     seal_cov = map(m_mat, function (x) {
#       expand.grid("time" = seq(1, ncol(x)),
#                   "timeseries" = seq(1, nrow(x)),
#                   "covariate" = 1) %>%
#         mutate(value = rep(seals$z_mean, times = nrow(x)))
#     })
#   )
# ss_surv_tbl$fit1 <- map(ss_surv_tbl$m_mat, function (x) {
#   fit_dfa(y = x, num_trends = 2, iter = 3500,  chains = 4, thin = 1, 
#           zscore = TRUE, control = list(adapt_delta = 0.97, max_treedepth = 20))
# })
# ss_surv_tbl$fit2 <- map2(ss_surv_tbl$m_mat, ss_surv_tbl$seal_cov, function (x, y) {
#   fit_dfa(y = x, obs_covar = y, num_trends = 2, iter = 3500, 
#           chains = 4, thin = 1, zscore = TRUE,
#           control = list(adapt_delta = 0.97, max_treedepth = 20))
# })

# saveRDS(ss_surv_tbl, here::here("data", "dfaBayesFits", "ss_seal_models.rds"))

ss_surv_tbl <- readRDS(here::here("data", "dfaBayesFits", "ss_seal_models.rds"))

pmap(list(ss_surv_tbl$group, ss_surv_tbl$fit1, ss_surv_tbl$fit2),
     function (group_name, x, x2) {
       file_name <- paste(group_name, "2trend_fit.pdf", sep = "_")
       file_name2 <- paste(group_name, "2trend_seals_fit.pdf", sep = "_")
       
       rot_dum <- rotate_trends(x)
       pdf(here::here("figs", "dfa", "bayes", file_name))
       print(plot_fitted(x))
       print(plot_trends(rot_dum))
       print(plot_loadings(rot_dum) +
               ylim(-2, 2))
       dev.off()
       
       rot_dum2 <- rotate_trends(x2)
       pdf(here::here("figs", "dfa", "bayes", file_name2))
       print(plot_fitted(x2))
       print(plot_trends(rot_dum2))
       print(plot_loadings(rot_dum2) +
               ylim(-2, 2))
       dev.off()
     }
  )
