## Fit Bayes DFA models to mean generation time data
# May 19, 2021
# Same as gen_length_dfa but with body mass covariate
# Note that data need to be trimmed because NAs can't be passed to covariate 
# dataframe

library(tidyverse)
library(bayesdfa)

source(here::here("R", "functions", "data_cleaning_functions.R"))

library(grid)
library(gridExtra)

# plotting functions
source(here::here("R", "functions", "plotting_functions.R"))



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
# covariates can't be included with NAs; i.e. need to exclude stocks/years 
# without average weight data
gen_trim <- gen %>% 
  filter(!is.na(avg_weight))

# number of stocks per group
kept_grps <- gen_trim %>%
  select(stock, j_group3b) %>% 
  distinct() %>% 
  group_by(j_group3b) %>%
  tally() %>%
  filter(n > 2)

#generate tbl by group
gen_trim_tbl <- tibble(group = levels(gen_trim$j_group3b)) %>% 
  mutate(
    gen_mat = gen_trim %>% 
      filter(!is.na(gen_length)) %>% 
      group_split(j_group3b) %>% 
      map(., make_mat)
  ) %>% 
  filter(group %in% kept_grps$j_group3b)
gen_trim_tbl$names <- map(gen_trim_tbl$gen_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., gen_trim %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
gen_trim_tbl$years <- map(gen_trim_tbl$gen_mat, function (x) {
  as.numeric(colnames(x))
})

gen_trim_tbl$cov_mat <- purrr::pmap(list(group = gen_trim_tbl$group, 
                                         mat = gen_trim_tbl$gen_mat),
                                    .f = make_cov_df, raw_dat = gen_trim,
                                    scale = "center")


# FIT MODEL -------------------------------------------------------------------- 

# furrr::future_pmap(
#   list(gen_trim_tbl$gen_mat,
#        gen_trim_tbl$group,
#        gen_trim_tbl$cov_mat),
#   .f = function (y, group, cov) {
#     fit <- fit_dfa(
#       y = y, num_trends = 2, zscore = FALSE,
#       obs_covar = cov,
#       # estimate_trend_ar = TRUE, estimate_trend_ma = TRUE,
#       iter = 3000, chains = 4, thin = 1,
#       control = list(adapt_delta = 0.99, max_treedepth = 20)
#     )
#     f_name <- paste(group, "two-trend", "mass", "bayesdfa_c.RDS", sep = "_")
#     saveRDS(fit, here::here("data", "generation_fits", f_name))
#   },
#   .progress = TRUE,
#   .options = furrr::furrr_options(seed = TRUE)
# )

# read outputs
dfa_fits <- map(gen_trim_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "mass", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "generation_fits", f_name))
})


# check diagnostics
map2(dfa_fits, gen_trim_tbl$group, function (x, y) {
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



# rotate trends and add to gen_tbl (keep DFA separate because they're huge)
gen_trim_tbl$rot_gen <- purrr::map(dfa_fits, rotate_trends)

# # hidden markov model for regimes
# hmm_list_g <- furrr::future_map(gen_tbl$rot_gen, regime_f)
# map(hmm_list_g, function(x) {
#   map(x, ~.$table)
# } )
# 
# map(hmm_list_g, function(x) {
#   map(x, function (y) plot_regime_model(y$best_model))
# } )
# 
# gen_tbl$regime_trend1 <- map(hmm_list_g, function(x) x[[1]]$best_model)
# gen_tbl$regime_trend2 <- map(hmm_list_g, function(x) x[[2]]$best_model)


group_labs <- gen_trim_tbl$group_labs <- c(
  "North\nYearling", 
  "Puget\nSubyearling", 
  "Puget\nYearling", 
  "SoG\nSubyearling", 
  "South\nSubyearling"
)

## FIT FIG ---------------------------------------------------------------------

# predicted gen length fits
gen_pred_list <- purrr::pmap(
  list(dfa_fits, gen_trim_tbl$names, gen_trim_tbl$years), 
  fitted_preds,
  descend_order = FALSE
)

# scale colors based on observed range over entire dataset
col_ramp_gen <- gen_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

x_axes <- c(F, F, F, F, T)
gen_fit <- pmap(list(gen_pred_list, group_labs, x_axes), 
                .f = function(x, y, z) {
                  plot_fitted_pred(x, print_x = z,
                                   col_ramp = col_ramp_gen,
                                   col_ramp_direction = 1,
                                   facet_col = 5
                  ) +
                    scale_y_continuous(
                      name = y, position = 'right', sec.axis = dup_axis(),
                      labels = scales::number_format(accuracy = 0.1)
                    ) 
                })

gen_fit_panel <- cowplot::plot_grid(
  gen_fit[[1]], gen_fit[[2]], gen_fit[[3]], gen_fit[[4]], gen_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(1/9, 3/9, 1/9, 2/9, 3/9),
  ncol=1
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Mean Age", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()

# make one version with a legend to use in panel fig
leg_plot_g <- plot_fitted_pred(gen_pred_list[[1]], col_ramp = col_ramp_gen, 
                               col_ramp_direction = 1,
                               leg_name = "5-year Mean of Centered Mean Age") +
  theme(legend.position = "top")  

gen_fit_panel2 <- cowplot::plot_grid(
  cowplot::get_legend(leg_plot_g),
  gen_fit_panel,
  ncol=1, rel_heights=c(.075, .925)
)

pdf(here::here("figs", "mass_cov_figs", "gen_fits.pdf"), height = 12, width = 8)
gen_fit_panel2
dev.off()


## PAR ESTIMATES ---------------------------------------------------------------

pull_par_f <- function(x, group) {
  as_tibble(x$samples, rownames = "iterations") %>% 
    pivot_longer(cols = -iterations,
                 names_sep = "\\.",
                 names_to = c("chains", "parameter")) %>% 
    filter(grepl("b_obs", parameter)) %>% 
    mutate(group = group)
}

gen_pars <- map2(dfa_fits, group_labs, pull_par_f) %>% 
  bind_rows() 

gen_par_plot <- ggplot(gen_pars, 
                       aes(x = parameter, y = value)) + 
  scale_fill_brewer(name = "", palette = "Set2") +
  geom_violin(position = position_dodge(0.3),
              draw_quantiles = 0.5) + 
  geom_hline(aes(yintercept = 0), lty = 2) + 
  coord_flip(ylim = c(-0.5, 0.5)) + 
  scale_x_discrete(limits = rev) +
  ylab("Posterior Estimates from Mean Age Model") +
  xlab("Stock Grouping") +
  ggsidekick::theme_sleek() +
  guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
  facet_wrap(~group, scales = "free_x")

pdf(here::here("figs", "mass_cov_figs", "gen_mass_beta_est.pdf"))
gen_par_plot
dev.off()


## TREND ESTIMATES ---------------------------------------------------------------

