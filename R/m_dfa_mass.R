## Fit Bayes DFA models to mortality rate data
# May 19, 2021
# Same as m_dfa but with body mass covariate
# Note that data need to be trimmed because NAs can't be passed to covariate 
# dataframe

library(tidyverse)
library(bayesdfa)
library(grid)
library(gridExtra)

# plotting functions
source(here::here("R", "functions", "data_cleaning_functions.R"))
source(here::here("R", "functions", "plotting_functions.R"))


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

surv_raw <- readRDS(here::here("data", "salmon_data", 
                              "cwt_indicator_surv_clean.RDS")) 

surv_trim <- surv_raw %>% 
  # covariates can't be included with NAs; i.e. need to exclude stocks/years 
  # without average weight data
  filter(!is.na(M),
         !is.na(avg_weight))

# number of stocks per group
kept_grps <- surv_trim %>%
  select(stock, j_group3b) %>% 
  distinct() %>% 
  group_by(j_group3b) %>%
  tally() %>%
  filter(n > 2)

#generate tbl by group
surv_trim_tbl <- tibble(group = levels(surv_trim$j_group3b)) %>% 
  mutate(
    surv_mat = surv_trim %>% 
      group_split(j_group3b) %>% 
      map(., make_mat, resp = "M")
  ) %>% 
  filter(group %in% kept_grps$j_group3b)
surv_trim_tbl$names <- map(surv_trim_tbl$surv_mat, function (x) {
  data.frame(stock = row.names(x)) %>% 
    left_join(., surv_trim %>% select(stock, stock_name) %>% distinct(),
              by = "stock")
})
surv_trim_tbl$years <- map(surv_trim_tbl$surv_mat, function (x) {
  as.numeric(colnames(x))
})

surv_trim_tbl$cov_mat <- purrr::pmap(list(group = surv_trim_tbl$group, 
                                         mat = surv_trim_tbl$surv_mat),
                                    .f = make_cov_df, raw_dat = surv_trim,
                                    scale = "center")


# FIT MODEL -------------------------------------------------------------------- 

purrr::pmap(
  .l = list(surv_trim_tbl$surv_mat[2:5],
       surv_trim_tbl$group[2:5],
       surv_trim_tbl$cov_mat[2:5]),
  .f = function (y, group, cov) {
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = FALSE,
      obs_covar = cov,
      # estimate_trend_ar = TRUE, estimate_trend_ma = TRUE,
      iter = 3000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    f_name <- paste(group, "two-trend", "mass", "bayesdfa_c.RDS", sep = "_")
    saveRDS(fit, here::here("data", "mortality_fits", f_name))
  }#,
  # .progress = TRUE,
  # .options = furrr::furrr_options(seed = TRUE)
)

# read outputs
dfa_fits <- map(surv_trim_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "mass", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "mortality_fits", f_name))
})


# check diagnostics
map2(dfa_fits, surv_trim_tbl$group, function (x, y) {
  data.frame(neff_ratio = bayesplot::neff_ratio(x$model),
             group = y)
}) %>% 
  bind_rows()  %>% 
  ggplot(.) +
  geom_histogram(aes(x = neff_ratio)) +
  facet_wrap(~group)

map(dfa_fits, function (x) {
  as.data.frame(summary(dfa_fits[[1]]$model)$summary) %>% 
    filter(n_eff < 500 | Rhat > 1.05)
})



# rotate trends and add to gen_tbl (keep DFA separate because they're huge)
surv_trim_tbl$rot_gen <- purrr::map(dfa_fits, rotate_trends)

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
# surv_trim_tbl$regime_trend1 <- map(hmm_list_g, function(x) x[[1]]$best_model)
# surv_trim_tbl$regime_trend2 <- map(hmm_list_g, function(x) x[[2]]$best_model)


group_labs <- surv_trim_tbl$group_labs <- c(
  "North\nYearling", 
  "Puget\nSubyearling", 
  "Puget\nYearling", 
  "SoG\nSubyearling", 
  "South\nSubyearling"
)

## FIT FIG ---------------------------------------------------------------------

# predicted fits
surv_pred_list <- purrr::pmap(
  list(dfa_fits, surv_trim_tbl$names, surv_trim_tbl$years), 
  fitted_preds,
  descend_order = FALSE
)

# scale colors based on observed range over entire dataset
col_ramp_surv <- surv_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

x_axes <- c(F, F, F, F, T)
surv_fit <- pmap(list(surv_pred_list, group_labs, x_axes), 
                .f = function(x, y, z) {
                  plot_fitted_pred(x, print_x = z,
                                   col_ramp = col_ramp_surv,
                                   facet_col = 5
                  ) +
                    scale_y_continuous(
                      name = y, position = 'right', sec.axis = dup_axis(),
                      labels = scales::number_format(accuracy = 0.1)
                    ) 
                })

surv_fit_panel <- cowplot::plot_grid(
  surv_fit[[1]], surv_fit[[2]], surv_fit[[3]], surv_fit[[4]], surv_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Mortality Rates", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()

# make one version with a legend to use in panel fig
surv_plot_g <- plot_fitted_pred(surv_pred_list[[1]], col_ramp = col_ramp_surv, 
                               leg_name = "5-year Mean of Centered Mortality Rates") +
  theme(legend.position = "top")  

surv_fit_panel2 <- cowplot::plot_grid(
  cowplot::get_legend(surv_plot_g),
  surv_fit_panel,
  ncol=1, rel_heights=c(.075, .925)
)

pdf(here::here("figs", "mass_cov_figs", "mort_fits.pdf"), height = 12, width = 8)
surv_fit_panel2
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

surv_pars <- map2(dfa_fits, group_labs, pull_par_f) %>% 
  bind_rows() 

surv_par_plot <- ggplot(surv_pars, 
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

pdf(here::here("figs", "mass_cov_figs", "surv_mass_beta_est.pdf"))
surv_par_plot
dev.off()


## TREND ESTIMATES -------------------------------------------------------------

surv_trends <- pmap(
  list(surv_trim_tbl$rot_gen, surv_trim_tbl$years, group_labs), 
  .f = prep_trends
) %>% 
  bind_rows() %>% 
  mutate(var = "Juvenile Mortality Rate",
         life_history = case_when(
           grepl("Sub", group) ~ "subyearling",
           TRUE ~ "yearling"
         ))

surv_t_one <- surv_trends %>% 
  filter(trend == "Trend 1", 
         var == "Juvenile Mortality Rate") %>% 
  plot_one_trend()
surv_t_two <- surv_trends %>% 
  filter(trend == "Trend 2", 
         var == "Juvenile Mortality Rate") %>% 
  plot_one_trend()

surv_trends <- cowplot::plot_grid(surv_t_one, surv_t_two, ncol = 2)

pdf(here::here("figs", "mass_cov_figs", "surv_trends.pdf"),
    height = 8.5, width = 6)
surv_trends
dev.off()


## LOADINGS ESTIMATES ----------------------------------------------------------

# generation length loadings
surv_load_dat <- pmap(list(surv_trim_tbl$rot_gen, surv_trim_tbl$names, group_labs),
                     .f = prep_loadings) 

#make list of figures
surv_load <- map2(surv_load_dat, group_labs, plot_load, guides = FALSE) 

# make single figure to steal legend from
leg_surv_load <- plot_load(surv_load_dat[[1]], guides = TRUE) +
  ggsidekick::theme_sleek()

surv_load_panel <- 
  cowplot::plot_grid(
    surv_load[[1]], surv_load[[2]], surv_load[[3]], surv_load[[4]], surv_load[[5]],
    cowplot::get_legend(leg_surv_load),
    axis = c("lr"), align = "hv", 
    nrow = 2
  ) %>% 
  arrangeGrob(., 
              bottom = textGrob("Juvenile Mortality Rate Model Loadings",
                                gp = gpar(col = "grey30", fontsize = 12))) 

pdf(here::here("figs", "mass_cov_figs", "surv_loadings.pdf"), height = 7,
    width = 10)
grid.arrange(surv_load_panel)
dev.off()