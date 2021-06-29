## Bayesian DFA Figures
# March 15, 2021
# Use juvenile mortality and mean generation length DFAs to generate manuscript
# figures
# Updated June 2 to use logit survival
# Updated June 28 to use AR only model

library(tidyverse)
library(bayesdfa)
library(grid)
library(gridExtra)

raw_data <- readRDS(here::here("data/salmon_data/cwt_indicator_surv_clean.RDS"))

# import juvenile mortality data
surv_tbl <- readRDS(here::here("data", "survival_fits", "surv_tbl.RDS"))
surv_dfa <- map(surv_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_") 
  # f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "survival_fits", f_name))
})

# import mean gen length
gen_tbl <- readRDS(here::here("data", "generation_fits", "gen_tbl.RDS"))
gen_dfa <- map(gen_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_")
  # f_name <- paste(y, "one-trend", "ar", "bayesdfa_c.RDS", sep = "_")
  readRDS(here::here("data", "generation_fits", f_name))
})[c(1, 4, 2, 3, 5)]

gen_dfa <- gen_dfa1 

# make labels for stock groupings and edit stock names to fit
group_labs <- c(
  "North\nYearling", 
  "SoG\nSubyearling", 
  "Puget\nSubyearling", 
  "Puget\nYearling", 
  "South\nSubyearling"
)

surv_names <- surv_tbl %>% 
  select(group, names) %>% 
  unnest(cols = c(names)) %>%
  mutate(
    stock_name = fct_recode(
      as.factor(stock_name), 
      "Behm Canal Spring" = "Alaska Hatchery Springs - Behm Canal",
      "U.W. Accelerated" = "University of Washington Accelerated",
      "G. Adams Fall Fingerling" = "George Adams Fall Fingerling",
      "South P.S. Fall Fingerling" = "South Puget Sound Fall Fingerling",
      "South P.S. Fall Yearling" = "South Puget Sound Fall Yearling",
      "Lower Shuswap R. Summer" = "Lower Shuswap River Summer",
      "Middle Shuswap R. Summer" = "Middle Shuswap River Summer",
      "Ok.-Sim. Summer" = "Okanagan-Similkameen Summer",
      "Columbia Upriver Bright" = "Columbia River Upriver Bright"
    )
  ) %>%
  nest(names = c(stock, stock_name))

gen_names <- gen_tbl %>% 
  select(group, names) %>% 
  unnest(cols = c(names)) %>%
  mutate(
    stock_name = fct_recode(
      as.factor(stock_name), 
      "U.W. Accelerated" = "University of Washington Accelerated",
      "G. Adams Fall Fingerling" = "George Adams Fall Fingerling",
      "South P.S. Fall Fingerling" = "South Puget Sound Fall Fingerling",
      "South P.S. Fall Yearling" = "South Puget Sound Fall Yearling",
      "Lower Shuswap R. Summer" = "Lower Shuswap River Summer",
      "Middle Shuswap R. Summer" = "Middle Shuswap River Summer",
      "Ok.-Sim. Summer" = "Okanagan-Similkameen Summer",
      "Columbia Upriver Bright" = "Columbia River Upriver Bright"
    )
  ) %>%
  nest(names = c(stock, stock_name))


# plotting functions
source(here::here("R", "functions", "plotting_functions.R"))


# PREDICTED FITS ---------------------------------------------------------------

#set seed to ensure same stocks are sampled
set.seed(345)

#remove x_axes except for last plot
x_axes <- c(F, F, F, F, T)

# predicted survival fits
surv_pred_list <- pmap(list(surv_dfa, surv_names$names, surv_tbl$years), 
                  fitted_preds,
                  descend_order = FALSE)
# reorder to match group_labs
# surv_pred_list <- surv_pred_list1[c(1, 4, 2, 3, 5)]

# scale colors based on observed range over entire dataset
col_ramp_surv <- surv_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

surv_fit <- pmap(list(surv_pred_list, group_labs, x_axes), 
                 .f = function(x, y, z) {
                   plot_fitted_pred(x, print_x = z,
                                    col_ramp = col_ramp_surv,
                                    facet_col = 5
                   ) +
                     scale_y_continuous(
                       name = y, position = 'right', sec.axis = dup_axis(),
                       labels = scales::number_format(accuracy = 1)
                     ) 
})

surv_fit_panel <- cowplot::plot_grid(
  surv_fit[[1]], surv_fit[[2]], surv_fit[[3]], surv_fit[[4]], surv_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 2/11, 3/11, 1/11, 3/11),
  ncol=1 
) %>% 
  arrangeGrob(
    ., 
    left = textGrob("Centered Juvenile Survival Rate (logit transformed)", 
                    gp = gpar(col = "grey30", fontsize = 12),
                    rot = 90)
  ) %>% 
  grid.arrange()

# real space fits
real_surv_pred_list <- purrr::map(surv_pred_list, function (x) {
  # calculate estimated uncentered survival rate in real space 
  x %>% 
    left_join(., 
              raw_data %>% 
                select(stock, Time = brood_year, survival),
              by = c("stock", "Time")) %>% 
    group_by(ID) %>% 
    mutate(obs_mean_logit = mean(qlogis(survival), na.rm = T)) %>%
    ungroup() %>% 
    # uncenter predictions and calculate in real space
    mutate(uncent_mean_logit = mean + obs_mean_logit,
           uncent_mean = plogis(uncent_mean_logit),
           uncent_lo = plogis(obs_mean_logit + lo),
           uncent_hi = plogis(obs_mean_logit + hi),
           last_mean = plogis(obs_mean_logit + last_mean))
})


# one version of figure with full ylimit range
# real_surv_ylims <- c(0, max(raw_data$survival, na.rm = T))
real_surv_fit <- pmap(list(real_surv_pred_list, group_labs, x_axes), 
                 .f = function(x, y, z) {
                   plot_fitted_pred_real(x, 
                                         print_x = z,
                                         facet_col = 5
                   ) +
                     scale_y_continuous(
                       name = y, position = 'right', sec.axis = dup_axis()) 
                 })
real_surv_fit_panel <- cowplot::plot_grid(
  real_surv_fit[[1]], real_surv_fit[[2]], real_surv_fit[[3]], real_surv_fit[[4]], 
  real_surv_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 2/11, 3/11, 1/11, 3/11),
  ncol=1 
) %>% 
  arrangeGrob(., 
              left = textGrob("Survival Rate", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()


# predicted gen length fits
gen_pred_list <- pmap(list(gen_dfa, gen_names$names, gen_tbl$years), 
                       fitted_preds,
                       descend_order = FALSE)

# scale colors based on observed range over entire dataset
col_ramp_gen <- gen_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

gen_fit <- pmap(list(gen_pred_list, group_labs, x_axes), 
                .f = function(x, y, z) {
                  plot_fitted_pred(x, print_x = z,
                                   col_ramp = col_ramp_gen,
                                   # col_ramp_direction = 1,
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
  rel_heights = c(2/11, 2/11, 3/11, 1/11, 3/11),
  ncol=1
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Mean Age", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()


png(here::here("figs", "ms_figs", "survival_fit_ar1.png"), height = 11, width = 9, 
    res = 300, units = "in")
cowplot::plot_grid(surv_fit_panel)
dev.off()

png(here::here("figs", "ms_figs", "survival_fit_real_ar1.png"), height = 11, width = 9, 
    res = 300, units = "in")
cowplot::plot_grid(real_surv_fit_panel)
dev.off()

png(here::here("figs", "ms_figs", "age_fit_ar1.png"), height = 11, width = 9, 
    res = 300, units = "in")
cowplot::plot_grid(gen_fit_panel)
dev.off()


# CALCULATE FINAL FIVE YEAR MEANS ----------------------------------------------

# generation length means in last five years
gen_prob <- map2(gen_dfa, gen_tbl$names, final_prob, n_years = 5) %>% 
  bind_rows() %>% 
  mutate(
    thresh = ifelse(prob_below_0 > 0.90, 1, 0)
  )
sum(gen_prob$thresh) / nrow(gen_prob)

surv_prob <- map2(surv_dfa, surv_tbl$names, final_prob, n_years = 5) %>% 
  bind_rows() %>% 
  mutate(
    thresh = ifelse(prob_below_0 > 0.90, 1, 0)
  )
sum(surv_prob$thresh) / nrow(surv_prob)


# ESTIMATES OF PARS ------------------------------------------------------------

# pull pars
pull_par_f <- function(x, group) {
  as_tibble(x$samples, rownames = "iterations") %>% 
    pivot_longer(cols = -iterations,
                 names_sep = "\\.",
                 names_to = c("chains", "parameter")) %>% 
    filter(#grepl("theta", parameter) | 
      grepl("phi", parameter) |
             grepl("nu", parameter)) %>% 
    mutate(group = group,
           trend = case_when(
             grepl("nu", parameter) ~ "Trend 1 and 2",
             grepl("[1]", parameter) ~ "Trend 1",
             grepl("[2]", parameter) ~ "Trend 2"
           ),
           parameter = case_when(
             # grepl("theta", parameter) ~ "theta",
             grepl("phi", parameter) ~ "phi",
             grepl("nu", parameter) ~ "nu",
             TRUE ~ parameter
           ),
           y_int = ifelse(parameter == "nu", 10, 0))
}

# clean parameter samples
surv_pars <- map2(surv_dfa, group_labs, pull_par_f) %>% 
  bind_rows() 
gen_pars <- map2(gen_dfa, group_labs, pull_par_f) %>% 
  bind_rows() 

surv_par_plot <- ggplot(surv_pars, 
                         aes(x = group, y = value, fill = trend)) + 
  scale_fill_brewer(name = "", palette = "Set2") +
  geom_violin(position = position_dodge(0.3),
              draw_quantiles = 0.5) + 
  geom_hline(aes(yintercept = y_int), lty = 2) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  ylab("Posterior Estimates from Mortality Model") +
  xlab("Stock Grouping") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggsidekick::theme_sleek() +
  guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
  facet_wrap(~parameter, scales = "free_x")

gen_par_plot <- ggplot(gen_pars, 
                       aes(x = group, y = value, fill = trend)) + 
  scale_fill_brewer(name = "", palette = "Set2") +
  geom_violin(position = position_dodge(0.3),
              draw_quantiles = 0.5) + 
  geom_hline(aes(yintercept = y_int), lty = 2) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  ylab("Posterior Estimates from Mean Age Model") +
  xlab("Stock Grouping") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggsidekick::theme_sleek() +
  guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
  facet_wrap(~parameter, scales = "free_x")

# pdf(here::here("figs", "pars_both_vars.pdf"))
# surv_par_plot
# gen_par_plot
# dev.off()

png(here::here("figs", "ms_figs", "surv_pars.png"), height = 7, width = 8, 
    res = 300, units = "in")
surv_par_plot
dev.off()

png(here::here("figs", "ms_figs", "age_pars.png"), height = 7, width = 8, 
    res = 300, units = "in")
gen_par_plot
dev.off()


# ESTIMATED TRENDS -------------------------------------------------------------

# prep dataframes for each
surv_trends <- pmap(
  list(surv_tbl$rot_surv, surv_tbl$years, group_labs), 
  .f = prep_trends
  ) %>% 
  bind_rows() %>% 
  mutate(var = "Juvenile Mortality Rate")

gen_trends <- pmap(
  list(gen_tbl$rot_gen, gen_tbl$years, group_labs), 
  .f = prep_trends
) %>% 
  bind_rows() %>% 
  mutate(var = "Mean Age")

trends <- rbind(surv_trends, gen_trends) %>%
  mutate(var = as.factor(var),
         trend = as.factor(trend),
         life_history = case_when(
           grepl("Sub", group) ~ "subyearling",
           TRUE ~ "yearling"
         ),
         group = fct_relevel(as.factor(group), "North\nYearling",
                             "SoG\nSubyearling", "Puget\nYearling",
                             "Puget\nSubyearling", "South\nSubyearling"),
         var = fct_relevel(as.factor(var), "Juvenile Mortality Rate",
                           "Mean Age")
         )


# new version combines trends and regimes
surv_regimes1 <- pmap(
  list(regime_model = surv_tbl$regime_trend1, years = surv_tbl$years, 
       group = group_labs), 
  .f = prep_regime
) %>% 
  bind_rows() %>%
  mutate(trend = "One") 
surv_regimes2 <- pmap(
  list(regime_model = surv_tbl$regime_trend2, years = surv_tbl$years, 
       group = group_labs), 
  .f = prep_regime
) %>% 
  bind_rows() %>% 
  mutate(trend = "Two")
surv_regimes <- rbind(surv_regimes1, surv_regimes2) %>% 
  mutate(var = "Juvenile Mortality Rate") 

gen_regimes1 <- pmap(
  list(regime_model = gen_tbl$regime_trend1, years = gen_tbl$years, 
       group = group_labs), 
  .f = prep_regime,
  flip_regimes = TRUE
) %>% 
  bind_rows() %>%
  mutate(trend = "One") 
gen_regimes2 <- pmap(
  list(regime_model = gen_tbl$regime_trend2, years = gen_tbl$years, 
       group = group_labs), 
  .f = prep_regime,
  flip_regimes = TRUE
) %>% 
  bind_rows() %>% 
  mutate(trend = "Two")
gen_regimes <- rbind(gen_regimes1, gen_regimes2) %>% 
  mutate(var = "Mean Age") 

regimes <- rbind(surv_regimes, gen_regimes) %>% 
  mutate(var = as.factor(var),
         state = as.factor(State), 
         life_history = case_when(
           grepl("Sub", group) ~ "subyearling",
           TRUE ~ "yearling"
         ),
         group = fct_relevel(as.factor(group), "North\nYearling", 
                             "SoG\nSubyearling", "Puget\nYearling",
                             "Puget\nSubyearling", "South\nSubyearling")#,
         # var = fct_relevel(as.factor(var), "Juvenile Mortality Rate", 
         #                   "Mean Age")
  ) 

# dummy version just for legend
leg_plot <- trends %>% 
  filter(trend == "Trend 1", 
         var == "Juvenile Mortality Rate") %>% 
  plot_one_trend() +
  theme(legend.position = "top")
  
# first survival trend
surv_t_one <- trends %>% 
  filter(trend == "Trend 1", 
         var == "Juvenile Mortality Rate") %>% 
  plot_one_trend()
surv_r_one <- regimes %>% 
  filter(trend == "One", 
         State == "State 2",
         var == "Juvenile Mortality Rate") %>% 
  plot_one_regime(y_lab = "Probability of High SurvivalRegime")
surv_one_panel <- cowplot::plot_grid(surv_t_one, surv_r_one, ncol = 2)

# second survival trend
surv_t_two <- trends %>% 
  filter(trend == "Trend 2", 
         var == "Juvenile Mortality Rate") %>% 
  plot_one_trend()
surv_r_two <- regimes %>% 
  filter(trend == "Two", 
         State == "State 2",
         var == "Juvenile Mortality Rate") %>% 
  plot_one_regime(y_lab = "Probability of High Survival Regime")
surv_two_panel <- cowplot::plot_grid(surv_t_two, surv_r_two, ncol = 2)

# first age trend
gen_t_one <- trends %>% 
  filter(trend == "Trend 1", 
         var == "Mean Age") %>% 
  plot_one_trend()
gen_r_one <- regimes %>% 
  filter(trend == "One", 
         State == "State 1",
         var == "Mean Age") %>% 
  plot_one_regime(y_lab = "Probability of Low Mean Age Regime")
gen_one_panel <- cowplot::plot_grid(gen_t_one, gen_r_one, ncol = 2)

# second age trend
gen_t_two <- trends %>% 
  filter(trend == "Trend 2", 
         var == "Mean Age") %>% 
  plot_one_trend()
gen_r_two <- regimes %>% 
  filter(trend == "Two", 
         State == "State 1",
         var == "Mean Age") %>% 
  plot_one_regime(y_lab = "Probability of Low Mean Age Regime")
gen_two_panel <- cowplot::plot_grid(gen_t_two, gen_r_two, ncol = 2)

#output plots
plot_out <- function(x) {
  cowplot::plot_grid(
    cowplot::get_legend(leg_plot),
    x,
    ncol=1, rel_heights=c(.075, .925)
  )
}

png(here::here("figs", "ms_figs", "trend_regime_surv1.png"), 
    height = 8.5, width = 6, res = 300, units = "in")
plot_out(surv_one_panel)
dev.off()

png(here::here("figs", "ms_figs", "trend_regime_surv2.png"), 
    height = 8.5, width = 6, res = 300, units = "in")
plot_out(surv_two_panel)
dev.off()

png(here::here("figs", "ms_figs", "trend_regime_gen1.png"), 
    height = 8.5, width = 6, res = 300, units = "in")
plot_out(gen_one_panel)
dev.off()

png(here::here("figs", "ms_figs", "trend_regime_gen2.png"), 
    height = 8.5, width = 6, res = 300, units = "in")
plot_out(gen_two_panel)
dev.off()


# ESTIMATED LOADINGS -----------------------------------------------------------

# generation length loadings
gen_load_dat <- pmap(list(gen_tbl$rot_gen, gen_tbl$names, gen_tbl$group),
                      .f = prep_loadings) 

#make list of figures
gen_load <- map2(gen_load_dat, group_labs, plot_load, guides = FALSE) 

# make single figure to steal legend from
leg_gen_load <- plot_load(gen_load_dat[[1]], guides = TRUE) +
  ggsidekick::theme_sleek()

gen_load_panel <- 
  cowplot::plot_grid(
  gen_load[[1]], gen_load[[2]], gen_load[[3]], gen_load[[4]], gen_load[[5]],
  cowplot::get_legend(leg_gen_load),
  axis = c("lr"), align = "hv", 
  nrow = 2
) %>% 
  arrangeGrob(., 
              bottom = textGrob("Mean Age Model Loadings",
                                gp = gpar(col = "grey30", fontsize = 12))) 
  
# survival loadings
surv_load_dat <- pmap(list(surv_tbl$rot_surv, surv_tbl$names, 
                           group_labs),
                     .f = prep_loadings) 

#make list of figures
surv_load <- map2(surv_load_dat, group_labs, plot_load, guides = FALSE, 
                  y_lims = c(-1.5, 1.5))

# make single figure to steal legend from
leg_surv_load <- plot_load(surv_load_dat[[1]], guides = TRUE) +
  ggsidekick::theme_sleek()

surv_load_panel <- 
  cowplot::plot_grid(
    surv_load[[1]], surv_load[[2]], surv_load[[3]], surv_load[[4]], 
    surv_load[[5]],
    cowplot::get_legend(leg_surv_load),
    axis = c("lr"), align = "hv", 
    nrow = 2
  ) %>% 
  arrangeGrob(., 
              bottom = textGrob("Mortality Model Loadings",
                                gp = gpar(col = "grey30", fontsize = 12)))

png(here::here("figs", "ms_figs", "surv_loadings.png"), height = 7, width = 10, 
    res = 300, units = "in")
grid.arrange(surv_load_panel)
dev.off()

png(here::here("figs", "ms_figs", "age_loadings.png"), height = 7, width = 10, 
    res = 300, units = "in")
grid.arrange(gen_load_panel)
dev.off()

