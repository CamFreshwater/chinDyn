## Bayesian DFA Figures
# March 15, 2021
# Use juvenile mortality and mean generation length DFAs to generate manuscript
# figures

library(tidyverse)
library(bayesdfa)
library(grid)
library(gridExtra)

# import juvenile mortality data
surv_tbl <- readRDS(here::here("data", "mortality_fits", "surv_tbl.RDS"))
surv_dfa <- map(surv_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_") 
  # f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "mortality_fits", f_name))
})

# import mean gen length
gen_tbl <- readRDS(here::here("data", "generation_fits", "gen_tbl.RDS"))
gen_dfa <- map(gen_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_")
  # f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_")
  readRDS(here::here("data", "generation_fits", f_name))
})

# import escapement data
# esc_tbl <- readRDS(here::here("data", "escapement_fits", "esc_tbl.RDS"))
# esc_dfa <- map(esc_tbl$group, function(y) {
#   f_name <- paste(y, "two-trend", "bayesdfa_scaled.RDS", sep = "_") 
#   readRDS(here::here("data", "escapement_fits", f_name))
# }) 


group_labs <- surv_tbl$group_labs <- gen_tbl$group_labs <- c(
  "North\nYearling", 
  "Puget\nSubyearling", 
  "Puget\nYearling", 
  "SoG\nSubyearling", 
  "South\nSubyearling"
)


# plotting functions
source(here::here("R", "functions", "plotting_functions.R"))


# PREDICTED FITS ---------------------------------------------------------------

#set seed to ensure same stocks are sampled
set.seed(345)

# predicted survival fits
surv_pred_list <- pmap(list(surv_dfa, surv_tbl$names, surv_tbl$years), 
                  fitted_preds,
                  descend_order = TRUE)
# scale colors based on observed range over entire dataset
col_ramp_surv <- surv_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

# # make one version with a legend to use in panel fig
# leg_plot <- plot_fitted_pred(surv_pred_list %>% bind_rows, 
#                              col_ramp = col_ramp_surv,
#                              leg_name = "5-year Mean of Centered Juvenile M") +
#   theme(legend.position = "top")  

#remove x_axes except for last plot
x_axes <- c(F, F, F, F, T)

surv_fit <- pmap(list(surv_pred_list, surv_tbl$group_labs, x_axes), 
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
  # rel_heights = c(rep(.195, 4), .22), #to account for text on bottom
  rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1 
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Instantaneous Juvenile Mortality Rate", 
                                gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()

cowplot::plot_grid(surv_fit_panel)

# surv_fit_panel2 <- cowplot::plot_grid(
#   cowplot::get_legend(leg_plot),
#   surv_fit_panel,
#   ncol=1, rel_heights=c(.075, .925)
# )


# predicted gen length fits
gen_pred_list <- pmap(list(gen_dfa, gen_tbl$names, gen_tbl$years), 
                       fitted_preds,
                       descend_order = FALSE)
# scale colors based on observed range over entire dataset
col_ramp_gen <- gen_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

gen_fit <- pmap(list(gen_pred_list, gen_tbl$group_labs, x_axes), 
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
  # rel_heights = c(rep(.195, 4), .22), #to account for text on bottom  
  rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Mean Age", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()

# make one version with a legend to use in panel fig
# leg_plot_g <- plot_fitted_pred(gen_pred_list[[1]], col_ramp = col_ramp_surv, 
#                              # facet_col = 5, 
#                              col_ramp_direction = 1,
#                              leg_name = "5-year Mean of Centered Mean Age") +
#   theme(legend.position = "top")  
# 
# gen_fit_panel2 <- cowplot::plot_grid(
#   cowplot::get_legend(leg_plot_g),
#   gen_fit_panel,
#   ncol=1, rel_heights=c(.075, .925)
# )

# pdf(here::here("figs", "fits_both_vars.pdf"), height = 12, width = 8)
# surv_fit_panel2
# gen_fit_panel2
# dev.off()

png(here::here("figs", "ms_figs", "mortality_fit.png"), height = 10, width = 8, 
    res = 300, units = "in")
surv_fit_panel2
dev.off()

png(here::here("figs", "ms_figs", "age_fit.png"), height = 10, width = 8, 
    res = 300, units = "in")
gen_fit_panel2
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
    thresh = ifelse(prob_above_0 > 0.90, 1, 0)
  )
sum(surv_prob$thresh) / nrow(surv_prob)


# ESTIMATES OF PARS ------------------------------------------------------------

# pull pars
pull_par_f <- function(x, group) {
  as_tibble(x$samples, rownames = "iterations") %>% 
    pivot_longer(cols = -iterations,
                 names_sep = "\\.",
                 names_to = c("chains", "parameter")) %>% 
    filter(grepl("theta", parameter) | grepl("phi", parameter) |
             grepl("nu", parameter)) %>% 
    mutate(group = group,
           trend = case_when(
             grepl("nu", parameter) ~ "Trend 1 and 2",
             grepl("[1]", parameter) ~ "Trend 1",
             grepl("[2]", parameter) ~ "Trend 2"
           ),
           parameter = case_when(
             grepl("theta", parameter) ~ "theta",
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

png(here::here("figs", "ms_figs", "mort_pars.png"), height = 7, width = 8, 
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
  plot_one_regime(y_lab = "Probability of High Mortality Regime")
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
  plot_one_regime(y_lab = "Probability of High Mortality Regime")
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

# old version combining trends of different variables

# png(here::here("figs", "ms_figs", "trend1.png"), height = 7, width = 4, 
#     res = 300, units = "in")
# plot_one_trend(trends %>% filter(trend == "Trend 1"))
# dev.off()
# 
# png(here::here("figs", "ms_figs", "trend2.png"), height = 7, width = 4, 
#     res = 300, units = "in")
# plot_one_trend(trends %>% filter(trend == "Trend 2"))
# dev.off()
# 
# png(here::here("figs", "ms_figs", "regime_trend1.png"), height = 7, width = 4, 
#     res = 300, units = "in")
# plot_one_regime(regimes %>% filter(State == "State 1" & trend == "One"))
# dev.off()
# 
# png(here::here("figs", "ms_figs", "regime_trend2.png"), height = 7, width = 4, 
#     res = 300, units = "in")
# plot_one_regime(regimes %>% filter(State == "State 1" & trend == "Two"))
# dev.off()


# ESTIMATED LOADINGS -----------------------------------------------------------

# generation length loadings
gen_load_dat <- pmap(list(rot_gen, gen_tbl$names, gen_tbl$group),
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
surv_load_dat <- pmap(list(rot_surv, surv_tbl$names, group_labs),
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

# 
# pdf(here::here("figs", "loadings_both_vars.pdf"), height = 7, width = 10)
# grid.arrange(gen_load_panel)
# grid.arrange(surv_load_panel)
# dev.off()

png(here::here("figs", "ms_figs", "mort_loadings.png"), height = 7, width = 10, 
    res = 300, units = "in")
grid.arrange(surv_load_panel)
dev.off()

png(here::here("figs", "ms_figs", "age_loadings.png"), height = 7, width = 10, 
    res = 300, units = "in")
grid.arrange(gen_load_panel)
dev.off()

