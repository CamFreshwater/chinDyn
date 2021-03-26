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


surv_tbl$group_labs <- gen_tbl$group_labs <- c("North\nYearling", 
                                               "Puget\nSubyearling", 
                                               "Puget\nYearling", 
                                               "SoG\nSubyearling", 
                                               "South\nSubyearling")


# plotting functions
source(here::here("R", "functions", "plotting_functions.R"))


# PREDICTED FITS ---------------------------------------------------------------

#set seed to ensure same stocks are sampled
set.seed(345)

# predicted survival fits
surv_pred_list <- pmap(list(surv_dfa, surv_tbl$names, surv_tbl$years), 
                  fitted_preds,
                  descend_order = TRUE,
                  subset = 5)
# scale colors based on observed range over entire dataset
col_ramp_surv <- surv_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

# make one version with a legend to use in panel fig
leg_plot <- plot_fitted_pred(surv_pred_list[[1]], surv_tbl$group_labs[[1]],
                             col_ramp = col_ramp_surv, 
                             leg_name = "5-year Mean of Centered Juvenile M") +
  theme(legend.position = "top")  

#remove x_axes except for last plot
x_axes <- c(F, F, F, F, T)

surv_fit <- pmap(list(surv_pred_list, surv_tbl$group_labs, x_axes), 
                 .f = function(x, y, z) {
                   plot_fitted_pred(x, print_x = z,
                                    col_ramp = col_ramp_surv
                                    # , facet_col = 7
                   ) +
                     scale_y_continuous(
                       name = y, position = 'right', sec.axis = dup_axis(),
                       labels = scales::number_format(accuracy = 1)
                     ) 
})

surv_fit_panel <- cowplot::plot_grid(
  surv_fit[[1]], surv_fit[[2]], surv_fit[[3]], surv_fit[[4]], surv_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(rep(.195, 4), .22), #to account for text on bottom
  # rel_heights = c(2/9, 2/9, 1/9, 2/9, 2/9),
  ncol=1 
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Instantaneous Juvenile Mortality Rate", 
                                gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()

surv_fit_panel2 <- cowplot::plot_grid(
  cowplot::get_legend(leg_plot),
  surv_fit_panel,
  ncol=1, rel_heights=c(.075, .925)
)


# predicted gen length fits
gen_pred_list <- pmap(list(gen_dfa, gen_tbl$names, gen_tbl$years), 
                       fitted_preds,
                       descend_order = FALSE,
                      subset = 5)
# scale colors based on observed range over entire dataset
col_ramp_gen <- gen_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)

gen_fit <- map2(gen_pred_list, x_axes, .f = function(x, y) {
  plot_fitted_pred(x, print_x = y,
                   col_ramp = col_ramp_gen, facet_col = 5, 
                   col_ramp_direction = 1) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
})

gen_fit <- pmap(list(gen_pred_list, gen_tbl$group_labs, x_axes), 
                .f = function(x, y, z) {
                  plot_fitted_pred(x, print_x = z,
                                   col_ramp = col_ramp_gen,
                                   col_ramp_direction = 1
                                   # , facet_col = 7
                  ) +
                    scale_y_continuous(
                      name = y, position = 'right', sec.axis = dup_axis(),
                      labels = scales::number_format(accuracy = 0.1)
                    ) 
                })

gen_fit_panel <- cowplot::plot_grid(
  gen_fit[[1]], gen_fit[[2]], gen_fit[[3]], gen_fit[[4]], gen_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(rep(.195, 4), .22), #to account for text on bottom  
  # rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1
) %>% 
  arrangeGrob(., 
              left = textGrob("Centered Mean Age", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()

# make one version with a legend to use in panel fig
leg_plot_g <- plot_fitted_pred(gen_pred_list[[1]], col_ramp = col_ramp_surv, 
                             # facet_col = 5, 
                             col_ramp_direction = 1,
                             leg_name = "5-year Mean of Centered Mean Age") +
  theme(legend.position = "top")  

gen_fit_panel2 <- cowplot::plot_grid(
  cowplot::get_legend(leg_plot_g),
  gen_fit_panel,
  ncol=1, rel_heights=c(.075, .925)
)

# pdf(here::here("figs", "fits_both_vars.pdf"), height = 12, width = 8)
# surv_fit_panel2
# gen_fit_panel2
# dev.off()

png(here::here("figs", "ms_figs", "mortality_fit.png"), height = 7, width = 8, 
    res = 300, units = "in")
surv_fit_panel2
dev.off()

png(here::here("figs", "ms_figs", "age_fit.png"), height = 7, width = 8, 
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
rot_surv <- map(surv_dfa, rotate_trends)
surv_trends <- pmap(
  list(rot_surv, surv_tbl$years, group_labs), 
  .f = prep_trends
  ) %>% 
  bind_rows() %>% 
  mutate(var = "Juvenile Mortality Rate")

rot_gen <- map(gen_dfa, rotate_trends)
gen_trends <- pmap(
  list(rot_gen, gen_tbl$years, group_labs), 
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


# pdf(here::here("figs", "trends_both_vars.pdf"), height = 7, width = 4)
# plot_one_trend(trends %>% filter(trend == "Trend 1"))
# plot_one_trend(trends %>% filter(trend == "Trend 2"))
# dev.off()

png(here::here("figs", "ms_figs", "trend1.png"), height = 7, width = 4, 
    res = 300, units = "in")
plot_one_trend(trends %>% filter(trend == "Trend 1"))
dev.off()

png(here::here("figs", "ms_figs", "trend2.png"), height = 7, width = 4, 
    res = 300, units = "in")
plot_one_trend(trends %>% filter(trend == "Trend 2"))
dev.off()


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


# CHECK REGIMES ----------------------------------------------------------------

glimpse(rot_surv[[4]])
r1 <- rot_surv[[4]]
f1 <- find_regimes(r1$trends_mean[1, ], 
             sds = (r1$trends_upper - r1$trends_mean)[1, ] / 1.96)
f2 <- find_regimes(r1$trends_mean[2, ], 
             sds = (r1$trends_upper - r1$trends_mean)[2, ] / 1.96)
plot_regime_model(f1$best_model)
plot_regime_model(f2$best_model)