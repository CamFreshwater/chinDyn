## Bayesian DFA Figures
# March 15, 2021
# Use juvenile mortality and mean generation length DFAs to generate manuscript
# figures

library(tidyverse)
library(bayesdfa)

# plotting functions
source(here::here("R", "functions", "plot_fitted_bayes.R"))

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


# PREDICTED FITS ---------------------------------------------------------------

# predicted survival fits
surv_pred_list <- pmap(list(surv_dfa, surv_tbl$names, surv_tbl$years), 
                  fitted_preds,
                  descend_order = TRUE)
# scale colors based on observed range over entire dataset
col_ramp_in <- surv_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range()

surv_fit <- map(surv_pred_list, plot_fitted_pred, 
                col_ramp = col_ramp_in, facet_col = 5)

surv_fit_panel <- cowplot::plot_grid(
  surv_fit[[1]], surv_fit[[2]], surv_fit[[3]], surv_fit[[4]], surv_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1 
)

# histogram of last_means
last_means <- map2(surv_pred_list, surv_tbl$group, function (x, y) {
  x %>% 
    select(ID, last_mean) %>% 
    distinct() %>% 
    mutate(group = y)
}) %>% 
  bind_rows() %>% 
  glimpse()
hist(last_means$last_mean)
nrow(last_means[which(last_means$last_mean > 0), ]) / nrow(last_means)


# predicted gen length fits
gen_pred_list <- pmap(list(gen_dfa, gen_tbl$names, gen_tbl$years), 
                       fitted_preds,
                       descend_order = FALSE)
# scale colors based on observed range over entire dataset
col_ramp_gen <- gen_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range()

gen_fit <- map(gen_pred_list, plot_fitted_pred, 
                col_ramp = col_ramp_gen, facet_col = 5, col_ramp_direction = 1)

gen_fit_panel <- cowplot::plot_grid(
  gen_fit[[1]], gen_fit[[2]], gen_fit[[3]], gen_fit[[4]], gen_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 3/11, 1/11, 2/11, 3/11),
  ncol=1 
)

pdf(here::here("figs", "fits_both_vars.pdf"), height = 12, width = 8)
surv_fit_panel
gen_fit_panel
dev.off()


# CALCULATE FINAL FIVE YEAR MEANS ----------------------------------------------

# generation length means in last five years
gen_prob <- map2(gen_dfa, gen_tbl$names, final_prob, n_years = 5) %>% 
  bind_rows() %>% 
  mutate(
    thresh = ifelse(prob_below_0 > 0.90, 1, 0)
  )

surv_prob <- map2(surv_dfa, surv_tbl$names, final_prob, n_years = 5) %>% 
  bind_rows() %>% 
  mutate(
    thresh = ifelse(prob_above_0 > 0.90, 1, 0)
  )


# ESTIMATED TRENDS -------------------------------------------------------------

# prep dataframes for each
rot_surv <- map(surv_dfa, rotate_trends)
surv_trends <- pmap(
  list(rot_surv, surv_tbl$years, surv_tbl$group), 
  .f = prep_trends
  ) %>% 
  bind_rows() %>% 
  mutate(var = "Juvenile M")

rot_gen <- map(gen_dfa, rotate_trends)
gen_trends <- pmap(
  list(rot_gen, gen_tbl$years, gen_tbl$group), 
  .f = prep_trends
) %>% 
  bind_rows() %>% 
  mutate(var = "Generation Length")

trends <- rbind(surv_trends, gen_trends) %>% 
  mutate(var = as.factor(var),
         trend = as.factor(trend), 
         life_history = case_when(
           grepl("stream", group) ~ "yearling",
           grepl("ocean", group) ~ "subyearling"
         ),
         group = fct_relevel(as.factor(group), "north_streamtype", 
                             "sog_oceantype", "puget_streamtype", 
                             "puget_oceantype", "south_oceantype"),
         var = fct_relevel(as.factor(var), "Juvenile M", "Generation Length"))


pdf(here::here("figs", "trends_both_vars.pdf"), height = 7, width = 4)
plot_one_trend(trends %>% filter(trend == "Trend 1"))
plot_one_trend(trends %>% filter(trend == "Trend 2"))
dev.off()



# ESTIMATED LOADINGS -----------------------------------------------------------

gen_loadings <- map2(rot_gen, gen_tbl$names, prep_loadings) %>% 
  bind_rows() %>% 
  glimpse()



ggplot(tt, aes_string(x = "name", y = "value", fill = "trend", 
                     alpha = "prob_diff0")) + 
  scale_alpha_continuous(name = "Probability\nDifferent") +
  scale_fill_discrete(name = "") +
  geom_violin(color = NA, position = position_dodge(0.3)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  coord_flip() + 
  xlab("Time Series") + 
  ylab("Loading") +
  lims(y = c(-0.75, 0.75)) +
  ggsidekick::theme_sleek() +
  guides(alpha = guide_legend(override.aes = list(fill = "grey")))

loadings_list <- pmap(list(rot_list, gen_tbl$names, gen_tbl$group), 
                      function(x, names, title) {
                        plot_loadings(x, names = names$stock_name) +
                          lims(y = c(-1, 1)) +
                          labs(title = title) +
                          theme(axis.text.y = element_text(angle = 45, 
                                                           vjust = -1))
                      }
) 