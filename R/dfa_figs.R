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
library(RColorBrewer)

raw_data <- readRDS(here::here("data/salmon_data/cwt_indicator_surv_clean.RDS"))

# import juvenile mortality data
surv_tbl <- readRDS(here::here("data", "survival_fits", "surv_tbl.RDS"))
surv_dfa <- purrr::map(surv_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_") 
  # f_name <- paste(y, "two-trend", "bayesdfa_c.RDS", sep = "_") 
  readRDS(here::here("data", "survival_fits", f_name))
})

# import mean gen length
gen_tbl <- readRDS(here::here("data", "generation_fits", "gen_tbl.RDS"))
gen_dfa <- purrr::map(gen_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "ar", "bayesdfa_c.RDS", sep = "_")
  # f_name <- paste(y, "one-trend", "ar", "bayesdfa_c.RDS", sep = "_")
  readRDS(here::here("data", "generation_fits", f_name))
})

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

mean_data <- raw_data %>% 
  group_by(j_group3b, year) %>% 
  summarize(mean_surv = mean(survival, na.rm = T),
            mean_age = mean(gen_length, na.rm = T))


# PREDICTED FITS ---------------------------------------------------------------


# plotting functions
source(here::here("R", "functions", "plotting_functions.R"))

#remove x_axes except for last plot
x_axes <- c(F, F, F, F, T)

# predicted survival fits
surv_pred_list <- pmap(list(surv_dfa, surv_names$names, surv_tbl$years), 
                  fitted_preds,
                  descend_order = FALSE, year1_last_mean = 2011)


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
                                    facet_col = 5,
                                    #where should vert line be drawn
                                    year1_last_mean = 2011
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


## predicted survival in real space 
real_surv_pred_list <- purrr::map(surv_pred_list, function (x) {
  # calculate estimated uncentered survival rate in real space 
  x %>% 
    left_join(., 
              raw_data %>% 
                select(stock, Time = year, survival),
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

real_surv_fit <- pmap(list(real_surv_pred_list, group_labs, x_axes), 
                 .f = function(x, y, z) {
                   plot_fitted_pred_real(x, 
                                         print_x = z,
                                         facet_col = 5,
                                         year1_last_mean = 2011
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


## predicted gen length fits
gen_pred_list <- pmap(list(gen_dfa, gen_names$names, gen_tbl$years), 
                       fitted_preds,
                       descend_order = FALSE, year1_last_mean = 2011)

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
                                   facet_col = 5, year1_last_mean = 2011
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
              left = textGrob("Centered Mean Age-At-Maturity", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()


## predicted uncentered generation length 
uncent_gen_pred_list <- purrr::map(gen_pred_list, function (x) {
  # calculate estimated uncentered survival rate in real space 
  x %>% 
    left_join(., 
              raw_data %>% 
                select(stock, Time = year, gen_length),
              by = c("stock", "Time")) %>% 
    group_by(ID) %>% 
    mutate(obs_mean_age = mean(gen_length, na.rm = T)) %>%
    ungroup() %>% 
    # uncenter predictions and calculate in real space
    mutate(
      obs_y = gen_length,
      mean = mean + obs_mean_age,
      lo = obs_mean_age + lo,
      hi = obs_mean_age + hi,
      uncent_last_mean = obs_mean_age + last_mean
    )
  })

uncent_gen_fit <- pmap(
  list(uncent_gen_pred_list, group_labs, x_axes), 
  .f = function(x, y, z) {
    plot_fitted_pred_uncent(x, print_x = z,
                            col_ramp = col_ramp_gen,
                            facet_col = 5, year1_last_mean = 2011
    ) +
      scale_y_continuous(
        name = y, position = 'right', sec.axis = dup_axis(),
        labels = scales::number_format(accuracy = 0.1)
      ) 
  }
)

uncent_gen_fit_panel <- cowplot::plot_grid(
  uncent_gen_fit[[1]], uncent_gen_fit[[2]], uncent_gen_fit[[3]], 
  uncent_gen_fit[[4]], uncent_gen_fit[[5]],
  axis = c("r"), align = "v", 
  rel_heights = c(2/11, 2/11, 3/11, 1/11, 3/11),
  ncol=1
) %>% 
  arrangeGrob(., 
              left = textGrob("Mean Age-At-Maturity", 
                              gp = gpar(col = "grey30", fontsize = 12),
                              rot = 90)) %>% 
  grid.arrange()


# export all figures
png(here::here("figs", "ms_figs", "survival_fit_ar1.png"), height = 11.5, 
    width = 10, res = 300, units = "in")
cowplot::plot_grid(surv_fit_panel)
dev.off()

png(here::here("figs", "ms_figs", "survival_fit_real_ar1.png"), height = 11.5, 
    width = 10, res = 300, units = "in")
cowplot::plot_grid(real_surv_fit_panel)
dev.off()

png(here::here("figs", "ms_figs", "age_fit_ar1.png"), height = 11.5, width = 10, 
    res = 300, units = "in")
cowplot::plot_grid(gen_fit_panel)
dev.off()

png(here::here("figs", "ms_figs", "age_fit_uncent_ar1.png"), height = 11.5,
    width = 10,  res = 300, units = "in")
cowplot::plot_grid(uncent_gen_fit_panel)
dev.off()


# CALCULATE SUMMARY STATISTICS -------------------------------------------------

# palette for following figs
tri_pal <- c("#08519c", "#a50f15", "grey60")
names(tri_pal) <- c("above", "below", "average")

## survival
surv_iters <- pmap(
  list(surv_tbl$group, surv_dfa, surv_tbl$names, surv_tbl$years), 
  function(group, x, y, z) {
    tt <- reshape2::melt(predicted(x), 
                         varnames = c("iter", "chain", "time", "stock")) %>% 
      left_join(., 
                data.frame(year = z,
                           time = unique(.$time)),
                by = "time") %>% 
      mutate(iter_chain = paste(iter, chain, sep = "_") %>% 
               as.factor %>% 
               as.numeric,
             group = group)
    tt$stock <- as.factor(y$stock[tt$stock])
    
    return(tt)
  }
) %>% 
  bind_rows()

mean_surv <- raw_data %>% 
  group_by(stock) %>%
  summarize(mean_logit_surv = mean(qlogis(survival), na.rm = T))
surv_iters2 <- left_join(surv_iters, mean_surv, by = "stock") %>% 
  mutate(uncent_value = value + mean_logit_surv,
         uncent_real_value = plogis(uncent_value))

min_surv_yrs <- raw_data %>% 
  filter(!is.na(survival)) %>% 
  group_by(j_group3b) %>% 
  summarize(window_yrs = min(year) + 4) %>% 
  rename(group = j_group3b)

roll_surv <- surv_iters2 %>% 
  group_by(iter_chain, year, group) %>% 
  summarize(mean_value_link = mean(value),
            mean_value = mean(uncent_real_value),
            .groups = "drop") %>% 
  group_by(iter_chain, group) %>% 
  arrange(iter_chain, year) %>% 
  mutate(
    roll_mean = slider::slide_dbl(mean_value, mean, .before = 5),
    roll_mean_link = slider::slide_dbl(mean_value_link, mean, .before = 5),
  ) %>% 
  group_by(year, group) %>% 
  summarize(
    mean_surv = mean(roll_mean),
    up = quantile(roll_mean, probs = 0.95),
    low = quantile(roll_mean, probs = 0.05),
    mean_surv_link = mean(roll_mean_link),
    up_link = quantile(roll_mean_link, probs = 0.95),
    low_link = quantile(roll_mean_link, probs = 0.05)
  ) %>% 
  left_join(., min_surv_yrs, by = "group") %>% 
  filter(year >= window_yrs) 
roll_surv$state <- case_when(
  roll_surv$low_link > 0 ~ "above",
  roll_surv$up_link < 0 ~ "below",
  TRUE ~ "average"
)
roll_surv$state <- factor(roll_surv$state, 
                          levels = c("above", "below", "average"))
roll_surv$group2 <- factor(
  roll_surv$group, 
  labels = c("North Yearling", "SoG Subyearling", "Puget Subyearling",
             "Puget Yearling", "South Subyearling")
)
 
roll_surv_link_ribbon <- ggplot(roll_surv) +
  geom_pointrange(aes(x = year, y = mean_surv_link, ymin = low_link, 
                      ymax = up_link, fill = state),
                  shape = 21) +
  scale_fill_manual(values = tri_pal, name = "") +
  facet_wrap(~group2) +
  labs(y = "Rolling Mean Logit Survival") +
  ggsidekick::theme_sleek()
roll_surv_ribbon <- ggplot(roll_surv) +
  geom_pointrange(aes(x = year, y = mean_surv, ymin = low, ymax = up, 
                      fill = state),
                  shape = 21) +
  scale_fill_manual(values = tri_pal, name = "") +
  facet_wrap(~group2) +
  labs(y = "Rolling Mean Survival") +
  ggsidekick::theme_sleek()


roll_surv_ppn <- surv_iters2 %>% 
  group_by(group) %>%
  mutate(n_stock = length(unique(stock)),
         window_year = min(year) + 4) %>%
  group_by(iter_chain, year, group) %>% 
  summarize(
    prob_below_0 = sum(value < 0) / n_stock,
    prob_above_0 = sum(value > 0) / n_stock,
    .groups = "drop"
  ) %>% 
  group_by(iter_chain, group) %>% 
  arrange(iter_chain, year) %>% 
  mutate(
    roll_ppn = slider::slide_dbl(prob_above_0, mean, .before = 5)
    ) %>% 
  group_by(year, group) %>% 
  summarize(
    mean_ppn = mean(roll_ppn),
    up = quantile(roll_ppn, probs = 0.95),
    low = quantile(roll_ppn, probs = 0.05),
    .groups = "drop"
  ) %>% 
  left_join(., min_surv_yrs, by = "group") %>% 
  filter(year >= window_yrs) 
roll_surv_ppn$state <- case_when(
      roll_surv_ppn$low > 0.5 ~ "above",
      roll_surv_ppn$up < 0.5 ~ "below",
      TRUE ~ "average"
    )
roll_surv_ppn$state <- factor(roll_surv_ppn$state, 
                              levels = c("above", "below", "average"))


## age
age_iters <- pmap(list(gen_tbl$group, gen_dfa, gen_tbl$names, gen_tbl$years), 
                  function(group, x, y, z) {
                    tt <- reshape2::melt(predicted(x), 
                                         varnames = c("iter", "chain", "time", "stock")) %>% 
                      left_join(., 
                                data.frame(year = z,
                                           time = unique(.$time)),
                                by = "time") %>% 
                      mutate(iter_chain = paste(iter, chain, sep = "_") %>% 
                               as.factor %>% 
                               as.numeric,
                             group = group)
                    tt$stock <- as.factor(y$stock[tt$stock])
                    
                    return(tt)
                  }) %>% 
  bind_rows()

mean_gen <- raw_data %>% 
  group_by(stock) %>%
  summarize(mean_age = mean(gen_length, na.rm = T))
age_iters2 <- left_join(age_iters, mean_gen, by = "stock") %>% 
  mutate(uncent_value = value + mean_age)

min_age_yrs <- raw_data %>% 
  filter(!is.na(gen_length)) %>% 
  group_by(j_group3b) %>% 
  summarize(window_yrs = min(year) + 4) %>% 
  rename(group = j_group3b)

roll_age <- age_iters2 %>% 
  group_by(iter_chain, year, group) %>% 
  summarize(mean_value_cent = mean(value),
            mean_value = mean(uncent_value),
            .groups = "drop") %>% 
  group_by(iter_chain, group) %>% 
  arrange(iter_chain, year) %>% 
  mutate(
    roll_mean = slider::slide_dbl(mean_value, mean, .before = 5),
    roll_mean_cent = slider::slide_dbl(mean_value_cent, mean, .before = 5),
  ) %>% 
  group_by(year, group) %>% 
  summarize(
    mean_age = mean(roll_mean),
    up = quantile(roll_mean, probs = 0.95),
    low = quantile(roll_mean, probs = 0.05),
    mean_age_cent = mean(roll_mean_cent),
    up_cent = quantile(roll_mean_cent, probs = 0.95),
    low_cent = quantile(roll_mean_cent, probs = 0.05)
  ) %>% 
  left_join(., min_age_yrs, by = "group") %>% 
  # remove incomplete years
  filter(year >= window_yrs) 

roll_age$state <- case_when(
  roll_age$low_cent > 0 ~ "above",
  roll_age$up_cent < 0 ~ "below",
  TRUE ~ "average"
)
roll_age$state <- factor(roll_age$state, levels = c("above", "below", "average"))
roll_age$group2 <- factor(
  roll_age$group, 
  labels = c("North Yearling", "SoG Subyearling", "Puget Subyearling",
             "Puget Yearling", "South Subyearling")
)

roll_age_cent_ribbon <- ggplot(roll_age) +
  geom_pointrange(aes(x = year, y = mean_age_cent, ymin = low_cent, 
                      ymax = up_cent, fill = state),
                  shape = 21) +
  scale_fill_manual(values = tri_pal) +
  facet_wrap(~group2) +
  labs(y = "Rolling Mean Centered Age") +
  ggsidekick::theme_sleek()
roll_age_ribbon <- ggplot(roll_age) +
  geom_pointrange(aes(x = year, y = mean_age, ymin = low, 
                      ymax = up, fill = state),
                  shape = 21) +
  scale_fill_manual(values = tri_pal, name = "") +
  facet_wrap(~group2) +
  labs(y = "Rolling Mean Age") +
  ggsidekick::theme_sleek()


roll_age_ppn <- age_iters2 %>% 
  group_by(group) %>%
  mutate(n_stock = length(unique(stock))) %>%
  group_by(iter_chain, year, group) %>% 
  summarize(
    prob_below_0 = sum(value < 0) / n_stock,
    prob_above_0 = sum(value > 0) / n_stock,
    .groups = "drop"
  ) %>% 
  group_by(iter_chain, group) %>% 
  arrange(iter_chain, year) %>% 
  mutate(
    roll_ppn = slider::slide_dbl(prob_above_0, mean, .before = 5, 
                                 complete = TRUE)
  ) %>% 
  group_by(year, group) %>% 
  summarize(
    mean_ppn = mean(roll_ppn),
    up = quantile(roll_ppn, probs = 0.95),
    low = quantile(roll_ppn, probs = 0.05)
  )  %>% 
  left_join(., min_age_yrs, by = "group") %>% 
  filter(year >= window_yrs) 

roll_age_ppn$state <- case_when(
  roll_age_ppn$low > 0.5 ~ "above",
  roll_age_ppn$up < 0.5 ~ "below",
  TRUE ~ "average"
)
roll_age_ppn$state <- factor(roll_age_ppn$state, 
                              levels = c("above", "below", "average"))


## combine survival and age proportion of stock estimates to look at last moving
# window estimate
roll_ppns <- rbind(roll_surv_ppn %>% mutate(dat = "survival"),
                   roll_age_ppn %>% mutate(dat = "mean age")) %>% 
  ungroup() %>% 
  filter(year == "2016") %>% 
  mutate(group = factor(group, labels = group_labs))


final_ppns <- ggplot(roll_ppns) +
  geom_pointrange(aes(x = group, y = mean_ppn, ymin = low, ymax = up, 
                      fill = dat), 
                  shape = 21, position = position_dodge(0.3)) +
  scale_fill_brewer(type = "div", palette = 1, name = "") +
  geom_hline(yintercept = 0.5, lty = 2) +
  labs(y = "Proportion of Stocks Above Long-Term Mean") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank())
  

# roll_surv_ppn_ribbon <- ggplot(roll_surv_ppn) +
#   geom_pointrange(aes(x = year, y = mean_ppn, ymin = low, ymax = up, 
#                       fill = state), shape = 21) +
#   scale_fill_manual(values = tri_pal) +
#   facet_wrap(~group) +
#   labs(y = "Rolling Mean Proportion of Stocks\nWith Above Average Survival") +
#   ggsidekick::theme_sleek()
# 
# roll_age_ppn_ribbon <- ggplot(roll_age_ppn) +
#   geom_pointrange(aes(x = year, y = mean_ppn, ymin = low, ymax = up, 
#                       fill = state), shape = 21) +
#   scale_fill_manual(values = tri_pal) +
#   facet_wrap(~group) +
#   labs(y = "Rolling Mean Proportion of Stocks\nWith Above Average Mean Age") +
#   ggsidekick::theme_sleek()


# pdf(here::here("figs", "supp_figs", "roll_window_surv_summary_stats_pts.pdf"))
# roll_surv_link_ribbon
# roll_surv_ribbon
# roll_surv_ppn_ribbon
# dev.off()
# 
# pdf(here::here("figs", "supp_figs", "roll_window_age_summary_stats_pts.pdf"))
# roll_age_cent_ribbon
# roll_age_ribbon
# roll_age_ppn_ribbon
# dev.off()

png(here::here("figs", "ms_figs", "surv_roll_pts.png"), height = 5, width = 7, 
    res = 300, units = "in")
roll_surv_ribbon
dev.off()

png(here::here("figs", "ms_figs", "age_roll_pts.png"), height = 5, width = 7, 
    res = 300, units = "in")
roll_age_ribbon
dev.off()


png(here::here("figs", "ms_figs", "final_ppns.png"), height = 4, width = 5.5, 
    res = 300, units = "in")
final_ppns
dev.off()


# slightly tweaked versions for presentations 
png(here::here("figs", "ms_figs", "surv_roll_pts_wide.png"), 
    height = 2, width = 8.5, res = 300, units = "in")
roll_surv_ribbon +
  facet_wrap(~group2, nrow = 1)
dev.off()

png(here::here("figs", "ms_figs", "age_roll_pts_wide.png"), 
    height = 2, width = 8.5, res = 300, units = "in")
roll_age_ribbon +
  facet_wrap(~group2, nrow = 1)
dev.off()



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
  bind_rows() %>% 
  mutate(model = "Survival Rate")
gen_pars <- map2(gen_dfa, group_labs, pull_par_f) %>% 
  bind_rows() %>% 
  mutate(model = "Age-At-Maturity")

phi_plot <- rbind(surv_pars, gen_pars) %>% 
  mutate(group = fct_relevel(as.factor(group), "SoG\nSubyearling", 
                             after = 1)) %>% 
  ggplot(., 
         aes(x = group, y = value, fill = trend)) + 
  scale_fill_brewer(name = "", palette = "Paired") +
  geom_violin(position = position_dodge(0.3),
              draw_quantiles = 0.5) + 
  geom_hline(aes(yintercept = y_int), lty = 2) + 
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  ylab("Posterior Probability Density") +
  xlab("Stock Grouping") +
  ggsidekick::theme_sleek() +
  guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
  facet_wrap(~model, scales = "free_x")


png(here::here("figs", "ms_figs", "phi_pars.png"), height = 7, width = 8, 
    res = 300, units = "in")
phi_plot
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
                             "SoG\nSubyearling", "Puget\nSubyearling", 
                             "Puget\nYearling","South\nSubyearling"),
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
         ) %>% 
           as.factor(),
         group = fct_relevel(as.factor(group), "North\nYearling", 
                             "SoG\nSubyearling", "Puget\nSubyearling", 
                             "Puget\nYearling","South\nSubyearling")) 

dum1 <- trends %>% 
  mutate(trend = as.numeric(trend),
         life_history = as.factor(life_history),
         data = "trend") %>% 
  select(median = x, lo, hi, trend, time, group, var, life_history,
         data)
dum2 <- regimes %>% 
  mutate(trend = as.numeric(as.factor(trend)),
         data = "regime",
         keep = ifelse(
           State == "State 2" & var == "Juvenile Mortality Rate" |
             State == "State 1" & var == "Mean Age",
           "yes",
           "no"
         )) %>%
  filter(keep == "yes") %>% 
  select(median, lo = lwr, hi = upr, trend, time, group, var, life_history,
         data)

dum3 <- rbind(dum1, dum2) %>% 
  mutate(
    facet_var = paste(trend, data, sep = "_") %>% 
      factor(., 
             levels = c("1_trend", "1_regime", "2_trend", "2_regime"),
             labels = c("Trend One", "Mean State One", "Trend Two", 
                        "Mean State Two")),
    color_var = paste(trend, life_history, sep = "_") %>% 
      factor(., 
             levels = c("1_yearling", "1_subyearling", "2_yearling",
                        "2_subyearling"))
  ) 

trend_pal <- c("#E08214", "#FDB863", "#8073AC", "#B2ABD2")
names(trend_pal) <- unique(dum3$color_var)

surv_trend_regime <- ggplot(
  data = dum3 %>% filter(var == "Juvenile Mortality Rate"), 
  aes(x = time, y = median)) +
  geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "color_var",
                         fill = "color_var", lty = "data"), 
              alpha = 0.4) + 
  geom_line(aes_string(colour = "color_var", lty = "data"), size = 1.2) + 
  scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
  facet_grid(facet_var~group, scales = "free_y") +
  ggsidekick::theme_sleek() + 
  labs(y = "Juvenile Survival Rate") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = unit(c(0.5,1,0.5,0.5), "cm") #t,r,b,l
    )+
  scale_colour_manual(values = trend_pal) +
  scale_fill_manual(values = trend_pal)


gen_trend_regime <- ggplot(
  data = dum3 %>% filter(var == "Mean Age"), 
  aes(x = time, y = median)
) +
  geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "color_var",
                         fill = "color_var"), 
              alpha = 0.4) + 
  geom_line(aes_string(colour = "color_var"), size = 1.2) + 
  scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
  facet_grid(facet_var~group, scales = "free_y") +
  ggsidekick::theme_sleek() + 
  labs(y = "Mean Age-At-Maturity") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = unit(c(0.5,1,0.5,0.5), "cm") #t,r,b,l
  ) +
  scale_colour_manual(values = trend_pal) +
  scale_fill_manual(values = trend_pal)

png(here::here("figs", "ms_figs", "surv_trend_regime.png"), 
    width = 8.5, height = 6, res = 300, units = "in")
surv_trend_regime
dev.off()

png(here::here("figs", "ms_figs", "age_trend_regime.png"), 
    width = 8.5, height = 6, res = 300, units = "in")
gen_trend_regime
dev.off()


# MEANS AFTER REGIMES ----------------------------------------------------------

surv_regime_years <- dum2 %>% 
  filter(var == "Juvenile Mortality Rate",
         trend == "1") %>%
  mutate(
    state = case_when(
      group %in% c("North\nYearling", "Puget\nSubyearling", 
                   "Puget\nYearling") ~ "stationary",
      median >= 0.5 ~ "high",
      median < 0.5 ~ "low"
    )
  ) %>% 
  select(year = time, group, state)
levels(surv_regime_years$group) <- levels(surv_iters2$group)

age_regime_years <- dum2 %>% 
  filter(var == "Mean Age",
         trend == "1") %>%
  mutate(
    state = case_when(
      group %in% c("North\nYearling", "Puget\nYearling") ~ "stationary",
      median >= 0.5 ~ "high",
      median < 0.5 ~ "low"
    )
  ) %>% 
  select(year = time, group, state)
levels(age_regime_years$group) <- levels(age_iters2$group)


## point intervals
reg_est <- left_join(surv_iters2, surv_regime_years, 
                     by = c("group", "year")) %>%
  group_by(group, state) %>% 
  summarize(
    mean_surv = mean(uncent_real_value),
    up = quantile(uncent_real_value, probs = 0.95),
    low = quantile(uncent_real_value, probs = 0.05)
  )
reg_est_age <- left_join(age_iters2, age_regime_years, 
                     by = c("group", "year")) %>%
  group_by(group, state) %>% 
  summarize(
    mean_age = mean(uncent_value),
    up = quantile(uncent_value, probs = 0.95),
    low = quantile(uncent_value, probs = 0.05)
  )


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
              bottom = textGrob("Mean Age-At-Maturity Model Loadings",
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
              bottom = textGrob("Survival Rate Model Loadings",
                                gp = gpar(col = "grey30", fontsize = 12)))

png(here::here("figs", "ms_figs", "surv_loadings.png"), height = 7, width = 10, 
    res = 300, units = "in")
grid.arrange(surv_load_panel)
dev.off()

png(here::here("figs", "ms_figs", "age_loadings.png"), height = 7, width = 10, 
    res = 300, units = "in")
grid.arrange(gen_load_panel)
dev.off()




