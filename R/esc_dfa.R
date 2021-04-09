## Fit Bayes DFA models to escapement data
# Apr 8 2021
# Similar to m_dfa and gen_length_dfa, but does not include a ML model selection
# framework because groupings are so different. I.e. assumes juv3b is appropriate
# split


library(tidyverse)
library(bayesdfa)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)


esc <- read.csv(here::here("data", "salmonData", 
                           "clean_escapement_data.csv")) %>% 
  group_by(stock) %>% 
  mutate(
    j_group3b = as.factor(j_group3b),
    esc_z = as.numeric(scale(esc)),
    esc_cent = as.numeric(scale(esc, center = TRUE, scale = FALSE)),
    # convert back to thousands before calculating log
    log_esc = log(esc * 1000)
  ) %>% 
  ungroup()


# plot raw
esc %>%
  filter(!is.na(esc_z)) %>%
  ggplot(.) +
  geom_point(aes(x = year, y = log(escapement), fill = j_group3b), shape = 21) +
  facet_wrap(~ fct_reorder(stock, as.numeric(j_group3b))) +
  theme(legend.position = "top") +
  labs(y = "Centered Escapement") +
  ggsidekick::theme_sleek()

# check distribution
esc %>%
  group_by(stock) %>% 
  ggplot(.) +
  geom_histogram(aes(x = log(escapement), fill = j_group3b)) +
  facet_wrap(~ fct_reorder(stock, as.numeric(j_group3b))) +
  theme(legend.position = "top") +
  ggsidekick::theme_sleek()


## BAYES DFA -------------------------------------------------------------------


#helper function to spread and label input matrices for bayesdfa
make_mat <- function(x) {
  mat1 <- x %>%
    select(year, stock, log_esc) %>%
    spread(key = stock, value = log_esc) %>%
    as.matrix() 
  out_mat <- t(mat1[, 2:ncol(mat1)])
  colnames(out_mat) <- mat1[, "year"]
  return(out_mat)
}

# number of stocks per group
kept_grps <-  esc %>% 
  select(stock, j_group3b) %>% 
  distinct() %>%
  group_by(j_group3b) %>%
  tally() %>%
  filter(n > 2)

#generate tbl by group
esc_tbl <- tibble(group = levels(esc$j_group3b)) %>% 
  mutate(
    esc_mat = esc %>% 
      filter(!is.na(esc_z)) %>% 
      group_split(j_group3b) %>% 
      map(., make_mat)
  ) %>% 
  filter(group %in% kept_grps$j_group3b)
esc_tbl$names <- map(esc_tbl$esc_mat, function (x) {
  data.frame(stock_name = row.names(x)) 
})
esc_tbl$years <- map(esc_tbl$esc_mat, function (x) {
  as.numeric(colnames(x))
})


## fit DFA
furrr::future_map2(
  esc_tbl$esc_mat[1],
  esc_tbl$group[1],
  .f = function (y, group) {
    fit <- fit_dfa(
      y = y, num_trends = 2, zscore = TRUE,
      estimate_trend_ar = TRUE, #estimate_nu = TRUE,
      iter = 4000, chains = 4, thin = 1,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
    f_name <- paste(group, "two-trend", "bayesdfa_scaled.RDS", sep = "_")
    saveRDS(fit, here::here("data", "escapement_fits", f_name))
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)

# import
dfa_fits <- map(esc_tbl$group, function(y) {
  f_name <- paste(y, "two-trend", "bayesdfa_scaled.RDS", sep = "_") 
  readRDS(here::here("data", "escapement_fits", f_name))
}) 


## DIAGNOSTICS -----------------------------------------------------------------

# check diagnostics
map2(dfa_fits, esc_tbl$group, function (x, y) {
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


# chain plots for key pars
post <- map(dfa_fits, function (x) {
  as.array(x$samples)
})
np <- map(dfa_fits, function (x) {
  bayesplot::nuts_params(x$model)
})

all_pars <- dimnames(post[[1]])[3] %>% unlist() %>% as.character()
to_match <- c("phi", "xstar", "sigma")
pars <- all_pars[grepl(paste(to_match, collapse = "|"), all_pars)]

trace_list <- pmap(list(post, np, esc_tbl$group),
                   .f = function(x, y, z) {
                     bayesplot::mcmc_trace(x, pars = pars, np = y) +
                       labs(title = z)
                   })


pdf(here::here("figs", "dfa", "bayes", "escapement", "trace_plots.pdf"))
trace_list
dev.off()


## EXPORT ----------------------------------------------------------------------

# rotate trends and add to surv_tbl (keep DFA separate because they're huge)
esc_tbl$rot_esc <- map(dfa_fits, rotate_trends)

# test for evidence of regimes 
regime_f <- function(rots_in) {
  dum <- vector(nrow(rots_in$trends_mean), mode = "list")
  for(i in 1:nrow(rots_in$trends_mean)) {
    dum[[i]] <- find_regimes(
      rots_in$trends_mean[i, ], 
      sds = (rots_in$trends_upper - rots_in$trends_mean)[i, ] / 1.96,
      max_regimes = 3,
      iter = 3000,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
  }
  return(dum)
}

hmm_list_g[[1]] <- regime_f(esc_tbl$rot_gen[[1]])#furrr::future_map(esc_tbl$rot_gen, regime_f)
map(hmm_list_g, function(x) {
  map(x, ~.$table)
} )
#uniformly greatest support for two regimes


map(hmm_list_g, function(x) {
  map(x, function (y) plot_regime_model(y$best_model))
} )

esc_tbl$regime_trend1 <- map(hmm_list_g, function(x) x[[1]]$best_model)
esc_tbl$regime_trend2 <- map(hmm_list_g, function(x) x[[2]]$best_model)

saveRDS(esc_tbl, here::here("data", "escapement_fits", "esc_tbl.RDS"))



## TEMP PLOTS ------------------------------------------------------------------

# plotting functions
source(here::here("R", "functions", "plotting_functions.R"))

#fig labels
esc_tbl$group_labs <- c(
  "North\nYearling", 
  "South\nSubyearling",
  "North\nSubyearling",
  "SoG\nSubyearling", 
  "SoG\nYearling", 
  "Puget\nSubyearling", 
  "Puget\nYearling"
)

## Fits

# predicted survival fits
esc_pred_list <- pmap(list(dfa_fits, esc_tbl$names, esc_tbl$years), 
                       fitted_preds,
                       descend_order = FALSE)

# scale colors based on observed range over entire dataset
col_ramp_esc <- esc_pred_list %>% 
  bind_rows() %>% 
  pull(last_mean) %>% 
  range() %>% 
  abs() %>% 
  max() * c(-1, 1)


esc_fit <- pmap(list(esc_pred_list, esc_tbl$group_labs), 
                 .f = function(x, y) {
                   plot_fitted_pred(x, print_x = TRUE,
                                    col_ramp = col_ramp_esc,
                                    col_ramp_direction = 1
                   ) +
                     scale_y_continuous(
                       name = y, position = 'right', sec.axis = dup_axis(),
                       labels = scales::number_format(accuracy = 0.1)
                     ) +
                     theme(legend.position = "top")  
                 })

pdf(here::here("figs", "ms_figs", "esc_fits.pdf"))
esc_fit
dev.off()


## Trends
esc_trend_dat <- pmap(
  list(esc_tbl$rot_esc, esc_tbl$years, esc_tbl$group_labs), 
  .f = prep_trends
) %>% 
  bind_rows() %>% 
  mutate(var = "Log Escapement",
         trend = as.factor(trend), 
         life_history = case_when(
           grepl("Sub", group) ~ "subyearling",
           TRUE ~ "yearling"
         ),
         group = fct_relevel(as.factor(group), "North\nYearling",
                             "North\nSubyearling", "SoG\nYearling",
                             "SoG\nSubyearling", "Puget\nYearling",
                             "Puget\nSubyearling", "South\nSubyearling")
         )


esc_trends <- ggplot(esc_trend_dat, 
                     aes_string(x = "time", y = "x")) + 
  geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "life_history",
                         fill = "life_history"), 
              alpha = 0.4) + 
  geom_line(aes_string(colour = "life_history"), size = 1.2) + 
  scale_colour_brewer(type = "qual", name = "") +
  scale_fill_brewer(type = "qual", name = "") +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Return Year") + 
  ylab("Estimated Trend") +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(group ~ trend) +
  ggsidekick::theme_sleek() + 
  theme(legend.position = "top")


png(here::here("figs", "ms_figs", "esc_trends.png"), height = 8.5, width = 4, 
    res = 300, units = "in")
esc_trends
dev.off()


## Regimes
esc_regimes1 <- pmap(
  list(regime_model = esc_tbl$regime_trend1, years = esc_tbl$years, 
       group = esc_tbl$group_labs), 
  .f = prep_regime
) %>% 
  bind_rows() %>%
  mutate(trend = "One") 
esc_regimes2 <- pmap(
  list(regime_model = esc_tbl$regime_trend2, years = esc_tbl$years, 
       group = esc_tbl$group_labs), 
  .f = prep_regime
) %>% 
  bind_rows() %>% 
  mutate(trend = "Two")
esc_regime_dat <- rbind(esc_regimes1, esc_regimes2) %>% 
  mutate(var = "Log Escapement",
         state = as.factor(State), 
         life_history = case_when(
           grepl("Sub", group) ~ "subyearling",
           TRUE ~ "yearling"
         ),
         group = fct_relevel(as.factor(group), "North\nYearling",
                             "North\nSubyearling", "SoG\nYearling",
                             "SoG\nSubyearling", "Puget\nYearling",
                             "Puget\nSubyearling", "South\nSubyearling")) 
  

esc_regime <- ggplot(esc_regime_dat %>% filter(State == "State 1"), 
                     aes_string(x = "time", y = "median")) + 
  geom_ribbon(aes_string(ymin = "lwr", ymax = "upr", colour = "life_history",
                         fill = "life_history"), 
              alpha = 0.4, lty = 6) + 
  geom_line(aes_string(colour = "life_history"), size = 1.2, lty = 6) + 
  scale_colour_brewer(type = "qual", name = "") +
  scale_fill_brewer(type = "qual", name = "") +
  xlab("Return Year") + 
  ylab("Probability of Regime 1") +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(group ~ trend) +
  ggsidekick::theme_sleek() + 
  theme(legend.position = "top")

png(here::here("figs", "ms_figs", "esc_regime.png"), height = 8.5, width = 4, 
    res = 300, units = "in")
esc_regime
dev.off()


## Loadings

esc_load_dat <- pmap(list(esc_tbl$rot_esc, esc_tbl$names, esc_tbl$group),
                     .f = prep_loadings) 

#make list of figures
esc_load <- map2(esc_load_dat, esc_tbl$group_labs, plot_load, guides = FALSE) 

# make single figure to steal legend from
leg_esc_load <- plot_load(esc_load_dat[[1]], guides = TRUE) +
  ggsidekick::theme_sleek()

esc_load_panel <- 
  cowplot::plot_grid(
    esc_load[[1]], esc_load[[2]], esc_load[[3]], esc_load[[4]], esc_load[[5]],
    esc_load[[6]], esc_load[[7]],
    cowplot::get_legend(leg_esc_load),
    axis = c("lr"), align = "hv", 
    nrow = 2
  ) %>% 
  arrangeGrob(., 
              bottom = textGrob("Log Escapement Model Loadings",
                                gp = gpar(col = "grey30", fontsize = 12)))

png(here::here("figs", "ms_figs", "esc_loadings.png"), height = 10, width = 12, 
    res = 300, units = "in")
grid.arrange(esc_load_panel)
dev.off()
