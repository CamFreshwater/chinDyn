## Supplementary analyses
# All data originates with PST Chinook technical committee and was provided by
# C. Parken

library(tidyverse)

dat <- readRDS(here::here("data/salmon_data/cwt_indicator_surv_clean.RDS")) %>% 
  mutate(logit_surv = qlogis(survival)) %>% 
  group_by(stock) %>% 
  mutate(
    age_z = scale(gen_length, scale = FALSE, center = TRUE) %>% as.numeric,
    logit_surv_z = scale(logit_surv, scale = FALSE, center = TRUE) %>% as.numeric
  ) %>% 
  ungroup()


# RAW BOX PLOTS ----------------------------------------------------------------

raw_data <- dat %>% 
  mutate(
    group = factor(
      j_group3b, 
      levels = c("north_oceantype", "north_streamtype", "sog_oceantype",
                 "sog_streamtype", "puget_oceantype", "puget_streamtype",    
                 "south_oceantype", "south_streamtype"),
      labels = c("North\nSubyearling", "North\nYearling", "SoG\nSubyearling", 
                 "SoG\nYearling", "Puget\nSubyearling", "Puget\nYearling", 
                 "South\nSubyearling",  "South\nYearling")
    ),
    smolt = factor(smolt, labels = c("Subyearling", "Yearling")),
    region_juv = factor(j_group3, levels = c("north", "sog", "puget", "south"),
                        labels = c("North", "SoG", "Puget", "South"))
  ) %>% 
  filter(!stock %in% c("TST", "AKS")) 

scale_1 <- c("#8073AC", "#E08214")
names(scale_1) <- levels(raw_data$smolt)


png(here::here("figs", "supp_figs", "survival_box.png"), height = 7, 
    width = 7.5, res = 300, units = "in")
ggplot() +
  geom_boxplot(data = raw_data, aes(x = as.factor(year), y = survival,
                                    fill = smolt)) +
  facet_grid(smolt ~ region_juv,
             scales = "free_y") +
  labs(x = "Ocean Entry Year", y = "Juvenile Marine Survival Rate") +
  scale_x_discrete(
    breaks = seq(1970, 2020, by = 10)
  ) +
  scale_fill_manual(values = scale_1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
dev.off()

png(here::here("figs", "supp_figs", "age_box.png"), height = 7, 
    width = 7.5, res = 300, units = "in")
ggplot() +
  geom_boxplot(data = raw_data, aes(x = as.factor(year), y = gen_length,
                                    fill = smolt)) +
  facet_grid(smolt ~ region_juv,
             scales = "free_y") +
  labs(x = "Ocean Entry Year", y = "Mean Age-At-Maturity") +
  scale_x_discrete(
    breaks = seq(1970, 2020, by = 10)
  ) +
  scale_fill_manual(values = scale_1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
dev.off()

png(here::here("figs", "supp_figs", "size_box.png"), height = 7, 
    width = 7.5, res = 300, units = "in")
ggplot(raw_data %>% filter(!is.na(avg_weight))) +
  geom_boxplot(aes(x = as.factor(year), y = avg_weight,
                                    fill = smolt)) +
  facet_grid(smolt ~ region_juv,
             scales = "free_y") +
  labs(x = "Ocean Entry Year", y = "Release Size") +
  scale_x_discrete(
    breaks = seq(1970, 2020, by = 10)
  ) +
  scale_fill_manual(values = scale_1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
dev.off()



## COVARIANCE ------------------------------------------------------------------

# function to visualize raw estimates
plot_foo <- function(x_in = "gen_length", y_in = "survival",
                     xlab = "Mean Age-At-Maturity", ylab = "Survival") {
  dat %>% 
    filter(!is.na(.data[[x_in]]),
           !is.na(.data[[y_in]])) %>% 
    droplevels() %>% 
    mutate(stock = fct_reorder(stock, as.numeric(j_group3b))) %>% 
    ggplot(.) +
    geom_point(aes_string(x = x_in, y = y_in, fill = "j_group3b"), shape = 21) +
    scale_fill_brewer(type = "qual", palette = 2, name = "Juvenile\nGrouping") +
    facet_wrap(~smolt, scales = "free") +
    labs(x = xlab, y = ylab) +
    ggsidekick::theme_sleek() 
}

pdf(here::here("figs", "supp_figs", "covs.pdf"))
plot_foo(x_in = "logit_surv", y_in = "gen_length")
plot_foo(x_in = "avg_weight", y_in = "survival",
         xlab = "Mean Mass", ylab = "Survival")
plot_foo(x_in = "avg_weight", y_in = "gen_length",
         xlab = "Mean Mass", ylab = "Mean Age-At-Maturity")
dev.off()


# fit gams to estimate effects
age_surv_mod <- mgcv::gam(gen_length ~ s(survival, m = 2, bs = "tp", k = 5) + 
                            s(survival, m = 1, bs = "tp", k = 5, by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat)

# wt_surv_mod <- mgcv::gam(logit_surv ~ s(avg_weight, m = 2, bs = "tp") + 
#                             s(avg_weight, m = 1, bs = "tp", by = j_group3b) +
#                             s(stock, bs = "re"),
#                           data = dat#,
#                           # family = mgcv::betar(link="logit")
#                          )
wt_surv_mod <- mgcv::gam(survival ~ s(avg_weight, m = 2, bs = "tp") + 
                           s(avg_weight, m = 1, bs = "tp", by = j_group3b) +
                           s(stock, bs = "re"),
                         data = dat,
                         family = mgcv::betar(link="logit")
)

age_wt_mod <- mgcv::gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 5) + 
                            s(avg_weight, m = 1, bs = "tp", k = 5, 
                              by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat#,
                          # family = Gamma(link="log")
                          )

# generate predictive dataframes
# cov_range <- dat %>% 
#   pivot_longer(cols = c(avg_weight, survival),
#                names_to = "cov") %>% 
#   group_by(j_group3b, cov) %>% 
#   mutate(min_val = min(value, na.rm = T),
#          max_val = max(value, na.rm = T)) %>%
#   select(j_group3b, cov, min_val, max_val) %>%
#   distinct() %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   mutate(
#     new_data = purrr::pmap(list(cov, min_val, max_val), 
#                            function (cov, x, y) {
#       cov = seq(from = x, to = y, length = 100)
#     })
#   ) 
# 
# surv_pred <- cov_range %>% 
#   filter(cov == "survival") %>% 
#   unnest(cols = c(new_data)) %>%
#   select(smolt, j_group3b, 
#          survival = new_data) %>% 
#   mutate(stock = unique(dat$stock)[3]) 
# 
# as_preds <- predict(age_surv_mod, newdata = surv_pred, se.fit = T,
#                     exclude = "s(stock)")
# # as_preds <- predict.lm(age_surv_lm_f, newdata = surv_pred, se.fit = T) 
# 
# surv_pred$fit <- as.numeric(as_preds$fit)
# surv_pred$se_fit <- as.numeric(as_preds$se.fit)
# 
# ggplot(surv_pred) +
#   geom_line(aes(x = survival, y = fit)) +
#   facet_wrap(~j_group3b)
# 
# 
# wt_pred <- cov_range %>% 
#   filter(cov == "avg_weight") %>% 
#   unnest(cols = c(new_data)) %>%
#   select(j_group3b, avg_weight = new_data) %>% 
#   mutate(stock = unique(dat$stock)[3]) 
# 
# sw_preds <- predict(wt_surv_mod2, newdata = wt_pred, se.fit = T,
#                     exclude = "s(stock)", type = "link")
# 
# wt_pred2 <- wt_pred %>% 
#   mutate(
#     fit = as.numeric(sw_preds$fit),
#     se_fit = as.numeric(as_preds$se.fit),
#     surv = plogis(fit),
#     surv_lo = plogis(fit + (qnorm(0.025) * se_fit)),
#     surv_up = plogis(fit + (qnorm(0.975) * se_fit))
#   )
# 
# ggplot(wt_pred2) +
#   geom_line(aes(x = avg_weight, y = surv)) +
#   geom_ribbon(aes(x = avg_weight, ymin = surv_lo, ymax = surv_up, alpha = 0.3)) +
#   facet_wrap(~j_group3b)
# 
# 
# aw_preds <- predict(age_wt_mod, newdata = wt_pred, se.fit = T,
#                     exclude = "s(stock)")
# 
# wt_pred$fit <- as.numeric(aw_preds$fit)
# wt_pred$se_fit <- as.numeric(aw_preds$se.fit)
# 
# ggplot(wt_pred) +
#   geom_line(aes(x = avg_weight, y = fit)) +
#   # geom_ribbon(aes(x = avg_weight, ymin = surv_lo, ymax = surv_up, alpha = 0.3)) +
#   facet_wrap(~j_group3b)


## time series of residuals
yr_preds <- expand.grid(
  year = seq(min(dat$year), max(dat$year), by = 1),
  j_group3b = unique(dat$j_group3b)
)


# age-survival residuals
as_dat <- dat %>% 
  filter(!is.na(gen_length), 
         !is.na(survival)) %>% 
  mutate(
    resids = resid(age_surv_mod),
  ) %>% 
  pivot_longer(
    cols = c(age_z, resids), names_to = "data_type"
  )
as_dat$data_type <- factor(as_dat$data_type, 
                           labels = c("observations", "residuals"))


as_obs <-  gam(value ~ s(year, m = 2, bs = "tp", k = 5) + 
                             s(year, m = 1, bs = "tp", k = 5, by = j_group3b), 
                           data = as_dat %>% filter(data_type == "observations"))
as_resid <-  gam(value ~ s(year, m = 2, bs = "tp", k = 5) + 
                       s(year, m = 1, bs = "tp", k = 5, by = j_group3b),
                     data = as_dat %>% filter(data_type == "residuals"))

as_obs_pred <- predict(as_obs, newdata = yr_preds, se.fit = T)
as_resid_pred <- predict(as_resid, newdata = yr_preds, se.fit = T)

yr_preds$as_obs <- as_obs_pred$fit %>% as.numeric
yr_preds$as_resid <- as_resid_pred$fit %>% as.numeric

yr_preds %>% 
  pivot_longer(cols = c(as_obs, as_resid), names_to = "data_type") %>%
  ggplot(., aes(x = year, y = value)) +
  geom_jitter(data = as_dat, aes(fill = data_type), shape = 21, alpha = 0.3) +
  geom_line(aes(colour = data_type), size = 1.5) +
  facet_wrap(~j_group3b) +
  ggsidekick::theme_sleek() +
  labs(x = "Release Year", y = "Centered Mean Age-At-Maturity")


# weight-survival 
ws_dat <- dat %>% 
  filter(!is.na(avg_weight), 
         !is.na(logit_surv)) %>% 
  mutate(
    resids = resid(wt_surv_mod),
  ) %>% 
  pivot_longer(
    cols = c(logit_surv_z, resids), names_to = "data_type"
  )
ws_dat$data_type <- factor(ws_dat$data_type, 
                           labels = c("observations", "residuals"))


ws_obs <-  gam(value ~ s(year, m = 2, bs = "tp", k = 5) + 
                 s(year, m = 1, bs = "tp", k = 5, by = j_group3b), 
               data = ws_dat %>% filter(data_type == "observations"))
ws_resid <-  gam(value ~ s(year, m = 2, bs = "tp", k = 5) + 
                   s(year, m = 1, bs = "tp", k = 5, by = j_group3b),
                 data = ws_dat %>% filter(data_type == "residuals"))

ws_obs_pred <- predict(ws_obs, newdata = yr_preds, se.fit = T)
ws_resid_pred <- predict(ws_resid, newdata = yr_preds, se.fit = T)

yr_preds$ws_obs <- ws_obs_pred$fit %>% as.numeric
yr_preds$ws_resid <- ws_resid_pred$fit %>% as.numeric

yr_preds %>% 
  pivot_longer(cols = c(ws_obs, ws_resid), names_to = "pred_type") %>%
  ggplot(., aes(x = year, y = value)) +
  geom_jitter(data = ws_dat, aes(fill = data_type), shape = 21, alpha = 0.3) +
  geom_line(aes(colour = pred_type), size = 2) +
  facet_wrap(~j_group3b) +
  ggsidekick::theme_sleek() +
  labs(x = "Release Year", y = "Centered Logit Juvenile Survival")


# age-weight residuals
aw_dat <- dat %>% 
  filter(!is.na(avg_weight), 
         !is.na(gen_length)) %>% 
  mutate(
    resids = resid(age_wt_mod),
  ) %>% 
  pivot_longer(
    cols = c(age_z, resids), names_to = "data_type"
  )
aw_dat$data_type <- factor(aw_dat$data_type, 
                           labels = c("observations", "residuals"))

aw_obs <-  gam(value ~ s(year, m = 2, bs = "tp", k = 5) + 
                 s(year, m = 1, bs = "tp", k = 5, by = j_group3b), 
               data = aw_dat %>% filter(data_type == "observations"))
aw_resid <-  gam(value ~ s(year, m = 2, bs = "tp", k = 5) + 
                   s(year, m = 1, bs = "tp", k = 5, by = j_group3b),
                 data = aw_dat %>% filter(data_type == "residuals"))

aw_obs_pred <- predict(aw_obs, newdata = yr_preds, se.fit = T)
aw_resid_pred <- predict(aw_resid, newdata = yr_preds, se.fit = T)

yr_preds$age_z <- aw_obs_pred$fit %>% as.numeric
yr_preds$resids <- aw_resid_pred$fit %>% as.numeric

yr_preds %>% 
  pivot_longer(cols = c(aw_obs, aw_resid), names_to = "pred_type") %>%
  ggplot(., aes(x = year, y = value)) +
  geom_jitter(data = aw_dat, aes(fill = data_type), shape = 21, alpha = 0.3) +
  geom_line(aes(colour = pred_type), size = 1.5) +
  facet_wrap(~j_group3b) +
  ggsidekick::theme_sleek() +
  labs(x = "Release Year", y = "Centered Mean Age-At-Maturity")



## SYNCHRONY -------------------------------------------------------------------

dat %>% 
  filter(!is.na(survival)) %>% 
  group_by(j_group3b, year) %>% 
  summarize(n_stocks = length(unique(stock)), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_line(aes(x = year, y = n_stocks, colour = j_group3b))
  glimpse()
