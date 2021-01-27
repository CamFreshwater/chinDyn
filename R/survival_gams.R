## Survival and covariate gams
# Quick mixed effects model analysis of correlations between seal density 
# and Chinook survival to age 2 
# For now focused only on Salish Sea individuals
# Oct. 1 2020
# Update Dec. 18: remove RKW component (that's now analyzed with age data) and
# incorporate herring effects

library(mgcv)
library(tidyverse)

file_path <- "data/salmonData/"

#survival data
surv_raw <- read.csv(here::here(file_path, "cwt_indicator_surv_clean.csv")) 
surv <- surv_raw %>% 
  filter(
    # j_group == "salish",
    j_group2 == "sog",
    !is.na(M)) %>%
  #add ocean entry year to match covariates based on life history
  mutate(year = ifelse(smolt == "streamtype", brood_year + 2, 
                       brood_year + 1)) %>% 
  select(year, stock, smolt, j_group, survival, M) %>% 
  mutate(smolt = as.factor(smolt),
         stock = as.factor(stock),
         j_group = as.factor(j_group)) %>% 
  droplevels()

stock_key <- surv_raw %>% 
  filter(j_group == "salish") %>%
  select(stock, stock_name, region, smolt) %>% 
  distinct()

# covariate data
cov <- readRDS(here::here(file_path, "cov_subset_juv.rds"))

# raw plot 
dat <- surv %>%
  left_join(., cov, by = "year") %>% 
  filter(!is.na(seal_abund)) %>% 
  mutate(herr_z = scale(herr_abund)[ , 1],
         seal_z = scale(seal_abund)[ , 1],
         sst_z = scale(prime_entrance_sst)[ , 1],
         salinity_z = scale(prime_entrance_salinity)[ , 1])

# correlation among covariates
dat %>%
  filter(!is.na(sst_z)) %>% 
  select(year, herr_z, seal_z, sst_z, salinity_z) %>% 
  as.matrix() %>% 
  cor() %>% 
  corrplot::corrplot(., method = "number", "upper")


# a_palette <- disco::disco("bright", n = length(unique(dat$a_group)))
# names(a_palette) <- unique(by_dat$a_group)
# j_palette <- disco::disco("muted", n = length(unique(dat$j_group3)))
# names(j_palette) <- unique(by_dat$j_group3)


# MODEL COMPARISON -------------------------------------------------------------

# Compare predictions of herring and seal model using structure identified in 
# following section
n_stks <- length(unique(dat$stock))
s_mod <- gam(survival ~ s(seal_z, m = 2, bs = "tp", k = 4) + 
               s(seal_z, by = stock, m = 1, bs = "tp", k = 4) + 
               s(stock, bs = "re"), 
             data = dat, family = betar(link="logit"), method = "REML")
h_mod <- gam(survival ~ s(herr_z, m = 2, bs = "tp", k = 4) + 
               s(herr_z, by = stock, m = 1, bs = "tp", k = 4) + 
               s(stock, bs = "re"), 
             data = dat, family = betar(link="logit"), method = "REML")
sst_mod <- gam(survival ~ s(sst_z, m = 2, bs = "tp", k = 4) + 
                 s(sst_z, by = stock, m = 1, bs = "tp", k = 4) + 
                 s(stock, bs = "re"), 
               data = dat, family = betar(link="logit"), method = "REML")
sal_mod <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
                 s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) + 
                 s(stock, bs = "re"), 
               data = dat, family = betar(link="logit"), method = "REML")
h_s_mod <- gam(survival ~ s(seal_z, m = 2, bs = "tp", k = 4) + 
                 s(seal_z, by = stock, m = 1, bs = "tp", k = 4) +
                 s(herr_z, m = 2, bs = "tp", k = 4) + 
                 s(herr_z, by = stock, m = 1, bs = "tp", k = 4) + 
                 s(stock, bs = "re"), 
               data = dat, family = betar(link="logit"), method = "REML")
h_s_int_mod <- gam(survival ~ te(herr_z, seal_z, 
                                 bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                     t2(herr_z, seal_z, stock, m = 1,
                        bs = c("tp", "tp", "re"), 
                        k = c(3, 3, n_stks)) + 
                     s(stock, bs = "re"), 
                   data = dat, family = betar(link="logit"), method = "REML")
h_sst_mod <- gam(survival ~ s(sst_z, m = 2, bs = "tp", k = 4) + 
                   s(sst_z, by = stock, m = 1, bs = "tp", k = 4) +
                   s(herr_z, m = 2, bs = "tp", k = 4) + 
                   s(herr_z, by = stock, m = 1, bs = "tp", k = 4) + 
                   s(stock, bs = "re"), 
                 data = dat, family = betar(link="logit"), method = "REML")
h_sal_mod <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
                   s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) +
                   s(herr_z, m = 2, bs = "tp", k = 4) + 
                   s(herr_z, by = stock, m = 1, bs = "tp", k = 4) + 
                   s(stock, bs = "re"), 
                 data = dat, family = betar(link="logit"), method = "REML")
s_sst_mod <- gam(survival ~ s(sst_z, m = 2, bs = "tp", k = 4) + 
                   s(sst_z, by = stock, m = 1, bs = "tp", k = 4) +
                   s(seal_z, m = 2, bs = "tp", k = 4) + 
                   s(seal_z, by = stock, m = 1, bs = "tp", k = 4) + 
                   s(stock, bs = "re"), 
                 data = dat, family = betar(link="logit"), method = "REML")
s_sal_mod <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
                   s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) +
                   s(seal_z, m = 2, bs = "tp", k = 4) + 
                   s(seal_z, by = stock, m = 1, bs = "tp", k = 4) + 
                   s(stock, bs = "re"), 
                 data = dat, family = betar(link="logit"), method = "REML")
sal_sst_mod <- gam(survival ~ s(sst_z, m = 2, bs = "tp", k = 4) + 
                     s(sst_z, by = stock, m = 1, bs = "tp", k = 4) +
                     s(salinity_z, m = 2, bs = "tp", k = 4) + 
                     s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) + 
                     s(stock, bs = "re"), 
                   data = dat, family = betar(link="logit"), method = "REML")
h_s_sst_mod <- gam(survival ~  s(herr_z, m = 2, bs = "tp", k = 4) + 
                     s(herr_z, by = stock, m = 1, bs = "tp", k = 4) +
                     s(sst_z, m = 2, bs = "tp", k = 4) + 
                     # s(sst_z, by = stock, m = 1, bs = "tp", k = 4) +
                     s(seal_z, m = 2, bs = "tp", k = 4) + 
                     s(seal_z, by = stock, m = 1, bs = "tp", k = 4) + 
                     s(stock, bs = "re"), 
                   data = dat, family = betar(link="logit"), method = "REML")
h_s_sal_mod <- gam(survival ~  s(herr_z, m = 2, bs = "tp", k = 4) + 
                     s(herr_z, by = stock, m = 1, bs = "tp", k = 4) +
                     s(salinity_z, m = 2, bs = "tp", k = 4) + 
                     s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) +
                     s(seal_z, m = 2, bs = "tp", k = 4) + 
                     s(seal_z, by = stock, m = 1, bs = "tp", k = 4) + 
                     s(stock, bs = "re"), 
                   data = dat, family = betar(link="logit"), method = "REML")
AIC(s_mod, h_mod, sst_mod, sal_mod, h_s_mod, h_sst_mod, h_sal_mod, s_sst_mod, 
    s_sal_mod, h_s_sst_mod, h_s_sal_mod, h_s_int_mod) %>% 
  arrange(-AIC)


# TOP MODEL PREDICTIONS --------------------------------------------------------

## Equivalent to h_s_mod (top ranked)
mod <- gam(survival ~ s(seal_z, m = 2, bs = "tp", k = 4) + 
             s(seal_z, by = stock, m = 1, bs = "tp", k = 4) +
             s(herr_z, m = 2, bs = "tp", k = 4) + 
             s(herr_z, by = stock, m = 1, bs = "tp", k = 4) + 
             s(stock, bs = "re"), 
           data = dat, family = betar(link="logit"), method = "REML")


## Check fit (note that fit is on scale of observations so unlikely to violate
# assumption in link space)
fit_preds <- predict(mod, newdata = dat, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
fit_dat <- dat %>% 
  mutate(link_fit = as.numeric(fit_preds$fit),
         link_se = as.numeric(fit_preds$se.fit),
         pred_surv = plogis(link_fit),
         pred_surv_lo = plogis(link_fit + (qnorm(0.025) * link_se)),
         pred_surv_up = plogis(link_fit + (qnorm(0.975) * link_se))
  )

ggplot(data = fit_dat, aes(x = survival, y = pred_surv)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(aes(color = year), size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = pred_surv_lo, 
                     ymax = pred_surv_up),
                 size = 1/2, color = "firebrick4") +
  labs(x = "Observed Survival", 
       y = "Predicted Survival") +
  theme_bw()

# no evidence of temporal autocorrelation
acf(resid(mod))


## Generate counterfactual predictions for additive model
excl_pars <- map(mod$smooth, function(x) x$label) %>% unlist
# sst_seq <- seq(min(dat$sst_z, na.rm = T), max(dat$sst_z, na.rm = T), 
#                length.out = 75)
# sal_seq <- seq(min(dat$salinity_z, na.rm = T), max(dat$salinity_z, na.rm = T), 
#                length.out = 75)
seal_seq <- seq(-2.5, 0.75, length.out = 75)
herr_seq <- seq(-1.7, 1.8, length.out = 75)
new_dat <- expand.grid(
  seal_z = seal_seq,
  herr_z = herr_seq
  # sst_z = sst_seq,
  # salinity_z = sal_seq
  )

# fixed effects predictions
preds <- predict(mod, new_dat, se.fit = TRUE, 
                 exclude = excl_pars[grepl("stock", excl_pars)],
                 newdata.guaranteed = TRUE)
new_dat2 <- new_dat %>% 
  mutate(link_fit = as.numeric(preds$fit),
         link_se = as.numeric(preds$se.fit),
         pred_surv = plogis(link_fit),
         pred_surv_lo = plogis(link_fit + (qnorm(0.025) * link_se)),
         pred_surv_up = plogis(link_fit + (qnorm(0.975) * link_se))
         )

# zero_sal <- new_dat2$salinity_z[which.min(abs(new_dat2$salinity_z - 0))]
# new_dat2 %>% 
#   filter(salinity_z == zero_sal) %>% 
#   ggplot(., aes(x = sst_z)) +
#   geom_line(aes(y = pred_surv)) +
#   geom_ribbon(aes(ymin = pred_surv_lo, ymax = pred_surv_up), alpha = 0.3) +
#   ggsidekick::theme_sleek()
# zero_sst <- new_dat2$sst_z[which.min(abs(new_dat2$sst_z - 0))]
# new_dat2 %>% 
#   filter(sst_z == zero_sst) %>% 
#   ggplot(., aes(x = salinity_z)) +
#   geom_line(aes(y = pred_surv)) +
#   geom_ribbon(aes(ymin = pred_surv_lo, ymax = pred_surv_up), alpha = 0.3) +
#   ggsidekick::theme_sleek()


heat_plot <- ggplot(new_dat2, aes(x = seal_z, y = herr_z)) +
  geom_raster(aes(fill = pred_surv)) +
  scale_fill_viridis_c(name = "Predicted\nSurvival") +
  labs(x = "Harbour Seal Abundance", y = "Age-0 Herring Abundance") +
  ggsidekick::theme_sleek()

# herring only fixed effects
zero_seal <- new_dat2$seal_z[which.min(abs(new_dat2$seal_z - 0))]
h_plot <- new_dat2 %>% 
  filter(seal_z == zero_seal) %>% 
  ggplot(., aes(x = herr_z)) +
  geom_line(aes(y = pred_surv)) +
  geom_ribbon(aes(ymin = pred_surv_lo, ymax = pred_surv_up), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  labs(x = "Age-0 Herring Abundance", y = "Predicted Survival")

# seal only fixed effects
zero_herr <- new_dat2$herr_z[which.min(abs(new_dat2$herr_z - 0))]
s_plot <- new_dat2 %>% 
  filter(herr_z == zero_herr) %>% 
  ggplot(., aes(x = seal_z)) +
  geom_line(aes(y = pred_surv)) +
  geom_ribbon(aes(ymin = pred_surv_lo, ymax = pred_surv_up), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  labs(x = "Harbour Seal Abundance", y = "Predicted Survival")


# Mixed effects predictions
# include stock-specific effects 
new_dat_stock <- expand.grid(seal_z = seal_seq, 
                             herr_z = herr_seq,
                             stock = unique(dat$stock))

preds <- predict(mod, new_dat_stock, se.fit = TRUE, 
                 newdata.guaranteed = TRUE)
new_dat_stock2 <- new_dat_stock %>% 
  mutate(link_fit = as.numeric(preds$fit),
         link_se = as.numeric(preds$se.fit),
         pred_surv = plogis(link_fit),
         pred_surv_lo = plogis(link_fit + (qnorm(0.025) * link_se)),
         pred_surv_up = plogis(link_fit + (qnorm(0.975) * link_se))
  )

# herring only stock-specific effects
herr_rand <- new_dat_stock2 %>% 
  filter(seal_z == zero_seal) %>% 
  ggplot(., aes(x = herr_z)) +
  geom_line(aes(y = pred_surv)) +
  geom_ribbon(aes(ymin = pred_surv_lo, ymax = pred_surv_up), alpha = 0.3) +
  geom_point(data = dat, aes(x = herr_z, y = survival)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock) +
  labs(x = "Age-0 Herring Abundance", y = "Predicted Age-2 Survival")

# seal only stock-specific effects
seal_rand <- new_dat_stock2 %>% 
  filter(herr_z == zero_herr) %>% 
  ggplot(., aes(x = seal_z)) +
  geom_line(aes(y = pred_surv)) +
  geom_ribbon(aes(ymin = pred_surv_lo, ymax = pred_surv_up), alpha = 0.3) +
  geom_point(data = dat, aes(x = seal_z, y = survival)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock) +
  labs(x = "Harbour Seal Abundance", y = "Predicted Age-2 Survival")


png(here::here("figs", "gam", "surv_gam", "herr_seal_fe_heatplot.png"), 
    height = 5, width = 6, res = 300, units = "in")
heat_plot
dev.off()

png(here::here("figs", "gam", "surv_gam", "herr_seal_fe_splines.png"), 
    height = 3, width = 5.5, res = 300, units = "in")
ggpubr::ggarrange(h_plot, s_plot, nrow = 1, ncol = 2)
dev.off()

png(here::here("figs", "gam", "surv_gam", "herr_RE.png"), 
    height = 4.25, width = 6.5, res = 300, units = "in")
herr_rand
dev.off()

png(here::here("figs", "gam", "surv_gam", "seal_RE.png"), 
    height = 4.25, width = 6.5, res = 300, units = "in")
seal_rand
dev.off()


# TEST STRUCTURE ---------------------------------------------------------------

## Fit seal (or herring) response with various model structures to isolate how
# stock-specific and life-history effects should be accounted for
# TAKE HOME: Global smooth, stock-specific intercepts and smooths, but no life
# history parameter are the most parsimonious

# Fit various GAM structures
m1 <- gam(survival ~ s(salinity_z, k = 4), data = dat, 
          family=betar(link="logit"), method = "REML")
m2 <- gam(survival ~ s(salinity_z, by = smolt, bs = "tp", m = 2, k = 4), 
          data = dat, family=betar(link="logit"), method = "REML")
m3 <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
            s(salinity_z, by = smolt, m = 1, bs = "tp", k = 4), 
          data = dat, family=betar(link="logit"), method = "REML")
m4 <- gam(survival ~ s(salinity_z, by = smolt, bs = "tp", m = 2, k = 4) +  
            s(stock, bs = "re"), 
          data = dat, family=betar(link="logit"), method = "REML")
m5 <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
            s(salinity_z, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = dat, family=betar(link="logit"), method = "REML")
m6 <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
            s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = dat, family=betar(link="logit"), method = "REML")
m7 <- gam(survival ~ s(salinity_z, m = 2, bs = "tp", k = 4) + 
            s(salinity_z, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(salinity_z, by = stock, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = dat, family=betar(link="logit"), method = "REML")
AIC(m1, m2, m3, m4, m5, m6, m7)

#vector of coefficients to exclude from main effects predictions
exclude_coefs <- map(m6$smooth, function(x) x$label) %>% unlist

# Compare predictions among model structures
temp_list <- tibble(model_name = c("group_sm", "group_global_sm", "group_sm_re",
                                   "group_global_sm_re", 
                                   "group_global_sm_re_splines"),
                    gam_in = list(m2, m3, m4, m5, m7)
                    ) %>%
  #add predictions
  mutate(pred_dat = map2(model_name, gam_in, function(model_name, gam_in) {
    seal_seq <- seq(-2.5, 0.75, length.out = 100)
    herr_seq <- seq(-1.7, 1.8, length.out = 100)
    new_dat <- expand.grid(seal_anom = seal_seq, 
                           # herr_anom = herr_seq,
                           smolt = unique(dat$smolt))
    # exclude REs from predictions if present
    if (model_name %in% c("group_sm_re", "group_global_sm_re")) {
      preds <- predict(gam_in, new_dat, se.fit = TRUE, 
                       exclude = "s(stock)",
                       newdata.guaranteed=TRUE, type = "response")
    } else if (model_name == "group_global_sm_re_splines") {
      preds <- predict(gam_in, new_dat, se.fit = TRUE, 
                       exclude = exclude_coefs[grepl("stock", exclude_coefs)],
                       newdata.guaranteed = TRUE, type = "response")
    } else {
      preds <- predict(gam_in, new_dat, se.fit = TRUE, type = "response")
    }
    new_dat$fit <- as.numeric(preds$fit)
    new_dat$up <- as.numeric(new_dat$fit + (1.96 * preds$se.fit))
    new_dat$low <- as.numeric(new_dat$fit - (1.96 * preds$se.fit))
    return(new_dat)
  })) 

#unnest and plot
temp_list %>% 
  select(-gam_in) %>% 
  unnest(pred_dat) %>% 
  mutate(model_name = fct_relevel(as.factor(model_name), "group_sm", 
                                  "group_global_sm", "group_sm_re", 
                                  "group_global_sm_re")) %>% 
  ggplot(., aes(x = seal_anom)) +
  # ggplot(., aes(x = herr_anom)) +
  geom_line(aes(y = fit, color = smolt)) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = smolt), alpha = 0.3) +
  facet_wrap(~model_name) +
  ggsidekick::theme_sleek()


# focus on random effects models and look at stock-specific predictions
seal_seq <- seq(min(temp$anomaly, na.rm = T), 1, length.out = 100)
new_dat <- expand.grid(anomaly = seal_seq, stock = unique(temp$stock)) %>% 
  left_join(., temp %>% select(stock, smolt) %>% distinct(), by = "stock")
pred_dat <- predict(m6, new_dat, se.fit = TRUE, type = "response")
new_dat$fit <- as.numeric(pred_dat$fit)
new_dat$up <- as.numeric(new_dat$fit + (qnorm(0.975) * pred_dat$se.fit))
new_dat$low <- pmax(0, as.numeric(new_dat$fit + (qnorm(0.025) * pred_dat$se.fit)))

m6_preds <- ggplot() +
  geom_line(data = new_dat, aes(x = anomaly, y = fit, color = smolt)) +
  geom_ribbon(data = new_dat, aes(x = anomaly, ymin = low, ymax = up, 
                                  fill = smolt), 
              alpha = 0.3) +
  geom_point(data = temp, aes(x = anomaly, y = survival, fill = smolt), 
             shape = 21) +
  facet_wrap(~fct_reorder(stock, as.numeric(smolt)))  +
  ggsidekick::theme_sleek()

pdf(here::here("figs", "gam", "beta_seal_gam_splineRE.pdf"), 
    height = 5, width = 7)
fixed_effs_only
m6_preds
dev.off()


