## Survival and covariate gams
# Exploratory survival analysis with full suite of environmental covariates
# Replaced by survival_gams
# Oct. 1 2020
# Updated Dec. 21 - split from survival_gams; needs to be cleaned and reran with 
# herring


library(mgcv)
library(tidyverse)

file_path <- "data/salmonData/"

#survival data
surv_raw <- read.csv(here::here(file_path, "cwt_indicator_surv_clean.csv")) 
surv <- surv_raw %>% 
  filter(
    j_group == "salish",
    !is.na(M)) %>%
  #add ocean entry year to match covariates based on life history
  mutate(year = ifelse(smolt == "streamtype", brood_year + 2, 
                       brood_year + 1)) %>% 
  select(year, stock, smolt, j_group, j_group2, survival, M) %>% 
  mutate(smolt = as.factor(smolt),
         stock = as.factor(stock),
         j_group = as.factor(j_group),
         j_group2 = as.factor(j_group2)) %>% 
  droplevels()

# merge covariate data then fit two GAMs (random intercepts and random splines)
cov <- read.csv(here::here("data/salmonData/survCovariateAnom.csv"),
                stringsAsFactors = F)

dat <- cov %>% 
  select(year, metric, env_anomaly = anomaly) %>% 
  nest(data = c(year, env_anomaly)) %>% 
  mutate(
    #for each covariate combine with survival data
    data = map(data, function (x) {
      surv %>% 
        left_join(., x, by = "year") %>% 
        filter(!is.na(env_anomaly))
    }),
    #fit gam w/ random intercepts
    gam1 = map(data, function (x) {
      gam(survival ~ s(env_anomaly, by = smolt, bs = "tp", m = 2, k = 4) +  
            s(stock, bs = "re"), 
          data = x, family=betar(link="logit"))
    }),
    #fit gam w/ stock-specific splines
    gam2 = map(data, function (x) {
      gam(survival ~ s(env_anomaly, m = 2, bs = "tp", k = 4) + 
            s(env_anomaly, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(env_anomaly, by = stock, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = x, family=betar(link="logit"))
    }),
    #save AIC for both
    aic = map2(gam1, gam2, function(x, y) {
      data.frame(
        model = c("rand_ints", "stock_splines"),
        aic = c(AIC(x), AIC(y))
      )
    })
  )

dat <- dat %>% 
  filter(!grepl("rkw", metric)) %>% 
  rbind(dat2)

saveRDS(dat, here::here("data", "gamFits", "salish_beta_gam.RDS"))

# generate fixed effects predictions from simpler model
dat$fixed_preds <- map2(dat$data, dat$gam1, function(dat_in, gam_in) {
  pred_seq <- seq(min(dat_in$env_anomaly, na.rm = T),
                  max(dat_in$env_anomaly, na.rm = T), 
                  length.out = 100)
  new_dat <- expand.grid(env_anomaly = pred_seq, smolt = unique(dat_in$smolt))
  preds <- predict(gam_in, new_dat, se.fit = TRUE, exclude="s(stock)",
                   newdata.guaranteed=TRUE, type = "link")
  new_dat %>%
    mutate(fit = boot::inv.logit(preds$fit),
           up = boot::inv.logit(preds$fit + (qnorm(0.975) * preds$se.fit)),
           low = boot::inv.logit(preds$fit + (qnorm(0.025) * preds$se.fit)),
           logit_fit = preds$fit,
           logit_up = preds$fit + (qnorm(0.975) * preds$se.fit),
           logit_low = preds$fit + (qnorm(0.025) * preds$se.fit)
    )
})

fe_splines <- dat %>% 
  select(metric, fixed_preds) %>% 
  unnest(fixed_preds) %>% 
  left_join(., cov %>% select(metric, class), by = "metric") %>%
  # glimpse()
  ggplot(.) +
  geom_line(aes(x = env_anomaly, y = fit, color = smolt)) +
  geom_ribbon(aes(x = env_anomaly, ymin = low, ymax = up, fill = smolt), 
              alpha = 0.3) +
  facet_wrap(~fct_reorder(metric, as.numeric(as.factor(class))), 
             scales = "free_x") +
  ggsidekick::theme_sleek() +
  labs(x = "Environmental Anomaly", y = "Predicted Survival")

pdf(here::here("figs", "gam", "fixed_effect_splines.pdf"), 
    height = 5, width = 7)
fe_splines
dev.off()

saveRDS(fe_splines, here::here("data", "gamFits", "fe_splines.rds"))

# FE predictions (as above) but in link space
pdf(here::here("figs", "gam", "fixed_effect_splines_link.pdf"), 
    height = 5, width = 7)
dat %>% 
  select(metric, fixed_preds) %>% 
  unnest(fixed_preds) %>% 
  left_join(., cov %>% select(metric, class), by = "metric") %>%
  ggplot(.) +
  geom_line(aes(x = env_anomaly, y = logit_fit, color = smolt)) +
  geom_ribbon(aes(x = env_anomaly, ymin = logit_low, ymax = logit_up, 
                  fill = smolt), alpha = 0.3) +
  facet_wrap(~fct_reorder(metric, as.numeric(as.factor(class))), 
             scales = "free_x") +
  ggsidekick::theme_sleek()
dev.off()


# generate stock-specific predictions from more complicated model
# observed data for plotting
dat$ss_preds <- map2(dat$data, dat$gam2, function(dat_in, gam_in) {
  pred_seq <- seq(min(dat_in$env_anomaly, na.rm = T),
                  max(dat_in$env_anomaly, na.rm = T), 
                  length.out = 100)
  new_dat <- expand.grid(env_anomaly = pred_seq, 
                         stock = unique(dat_in$stock)) %>% 
    left_join(., dat_in %>% select(stock, smolt) %>% distinct(), by = "stock")
  
  preds <- predict(gam_in, new_dat, se.fit = TRUE, type = "link")
  new_dat %>%
    mutate(
      stock = fct_reorder(stock, as.numeric(smolt)),
      fit = as.numeric(boot::inv.logit(preds$fit)),
      up = as.numeric(boot::inv.logit(preds$fit + 
                                        (qnorm(0.975) * preds$se.fit))),
      low = as.numeric(boot::inv.logit(preds$fit + 
                                         (qnorm(0.025) * preds$se.fit))),
      logit_fit = preds$fit,
      logit_up = preds$fit + (qnorm(0.975) * preds$se.fit),
      logit_low = preds$fit + (qnorm(0.025) * preds$se.fit)
    )
})

pdf(here::here("figs", "gam", "stock_specific_splines.pdf"), 
    height = 5, width = 7)
pmap(list(dat$ss_preds, dat$data, dat$metric), function (pred_dat, obs_dat, title) {
  ggplot() +
    geom_line(data = pred_dat, aes(x = env_anomaly, y = fit, color = smolt)) +
    geom_ribbon(data = pred_dat, aes(x = env_anomaly, ymin = low, ymax = up, 
                                     fill = smolt), 
                alpha = 0.3) +
    geom_point(data = obs_dat, aes(x = env_anomaly, y = survival, fill = smolt),
               shape = 21) +
    facet_wrap(~stock, scales = "free_y")  +
    ggsidekick::theme_sleek() +
    labs(title = title)
})
dev.off()

pdf(here::here("figs", "gam", "stock_specific_splines_link.pdf"), 
    height = 5, width = 7)
pmap(list(dat$ss_preds, dat$data, dat$metric), function (pred_dat, obs_dat, title) {
  ggplot() +
    geom_line(data = pred_dat, aes(x = env_anomaly, y = logit_fit, 
                                   color = smolt)) +
    geom_ribbon(data = pred_dat, aes(x = env_anomaly, ymin = logit_low, 
                                     ymax = logit_up, fill = smolt), 
                alpha = 0.3) +
    geom_point(data = obs_dat, aes(x = env_anomaly, y = boot::logit(survival),
                                   fill = smolt),
               shape = 21) +
    facet_wrap(~stock)  +
    ggsidekick::theme_sleek() +
    labs(title = title)
})
dev.off()


# ------------------------------------------------------------------------------

# Model selection with SST arc, seals and RKW 
trim_cov <- cov %>% 
  filter(metric %in% c("annual_anom_sst_arc", "rkw_anom", "seal_anom"),
         !is.na(anomaly),
         year > 1977, 
         year < 2013) %>% 
  group_by(metric) %>% 
  select(year, metric, anomaly) %>% 
  pivot_wider(names_from = "metric", values_from = "anomaly")

dat2 <- surv %>% 
  left_join(., trim_cov, by = "year") %>% 
  filter(!is.na(seal_anom))


# models
m1_sst <- gam(survival ~ s(annual_anom_sst_arc, m = 2, bs = "tp", k = 4) + 
                s(annual_anom_sst_arc, by = smolt, m = 1, bs = "tp", k = 4) + 
                s(stock, bs = "re"), 
              data = dat2, family=betar(link="logit"))
m1_seal <- gam(survival ~ s(seal_anom, m = 2, bs = "tp", k = 4) + 
                 s(seal_anom, by = smolt, m = 1, bs = "tp", k = 4) + 
                 s(seal_anom, by = stock, m = 1, bs = "tp", k = 4) + 
                 s(stock, bs = "re"), 
               data = dat2, family=betar(link="logit"))
m1_rkw <- gam(survival ~ s(rkw_anom, m = 2, bs = "tp", k = 4) + 
                s(rkw_anom, by = smolt, m = 1, bs = "tp", k = 4) + 
                s(rkw_anom, by = stock, m = 1, bs = "tp", k = 4) + 
                s(stock, bs = "re"), 
              data = dat2, family=betar(link="logit"))
m2_rkw <- gam(survival ~ s(annual_anom_sst_arc, m = 2, bs = "tp", k = 4) + 
                s(annual_anom_sst_arc, by = smolt, m = 1, bs = "tp", k = 4) + 
                s(rkw_anom, m = 2, bs = "tp", k = 4) + 
                s(rkw_anom, by = smolt, m = 1, bs = "tp", k = 4) + 
                s(rkw_anom, by = stock, m = 1, bs = "tp", k = 4) + 
                s(stock, bs = "re"), 
              data = dat2, family=betar(link="logit"))
m2_seal <- gam(survival ~ s(annual_anom_sst_arc, m = 2, bs = "tp", k = 4) + 
                 s(annual_anom_sst_arc, by = smolt, m = 1, bs = "tp", k = 4) +
                 s(seal_anom, m = 2, bs = "tp", k = 4) + 
                 s(seal_anom, by = smolt, m = 1, bs = "tp", k = 4) + 
                 s(seal_anom, by = stock, m = 1, bs = "tp", k = 4) + 
                 s(stock, bs = "re"), 
               data = dat2, family=betar(link="logit"))

AIC(m1_sst, m1_seal, m1_rkw, m2_rkw, m2_seal)


# Fit southern and Salish Sea --------------------------------------------------

surv2 <- surv_raw %>% 
  filter(
    !is.na(M),
    !group %in% c("north_streamtype", "north_oceantype")) %>% 
  #add ocean entry year to match covariates based on life history
  mutate(year = ifelse(smolt == "streamtype", brood_year + 2, 
                       brood_year + 1)) %>% 
  select(year, stock, smolt, group, survival, M) %>% 
  mutate(smolt = as.factor(smolt),
         stock = as.factor(stock),
         group = as.factor(group)) %>% 
  droplevels()

dat2 <- cov %>% 
  filter(!region == "sog") %>% 
  select(year, metric, env_anomaly = anomaly) %>% 
  nest(data = c(year, env_anomaly)) %>% 
  mutate(
    #for each covariate combine with survival data
    data = map(data, function (x) {
      surv2 %>% 
        left_join(., x, by = "year") %>% 
        filter(!is.na(env_anomaly))
    }),
    #fit gam w/ random intercepts
    gam1 = map(data, function (x) {
      gam(survival ~ s(env_anomaly, by = smolt, bs = "tp", m = 2, k = 4) +  
            s(stock, bs = "re"), 
          data = x, family=betar(link="logit"))
    }),
    #fit gam w/ stock-specific splines
    gam2 = map(data, function (x) {
      gam(survival ~ s(env_anomaly, m = 2, bs = "tp", k = 4) + 
            s(env_anomaly, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(env_anomaly, by = stock, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = x, family=betar(link="logit"))
    }),
    #save AIC for both
    aic = map2(gam1, gam2, function(x, y) {
      data.frame(
        model = c("rand_ints", "stock_splines"),
        aic = c(AIC(x), AIC(y))
      )
    })
  )


