## Survival and covariate gams
# Quick mixed effects model analysis of correlations between seal density 
# and Chinook survival to age 2 
# For now focused only on Salish Sea individuals
# Oct. 1 2020

library(mgcv)
library(tidyverse)

cov <- read.csv(here::here("data/salmonData/survCovariateAnom.csv"), 
                stringsAsFactors = F)
surv <- read.csv(here::here("data/salmonData/cwt_indicator_surv_clean.csv"), 
                 stringsAsFactors = FALSE) %>% 
  filter(agg_reg == "SS",
         !is.na(M)) %>% 
  select(year, stock, smolt, survival, M) %>% 
  mutate(smolt = as.factor(smolt),
         stock = as.factor(stock))

# # raw plot 
plot_dat <- surv %>%
  left_join(., cov, by = "year")
# 
plot_dat %>%
  filter(smolt == "streamtype") %>%
  ggplot(., aes(x = anomaly, y = M)) +
  geom_point() +
  geom_smooth(method = "gam") +
  facet_wrap(~metric)
  

# Test fit ---------------------------------------------------------------------

# Fit various GAM structures
seal_dat <- plot_dat %>% 
  filter(metric == "seal_anom") 
m1 <- gam(survival ~ s(anomaly, k = 4), data = seal_dat, 
          family=betar(link="logit"))
m2 <- gam(survival ~ s(anomaly, by = smolt, bs = "tp", m = 2, k = 4), 
          data = seal_dat, family=betar(link="logit"))
m3 <- gam(survival ~ s(anomaly, m = 2, bs = "tp", k = 4) + 
            s(anomaly, by = smolt, m = 1, bs = "tp", k = 4), 
          data = seal_dat, family=betar(link="logit"))
m4 <- gam(survival ~ s(anomaly, by = smolt, bs = "tp", m = 2, k = 4) +  
            s(stock, bs = "re"), 
          data = seal_dat, family=betar(link="logit"))
m5 <- gam(survival ~ s(anomaly, m = 2, bs = "tp", k = 4) + 
            s(anomaly, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = seal_dat, family=betar(link="logit"))
m6 <- gam(survival ~ s(anomaly, m = 2, bs = "tp", k = 4) + 
            s(anomaly, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(anomaly, by = stock, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = seal_dat, family=betar(link="logit"))
AIC(m1, m2, m3, m4, m5, m6)

#vector of coefficients to exclude from main effects predictions
tt <- map(m6$smooth, function(x) x$label) %>% unlist

# Compare predictions among model structures
temp_list <- tibble(model_name = c("group_sm", "group_global_sm", "group_sm_re",
                                   "group_global_sm_re", 
                                   "group_global_sm_re_splines"),
                    gam_in = list(m2, m3, m4, m5, m6)
                    ) %>%
  #add predictions
  mutate(pred_dat = map2(model_name, gam_in, function(model_name, gam_in) {
    seal_seq <- seq(-1.5, 1, length.out = 100)
    new_dat <- expand.grid(anomaly = seal_seq, smolt = unique(temp$smolt))
    # exclude REs from predictions if present
    if (model_name %in% c("group_sm_re", "group_global_sm_re")) {
      preds <- predict(gam_in, new_dat, se.fit = TRUE, exclude="s(stock)",
                       newdata.guaranteed=TRUE, type = "response")
    } else if (model_name == "group_global_sm_re_splines") {
      preds <- predict(gam_in, new_dat, se.fit = TRUE, exclude = tt[4:31],
                       newdata.guaranteed=TRUE, type = "response")
    } else {
      preds <- predict(gam_in, new_dat, se.fit = TRUE, type = "response")
    }
    new_dat$fit <- as.numeric(preds$fit)
    new_dat$up <- as.numeric(new_dat$fit + (1.96 * preds$se.fit))
    new_dat$low <- as.numeric(new_dat$fit - (1.96 * preds$se.fit))
    return(new_dat)
  })) 

#unnest and plot
fixed_effs_only <- temp_list %>% 
  select(-gam_in) %>% 
  unnest(pred_dat) %>% 
  mutate(model_name = fct_relevel(as.factor(model_name), "group_sm", 
                                  "group_global_sm", "group_sm_re", 
                                  "group_global_sm_re")) %>% 
  ggplot(., aes(x = anomaly)) +
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

# ------------------------------------------------------------------------------

# merge covariate data then fit two GAMs (random intercepts and random splines)
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
saveRDS(dat, here::here("data", "gamFits", "salish_beta_gam.RDS"))

# check AIC scores
dat %>% 
  select(metric, aic) %>% 
  unnest(aic) %>% 
  arrange(aic) %>% 
  print(n = Inf)

# generate fixed effects predictions from simpler model
dat$fixed_preds <- map2(dat$data, dat$gam1, function(dat_in, gam_in) {
  pred_seq <- seq(min(dat_in$env_anomaly, na.rm = T),
                  max(dat_in$env_anomaly, na.rm = T), 
                  length.out = 100)
  new_dat <- expand.grid(env_anomaly = pred_seq, smolt = unique(dat_in$smolt))
  preds <- predict(gam_in, new_dat, se.fit = TRUE, exclude="s(stock)",
                   newdata.guaranteed=TRUE, type = "response")
  new_dat %>%
    mutate(fit = as.numeric(preds$fit),
           up = as.numeric(fit + (1.96 * preds$se.fit)),
           low = as.numeric(fit - (1.96 * preds$se.fit))
           )
})

pdf(here::here("figs", "gam", "fixed_effect_splines.pdf"), 
    height = 5, width = 7)
dat %>% 
  select(metric, fixed_preds) %>% 
  unnest(fixed_preds) %>% 
  ggplot(.) +
  geom_line(aes(x = env_anomaly, y = fit, color = smolt)) +
  geom_ribbon(aes(x = env_anomaly, ymin = low, ymax = up, 
                                   fill = smolt), alpha = 0.3) +
  facet_wrap(~metric, scales = "free_x") +
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
  
  preds <- predict(gam_in, new_dat, se.fit = TRUE, type = "response")
  new_dat %>%
    mutate(fit = as.numeric(preds$fit),
           up = as.numeric(fit + (1.96 * preds$se.fit)),
           low = pmax(0, as.numeric(fit - (1.96 * preds$se.fit)))
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
    facet_wrap(~fct_reorder(stock, as.numeric(smolt)))  +
    ggsidekick::theme_sleek() +
    labs(title = title)
})
dev.off()