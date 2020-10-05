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
  select(year, stock, smolt, survival, M) 

# # raw plot 
plot_dat <- surv %>%
  left_join(., cov, by = "year")
# 
# plot_dat %>% 
#   filter(smolt == "streamtype") %>% 
#   ggplot(., aes(x = anomaly, y = M)) +
#   geom_point() +
#   geom_smooth(method = "gam") + 
#   facet_wrap(~metric) 
  

# Test fit ---------------------------------------------------------------------

# Fit various GAM structures
temp <- plot_dat %>% 
  filter(metric == "seal_anom") %>% 
  mutate(smolt = as.factor(smolt),
         stock = as.factor(stock))
m1 <- gam(survival ~ s(anomaly, k = 4), data = temp, family=betar(link="logit"))
m2 <- gam(survival ~ s(anomaly, by = smolt, bs = "tp", m = 2, k = 4), 
          data = temp, family=betar(link="logit"))
m3 <- gam(survival ~ s(anomaly, m = 2, bs = "tp", k = 4) + 
            s(anomaly, by = smolt, m = 1, bs = "tp", k = 4), 
          data = temp, family=betar(link="logit"))
m4 <- gam(survival ~ s(anomaly, by = smolt, bs = "tp", m = 2, k = 4) +  
            s(stock, bs = "re"), 
          data = temp, family=betar(link="logit"))
m5 <- gam(survival ~ s(anomaly, m = 2, bs = "tp", k = 4) + 
            s(anomaly, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = temp, family=betar(link="logit"))
m6 <- gam(survival ~ s(anomaly, m = 2, bs = "tp", k = 4) + 
            s(anomaly, by = smolt, m = 1, bs = "tp", k = 4) + 
            s(anomaly, by = stock, m = 1, bs = "tp", k = 4) + 
            s(stock, bs = "re"), 
          data = temp, family=betar(link="logit"))
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
temp_list %>% 
  select(-gam_in) %>% 
  unnest(pred_dat) %>% 
  mutate(model_name = fct_relevel(as.factor(model_name), "group_sm", 
                                  "group_global_sm", "group_sm_re", 
                                  "group_global_sm_re")) %>% 
  ggplot(., aes(x = anomaly)) +
  geom_line(aes(y = fit, color = smolt)) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = smolt), alpha = 0.3) +
  facet_wrap(~model_name)

# focus on random effects models and look at stock-specific predictions
seal_seq <- seq(min(temp$anomaly, na.rm = T), 1, length.out = 100)
new_dat <- expand.grid(anomaly = seal_seq, stock = unique(temp$stock)) %>% 
  left_join(., temp %>% select(stock, smolt) %>% distinct(), by = "stock")
pred_dat <- predict(m5, new_dat, se.fit = TRUE, type = "response")
new_dat$fit <- as.numeric(pred_dat$fit)
new_dat$up <- as.numeric(new_dat$fit + (qnorm(0.975) * pred_dat$se.fit))
new_dat$low <- pmax(0, as.numeric(new_dat$fit + (qnorm(0.025) * pred_dat$se.fit)))

ggplot() +
  geom_line(data = new_dat, aes(x = anomaly, y = fit, color = smolt)) +
  geom_ribbon(data = new_dat, aes(x = anomaly, ymin = low, ymax = up, 
                                  fill = smolt), 
              alpha = 0.3) +
  geom_point(data = temp, aes(x = anomaly, y = survival, fill = smolt), 
             shape = 21) +
  facet_wrap(~fct_reorder(stock, as.numeric(smolt))) 

# ------------------------------------------------------------------------------



# merge covariate data
dat <- cov %>% 
  select(year, metric, env_anomaly = anomaly) %>% 
  nest(data = c(year, env_anomaly)) %>% 
  mutate(
    #for each covariate combine with survival data
    data = map(data, function (x) {
      surv %>% 
        left_join(., x, by = "year")
    }),
    raw_plot = map2(metric, data, function (name, data) {
      ggplot(data) +
        geom_point(aes(x = M, y = env_anomaly))
    })
  )

# plot data
dat$raw_data <- map2(dat)