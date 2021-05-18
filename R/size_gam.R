## Generation length gams
# Hierarchical analysis of size at release on juvenile survival and mean age
# May 17, 2021

library(mgcv)
library(tidyverse)

chin_dat <- readRDS(here::here("data", "salmon_data", 
                              "cwt_indicator_surv_clean.RDS")) %>% 
  filter(!is.na(avg_weight)) %>% 
  droplevels()

gen_dat <- chin_dat %>%
  filter(!is.na(gen_length)) %>% 
  droplevels()
surv_dat <- chin_dat %>%
  filter(!is.na(M)) %>% 
  droplevels()


# raw data
p <- ggplot(chin_dat) +
  facet_wrap(~fct_reorder(stock, as.numeric(j_group3b)),
             scales = "free_x") +
  ggsidekick::theme_sleek()
p +
  geom_point(aes(x = avg_weight, y = gen_length, fill = j_group3b), shape = 21)
p +
  geom_point(aes(x = avg_weight, y = M, fill = j_group3b), shape = 21)

ggplot(chin_dat) +
  geom_point(aes(x = avg_weight, y = gen_length, fill = smolt), shape = 21)


# helper functions
fit_dat_f <- function(dat, preds) {
  dat %>% 
    mutate(link_fit = as.numeric(preds$fit),
           link_se = as.numeric(preds$se.fit),
           link_lo = link_fit + (qnorm(0.025) * link_se),
           link_up = link_fit + (qnorm(0.975) * link_se),
           pred_gen = exp(link_fit),
           pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
           pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
    )
} 

gen_pred <- function(mod_in, dat_in, random = FALSE) {
  if (random == FALSE) {
    preds <- predict(mod_in, dat_in, se.fit = TRUE, 
                     exclude = smooth_pars[grepl("stock", smooth_pars)],
                     newdata.guaranteed = TRUE)
  } else {
    preds <- predict(mod_in, dat_in, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
  }
  
  #remove zero columns 
  dat_in %>% 
    mutate(link_fit = as.numeric(preds$fit),
           link_se = as.numeric(preds$se.fit),
           pred = exp(link_fit),
           pred_lo = exp(link_fit + (qnorm(0.025) * link_se)),
           pred_up = exp(link_fit + (qnorm(0.975) * link_se))
    ) 
}

  
# AGE MODELS -------------------------------------------------------------------

# quick test for best model structure
a_mod1 <- gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = chin_dat, method = "REML", family=Gamma(link="log"))
a_mod2 <- gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                smolt + s(stock, bs = "re"), 
              data = chin_dat, method = "REML", family=Gamma(link="log"))
a_mod3 <- gam(gen_length ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = chin_dat, method = "REML", family=Gamma(link="log"))
a_mod4 <- gam(gen_length ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = chin_dat, method = "REML", family=Gamma(link="log"))
AIC(a_mod1, a_mod2, a_mod3, a_mod4)


## Check fit
fit_preds <- predict(a_mod4, newdata = gen_dat, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
fit_dat <- fit_dat_f(gen_dat, preds = fit_preds) 

ggplot(data = fit_dat, aes(x = log(gen_length), y = link_fit)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(aes(color = brood_year), size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = link_lo, 
                     ymax = link_up),
                 size = 1/2, color = "firebrick4") +
  labs(x = "Observed Generation Length", 
       y = "Predicted Generation Length") +
  theme_bw()

# minimal evidence of temporal autocorrelation
acf(resid(a_mod4))


## Generate predictions
smooth_pars <- map(a_mod4$smooth, function(x) x$label) %>% unlist
weight_seq <- seq(min(gen_dat$avg_weight), max(gen_dat$avg_weight), 
                  length.out = 50)
fe_preds <- expand.grid(avg_weight = weight_seq,
                        smolt = unique(gen_dat$smolt)) %>% 
  gen_pred(a_mod4, ., random = FALSE)
re_preds <- expand.grid(avg_weight = weight_seq, 
                        stock = unique(gen_dat$stock)) %>% 
  left_join(., gen_dat %>% select(stock, smolt) %>% distinct(), by = "stock") %>% 
  gen_pred(a_mod4, ., random = TRUE)

fe_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  ggsidekick::theme_sleek() 

re_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  coord_cartesian(ylim = c(3, 6)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)



# MORTALITY RATE MODELS --------------------------------------------------------

# quick test for best model structure
s_mod1 <- gam(survival ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = surv_dat, method = "REML", family=Gamma(link="log"))
s_mod2 <- gam(survival ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                smolt + s(stock, bs = "re"), 
              data = surv_dat, method = "REML", family=Gamma(link="log"))
s_mod3 <- gam(survival ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = surv_dat, method = "REML", family=Gamma(link="log"))
s_mod4 <- gam(survival ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = surv_dat, method = "REML", family=Gamma(link="log"))
AIC(s_mod1, s_mod2, s_mod3, s_mod4)


## Check fit
fit_preds <- predict(s_mod4, newdata = surv_dat, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)

fit_dat_f(surv_dat, preds = fit_preds) %>% 
  ggplot(., aes(x = log(survival), y = link_fit)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(aes(color = brood_year), size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = link_lo, 
                     ymax = link_up),
                 size = 1/2, color = "firebrick4") +
  labs(x = "Observed Survival", 
       y = "Predicted Survival") +
  theme_bw()

# minimal evidence of temporal autocorrelation
acf(resid(s_mod4))


## Generate predictions
smooth_pars <- map(s_mod4$smooth, function(x) x$label) %>% unlist

weight_seq <- seq(min(surv_dat$avg_weight), max(surv_dat$avg_weight), 
                  length.out = 50)
fe_preds <- expand.grid(avg_weight = weight_seq,
                        smolt = unique(surv_dat$smolt)) %>% 
  gen_pred(s_mod4, ., random = FALSE)
re_preds <- expand.grid(avg_weight = weight_seq, 
                        stock = unique(surv_dat$stock)) %>% 
  left_join(., surv_dat %>% select(stock, smolt) %>% distinct(), by = "stock") %>% 
  gen_pred(s_mod4, ., random = TRUE)

fe_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  ggsidekick::theme_sleek() 

re_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  coord_cartesian(y = c(0, 0.15)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)
