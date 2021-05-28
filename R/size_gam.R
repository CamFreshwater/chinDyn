## Generation length gams
# Hierarchical analysis of size at release on juvenile survival and mean age
# May 17, 2021

library(mgcv)
library(tidyverse)

chin_dat <- readRDS(here::here("data", "salmon_data", 
                              "cwt_indicator_surv_clean.RDS")) %>% 
  filter(!is.na(avg_weight)) %>% 
  droplevels() %>% 
  #scale within a stock 
  group_by(stock) %>% 
  mutate(avg_weight_z = as.numeric(scale(avg_weight))) %>% 
  ungroup() 

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
  geom_point(aes(x = avg_weight_z, y = gen_length, fill = j_group3b), shape = 21)
ggplot(chin_dat) +
  geom_point(aes(x = avg_weight_z, y = M, fill = j_group3b), shape = 21)


# helper functions
fit_dat_f <- function(dat, preds, resp_dist = c("gamma", "beta")) {
  dum <- dat %>% 
    mutate(link_fit = as.numeric(preds$fit),
           link_se = as.numeric(preds$se.fit),
           link_lo = link_fit + (qnorm(0.025) * link_se),
           link_up = link_fit + (qnorm(0.975) * link_se))
  
  if (resp_dist == "gamma") {
    out <- dum %>%  
      mutate(
        pred = exp(link_fit),
        pred_lo = exp(link_fit + (qnorm(0.025) * link_se)),
        pred_up = exp(link_fit + (qnorm(0.975) * link_se))
      )
  }
  if (resp_dist == "beta") {
    out <- dum %>%  
      mutate(
        pred = plogis(link_fit),
        pred_lo = plogis(link_fit + (qnorm(0.025) * link_se)),
        pred_up = plogis(link_fit + (qnorm(0.975) * link_se))
      )
  }
  return(out)
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

# quick test for best model structure (strong support for including 
# autocorrelation pars)
a_mod1 <- gamm(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3),
               random = list(stock = ~ 1),
               data = gen_dat, method = "REML", family=Gamma(link="log"),
               correlation = corAR1(form = ~ brood_year | stock))
a_mod2 <- gamm(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                smolt, 
               random = list(stock = ~ 1),
              data = gen_dat, method = "REML", family=Gamma(link="log"),
              correlation = corAR1(form = ~ brood_year | stock))
a_mod3 <- gamm(gen_length ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                smolt,
               random = list(stock = ~ 1),
              data = gen_dat, method = "REML", family=Gamma(link="log"),
              correlation = corAR1(form = ~ brood_year | stock))
a_mod4 <- gamm(gen_length ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                smolt,
               random = list(stock = ~ 1 + avg_weight),
              data = gen_dat, method = "REML", family=Gamma(link="log"),
              correlation = corAR1(form = ~ brood_year | stock))
AIC(a_mod1$lme, a_mod2$lme, a_mod3$lme, a_mod4$lme)

acf(resid(a_mod2$lme, type = "normalized"))
hist(resid(a_mod2$lme, type = "normalized"))


# second equivalent model with scaled data -- i.e. how does observed variability
# within a stock impact survival
a_mod2b <- gamm(gen_length ~ s(avg_weight_z, m = 2, bs = "tp", k = 3) +
                 smolt, 
               random = list(stock = ~ 1),
               data = gen_dat, method = "REML", family=Gamma(link="log"),
               correlation = corAR1(form = ~ brood_year | stock))


## Check fit
fit_preds <- predict(a_mod2$gam, newdata = gen_dat, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
fit_dat <- fit_dat_f(gen_dat, preds = fit_preds, resp_dist = "gamma") 

ggplot(data = fit_dat, aes(x = log(gen_length), y = link_fit)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(aes(color = brood_year), size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = link_lo, 
                     ymax = link_up),
                 size = 1/2, color = "firebrick4") +
  labs(x = "Observed Generation Length", 
       y = "Predicted Generation Length") +
  theme_bw()


## Generate predictions
smooth_pars <- map(a_mod2$gam$smooth, function(x) x$label) %>% unlist
weight_seq <- seq(min(gen_dat$avg_weight), max(gen_dat$avg_weight), 
                  length.out = 50)
fe_preds <- expand.grid(avg_weight = weight_seq,
                        smolt = unique(gen_dat$smolt)) %>% 
  gen_pred(a_mod2$gam, ., random = FALSE)
re_preds <- expand.grid(avg_weight = weight_seq, 
                        stock = unique(gen_dat$stock)) %>% 
  left_join(., gen_dat %>% select(stock, smolt) %>% distinct(), by = "stock") %>% 
  gen_pred(a_mod2$gam, ., random = TRUE)

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


# As above but with scaled data model 
smooth_pars <- map(a_mod2b$gam$smooth, function(x) x$label) %>% unlist
weight_seq_z <- seq(min(gen_dat$avg_weight_z), max(gen_dat$avg_weight_z), 
                  length.out = 50)
fe_predsB <- expand.grid(avg_weight_z = weight_seq_z,
                        smolt = unique(gen_dat$smolt)) %>% 
  gen_pred(a_mod2b$gam, ., random = FALSE)
re_predsB <- expand.grid(avg_weight_z = weight_seq_z, 
                        stock = unique(gen_dat$stock)) %>% 
  left_join(., gen_dat %>% select(stock, smolt) %>% distinct(), by = "stock") %>% 
  gen_pred(a_mod2b$gam, ., random = TRUE)

fe_predsB %>% 
  ggplot(., aes(x = avg_weight_z)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  ggsidekick::theme_sleek() 

re_predsB %>% 
  ggplot(., aes(x = avg_weight_z)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  coord_cartesian(ylim = c(3, 6)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)


# MORTALITY RATE MODELS --------------------------------------------------------

# quick test for best model structure
s_mod1 <- gamm(survival ~ s(avg_weight, m = 2, bs = "tp", k = 3),
               random = list(stock = ~ 1),
               data = surv_dat, method = "REML", family=betar(link="logit"),
               correlation = corAR1(form = ~ brood_year | stock))
s_mod2 <- gamm(survival ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                 smolt, 
               random = list(stock = ~ 1),
               data = surv_dat, method = "REML", family=betar(link="logit"),
               correlation = corAR1(form = ~ brood_year | stock))
s_mod3 <- gamm(survival ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                 smolt,
               random = list(stock = ~ 1),
               data = surv_dat, method = "REML", family=betar(link="logit"),
               correlation = corAR1(form = ~ brood_year | stock))
s_mod4 <- gamm(survival ~ s(avg_weight, by = smolt, m = 2, bs = "tp", k = 3) +
                 smolt,
               random = list(stock = ~ 1 + avg_weight),
               data = surv_dat, method = "REML", family=betar(link="logit"),
               correlation = corAR1(form = ~ brood_year | stock))
AIC(s_mod1$lme, s_mod2$lme, s_mod3$lme, s_mod4$lme)

acf(resid(s_mod3$lme, type = "normalized"))
hist(resid(s_mod3$lme, type = "normalized"))


s_mod3b <- gamm(survival ~ s(avg_weight_z, by = smolt, m = 2, bs = "tp", k = 3) +
                  smolt,
                random = list(stock = ~ 1),
                data = surv_dat, method = "REML", family=betar(link="logit"),
                correlation = corAR1(form = ~ brood_year | stock))


## Check fit
fit_preds <- predict(s_mod3$gam, newdata = surv_dat, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
fit_dat <- fit_dat_f(surv_dat, preds = fit_preds, resp_dist = "beta") 

ggplot(data = fit_dat, aes(x = qlogis(survival), y = link_fit)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(aes(color = brood_year), size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = link_lo, 
                     ymax = link_up),
                 size = 1/2, color = "firebrick4") +
  labs(x = "Observed Logit Survival", 
       y = "Predicted Logit Survival") +
  theme_bw()


## Generate predictions
smooth_pars <- map(s_mod3$gam$smooth, function(x) x$label) %>% unlist
weight_seq <- seq(min(surv_dat$avg_weight), max(surv_dat$avg_weight), 
                  length.out = 50)
fe_preds <- expand.grid(avg_weight = weight_seq,
                        smolt = unique(surv_dat$smolt)) %>% 
  gen_pred(s_mod3$gam, ., random = FALSE)
re_preds <- expand.grid(avg_weight = weight_seq, 
                        stock = unique(surv_dat$stock)) %>% 
  left_join(., surv_dat %>% select(stock, smolt) %>% distinct(), by = "stock") %>% 
  gen_pred(s_mod3$gam, ., random = TRUE)

fe_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  geom_point(data = surv_dat, aes(x = avg_weight, y = survival, fill = smolt),
             shape = 21, alpha = 0.2) +
  ggsidekick::theme_sleek() 

re_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)


# As above but with scaled data model 
smooth_pars <- map(s_mod3b$gam$smooth, function(x) x$label) %>% unlist
weight_seq_z <- seq(min(surv_dat$avg_weight_z), max(surv_dat$avg_weight_z), 
                    length.out = 50)
fe_predsB <- expand.grid(avg_weight_z = weight_seq_z,
                         smolt = unique(surv_dat$smolt)) %>% 
  gen_pred(s_mod3b$gam, ., random = FALSE)
re_predsB <- expand.grid(avg_weight_z = weight_seq_z, 
                         stock = unique(surv_dat$stock)) %>% 
  left_join(., surv_dat %>% select(stock, smolt) %>% distinct(), by = "stock") %>% 
  gen_pred(s_mod3b$gam, ., random = TRUE)

fe_predsB %>% 
  ggplot(., aes(x = avg_weight_z)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  geom_point(data = surv_dat, aes(x = avg_weight_z, y = survival, fill = smolt),
             shape = 21, alpha = 0.2) +
  ggsidekick::theme_sleek() 

re_predsB %>% 
  ggplot(., aes(x = avg_weight_z)) +
  geom_line(aes(y = pred, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_up, fill = smolt), 
              alpha = 0.3) +
  geom_point(data = surv_dat, aes(x = avg_weight_z, y = survival, fill = smolt),
             shape = 21, alpha = 0.2) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)

