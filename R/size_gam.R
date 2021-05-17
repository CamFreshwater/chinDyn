## Generation length gams
# Hierarchical analysis of size at release on juvenile survival and mean age
# May 17, 2021

library(mgcv)
library(tidyverse)

chin_dat <- readRDS(here::here("data", "salmon_data", 
                              "cwt_indicator_surv_clean.RDS")) %>% 
  filter(!is.na(avg_weight)) %>% 
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


# AGE MODELS -------------------------------------------------------------------

# quick test for best model structure
a_mod1 <- gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = chin_dat, method = "REML", family=Gamma(link="log"))
a_mod2 <- gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = chin_dat, method = "REML", family=Gamma(link="log"))
a_mod3 <- gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = smolt, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = chin_dat, method = "REML", family=Gamma(link="log"))
a_mod4 <- gam(gen_length ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = chin_dat, method = "REML", family=Gamma(link="log"))
AIC(a_mod1, a_mod2, a_mod3, a_mod4)


## Check fit
chin_dat_gen <- chin_dat %>% filter(!is.na(gen_length))
fit_preds <- predict(a_mod4, newdata = chin_dat_gen, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
fit_dat <- chin_dat %>% 
  filter(!is.na(gen_length)) %>% 
  mutate(link_fit = as.numeric(fit_preds$fit),
         link_se = as.numeric(fit_preds$se.fit),
         link_lo = link_fit + (qnorm(0.025) * link_se),
         link_up = link_fit + (qnorm(0.975) * link_se),
         pred_gen = exp(link_fit),
         pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
         pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
  )

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
acf(resid(a_mod1))


## Generate predictions

# function to generate predictions
excl_pars <- map(mod$smooth, function(x) x$label) %>% unlist
gen_pred <- function(mod_in, dat_in, random = FALSE) {
  if (random == FALSE) {
    preds <- predict(mod_in, dat_in, se.fit = TRUE, 
                     exclude = excl_pars[grepl("stock", excl_pars)],
                     newdata.guaranteed = TRUE)
  } else {
    preds <- predict(mod_in, dat_in, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
  }
  
  #remove zero columns 
  dat_in %>% 
    mutate(link_fit = as.numeric(preds$fit),
           link_se = as.numeric(preds$se.fit),
           pred_gen = exp(link_fit),
           pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
           pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
    ) 
}


weight_seq <- seq(min(chin_dat$avg_weight), max(chin_dat$avg_weight), 
                  length.out = 50)
fe_preds <- data.frame(avg_weight = weight_seq) %>% 
  gen_pred(a_mod1, ., random = FALSE)
re_preds <- expand.grid(avg_weight = weight_seq, 
                        stock = unique(chin_dat_gen$stock)) %>% 
  gen_pred(a_mod1, ., random = TRUE) %>% 
  left_join(., chin_dat %>% select(stock, smolt), by = "stock")


fe_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up), 
              alpha = 0.3) +
  ggsidekick::theme_sleek() 

re_preds %>% 
  ggplot(., aes(x = avg_weight)) +
  geom_line(aes(y = pred_gen, colour = smolt)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up, fill = smolt), 
              alpha = 0.3) +
  lims(y = c(3, 6)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)


a_mod1 <- gam(M ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = chin_dat, method = "REML")
a_mod2 <- gam(M ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = chin_dat, method = "REML")
a_mod3 <- gam(M ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = smolt, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = chin_dat, method = "REML")
a_mod4 <- gam(M ~ s(avg_weight, m = 2, bs = "tp", k = 3) +
                s(avg_weight, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re") +
                smolt,
              data = chin_dat, method = "REML")
AIC(a_mod1, a_mod2, a_mod3)