## Generation length gams
# Hierarchical analysis of RKW and herring abundance on mean generation length
# For now focused only on Strait of Georgia stocks to match overlap with 
# available pinniped data
# Jan. 14, 2021


library(mgcv)
library(tidyverse)

gen_raw <- readRDS(here::here("data", "salmonData", 
                              "cwt_indicator_surv_clean.RDS")) 

#remove stocks with no gen_length data
gen <- gen_raw %>% 
  filter(!is.na(gen_length)) %>% 
  mutate(gen_z = as.numeric(scale(gen_length)))

# dataframe of only stocks and adult groupings
stk_tbl <- gen %>% 
  group_by(stock) %>% 
  mutate(max_age = ceiling(max(gen_length)),
         max_ocean_age = ifelse(smolt == "oceantype", max_age - 1, 
                                max_age - 2)) %>% 
  ungroup() %>% 
  select(stock, stock_name, smolt, run, j_group2, a_group:a_group3,
         max_age, max_ocean_age) %>% 
  distinct() 

# prepare covariate data
adult_cov <- readRDS(here::here("data/salmonData/cov_subset_adult.rds"))


# PREP DATA --------------------------------------------------------------------

# Focus on Strait of Georgia stocks and add herring at two lags to test timing
# of prey availability effects: 
# 1) year of ocean entry (Strait of Georgia herring only)
# 2) year after ocean entry (WCVI + SoG for most SoG stocks, PRD + HG + CC for
# north migrating)
herr_year0 <- adult_cov$herring %>% 
  filter(stock == "SoGHerringR") %>% 
  mutate(
    ck_ocean_year0 = herring_age0_year) %>% 
  select(herr_age0_abund = herr_abund, ck_ocean_year0)

herr_year1 <- adult_cov$herring %>% 
  mutate(
    herr_reg = case_when(
      stock %in% c("WCVIHerringR", "SoGHerringR") ~ "south",
      stock %in% c("HGHerringR", "CCHerringR", "PRDHerringR") ~ "north"
    ),
    ck_ocean_year1 = herring_age0_year - 1) %>% 
  group_by(ck_ocean_year1, herr_reg) %>% 
  summarize(agg_abund = sum(herr_abund), .groups = "drop") %>% 
  pivot_wider(names_from = herr_reg, values_from = agg_abund,
              names_prefix = "herr_abund_OEY1_") 


# calculate rkw abundance for northern stocks (NRKW) and southern stocks (SRKW +
# NRKW) by brood year (rolling average of brood year + 1:4)
sog_stks <- stk_tbl %>% 
  filter(j_group2 == "sog") 

# add RKW data to each stock
rkw <- expand.grid(stock = sog_stks$stock,
                   year = unique(adult_cov$rkw$year)) %>% 
  arrange(stock) %>%
  left_join(., stk_tbl, by = "stock") %>% 
  select(stock, stock_name, smolt, a_group2, max_age, max_ocean_age, year) %>% 
  droplevels() %>% 
  left_join(., adult_cov$rkw, by = "year") %>% 
  mutate(
    #identify which rkw pop is relevant
    relevant_rkw = case_when(
      a_group2 == "south" ~ total_n,
      a_group2 == "offshore" ~ srkw_n,
      a_group2 == "north" ~ nrkw_n
    ),
    #calculate brood year based on smolt type to pair with generation data
    brood_year = year - max_age
  ) %>% 
  group_by(stock) %>% 
  mutate(
    #calculate rolling mean based on the maximum ocean age of each stock
    rollmean_rkw = zoo::rollmean(relevant_rkw, max_ocean_age, fill = NA, 
                                 align = "right")
  ) %>% 
  ungroup() 


# join relevant stock data, herring data and rkw
dat <- gen %>% 
  filter(j_group2 == "sog") %>% 
  # ocean entry year herring abundance
  left_join(.,
            herr_year0 %>% 
              select(year = ck_ocean_year0, herr_abund_OEY = herr_age0_abund),
            by = "year") %>% 
  # ocean entry year + 1 herring abundance
  left_join(., 
            herr_year1 %>% 
              rename(year = ck_ocean_year1),
            by = "year") %>%
  # define oey1 herring abundance based on adult distribution
  mutate(
    herr_abund_OEY1 = case_when(
      a_group2 == "south" ~ herr_abund_OEY1_south,
      a_group2 == "offshore" ~ herr_abund_OEY1_north,
      a_group2 == "north" ~ herr_abund_OEY1_north
    )
  ) %>% 
  # rkw abundance
  left_join(., 
            rkw %>% 
              select(stock, brood_year, rollmean_rkw), 
            by = c("brood_year", "stock")) %>% 
  # scale covariates
  mutate(
    herr_OEY_z = as.numeric(scale(herr_abund_OEY)),
    herr_OEY1_z = as.numeric(scale(herr_abund_OEY1)),
    rkw_z = as.numeric(scale(rollmean_rkw))
  ) %>% 
  select(stock, stock_name, brood_year, a_group2, gen_length, gen_z, 
         herr_abund_OEY, herr_abund_OEY1, rollmean_rkw, herr_OEY_z, herr_OEY1_z,
         rkw_z) %>% 
  # remove missing data
  filter(!is.na(rkw_z))


# look at correlations among covariates (minimal)
dat %>%
  select(brood_year, herr_OEY_z, herr_OEY1_z, rkw_z) %>% 
  as.matrix() %>% 
  cor() %>% 
  corrplot::corrplot(., method = "number", "upper")


# look at raw data
p <- ggplot(dat) +
  facet_wrap(~fct_reorder(stock, as.numeric(a_group2))) +
  ggsidekick::theme_sleek()
p +
  geom_point(aes(x = herr_OEY_z, y = gen_length, fill = a_group2), shape = 21)
p +
  geom_point(aes(x = herr_OEY1_z, y = gen_length, fill = a_group2), shape = 21)
p +
  geom_point(aes(x = rkw_z, y = gen_length, fill = a_group2), shape = 21)


# GAM COMPARISON ---------------------------------------------------------------

# Test various model structures including single, additive or interaction models
# between herring (various lags) and RKW with the potential for stratification
# by adult rearing location

# Herring at one year lags: model structure supports model 1
h_mod1 <- gam(gen_length ~ s(herr_OEY1_z, m = 2, bs = "tp", k = 4) +
                s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 4) +
                s(stock, bs = "re"),
              data = dat)
h_mod2 <- gam(gen_length ~ s(herr_OEY1_z, by = a_group2, m = 2, bs = "tp", k = 4) +
                a_group2 +
                s(stock, bs = "re"),
              data = dat)
h_mod3 <- gam(gen_length ~ s(herr_OEY1_z, by = a_group2, m = 2, bs = "tp", k = 4) +
                s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 4) +
                a_group2 +
                s(stock, bs = "re"),
              data = dat)
h_mod4 <- gam(gen_length ~ s(herr_OEY1_z, m = 2, bs = "tp", k = 4) +
                a_group2 +
                s(stock, bs = "re"),
              data = dat)
AIC(h_mod1, h_mod2, h_mod3, h_mod4)


# Herring at entry - global smooth for juveniles strongly supported
hj_mod1 <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 4) +
                  s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 4) +
                  s(stock, bs = "re"),
                data = dat)
hj_mod2 <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 4) +
                  s(stock, bs = "re"),
                data = dat)
AIC(hj_mod1, hj_mod2)


# RKW - greatest support for model 2 (global smooth and stock effects)
k_mod1 <- gam(gen_length ~  s(rkw_z, by = a_group2, m = 2, bs = "tp", k = 4) +
               s(rkw_z, by = stock, m = 1, bs = "tp", k = 4) +
               s(stock, bs = "re") +
               a_group2,
               data = dat)
k_mod2 <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 4) +
                s(rkw_z, by = stock, m = 1, bs = "tp", k = 4) +
                 s(stock, bs = "re"),
               data = dat)
k_mod3 <- gam(gen_length ~ s(rkw_z, by = a_group2, m = 2, bs = "tp", k = 4) +
                s(stock, bs = "re"),
              data = dat)
k_mod4 <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 4) +
                s(stock, bs = "re"),
              data = dat)
AIC(k_mod1, k_mod2, k_mod3, k_mod4)


# Combined models to resolve whether relationships are additive
n_stks <- length(unique(dat$stock))
h_hj_mod <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
                  s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
                  s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                  s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
                  s(stock, bs = "re"),
              data = dat,
              method = "REML")
k_hj_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                   s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                   s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
                   s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
                   s(stock, bs = "re"),
                data = dat,
                method = "REML")
k_hj_int_mod <- gam(gen_length ~ te(herr_OEY_z, rkw_z, 
                                   bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                     t2(herr_OEY_z, rkw_z, stock, m = 2,
                        bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)),
                   data = dat,
                   method = "REML")
k_h_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                  s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                  s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                  s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
                  s(stock, bs = "re"),
              data = dat,
              method = "REML")
k_h_int_mod <- gam(gen_length ~ te(herr_OEY1_z, rkw_z, 
                                  bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                     t2(herr_OEY1_z, rkw_z, stock, m = 2,
                       bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)),
               data = dat,
               method = "REML")
k_h_hj_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                    s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                    s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                    s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
                    s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
                    s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
                    s(stock, bs = "re"),
               data = dat,
               method = "REML")
k_h_hj_int_mod <- gam(gen_length ~ te(herr_OEY1_z, rkw_z, 
                                      bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                        t2(herr_OEY1_z, rkw_z, stock, m = 2,
                           bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)) +
                        te(herr_OEY_z, rkw_z, 
                           bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                        t2(herr_OEY_z, rkw_z, stock, m = 2,
                           bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)),
                      data = dat,
                      method = "REML")
AIC(h_mod1, hj_mod1, k_mod2, h_hj_mod, k_hj_mod, k_h_mod, k_h_int_mod, 
    k_h_hj_mod, k_hj_int_mod, k_h_hj_int_mod)
AIC(k_h_hj_mod, k_h_hj_modB)
# top supported model is most complex additive (i.e. no support for 
# interactions)


# GAM PREDICTIONS --------------------------------------------------------------

# function to generate predictions
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
  dat_out <- dat_in[c(-2, -3)] %>% 
    mutate(pred_gen = as.numeric(preds$fit),
           pred_gen_se = as.numeric(preds$se.fit),
           pred_gen_lo = pred_gen + (qnorm(0.025) * pred_gen_se),
           pred_gen_up = pred_gen + (qnorm(0.975) * pred_gen_se)
    ) 
  colnames(dat_out)[1] <- "cov1"
  
  return(dat_out)
}

# equivalent model to k_h_hj_mod
mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
             s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(stock, bs = "re"),
           data = dat,
           method = "REML")
mod2 <- gam(gen_length ~ s(rkw_z, by = a_group2, m = 2, bs = "tp", k = 3) +
             s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(stock, bs = "re"),
           data = dat,
           method = "REML")

# input vectors
excl_pars <- map(mod$smooth, function(x) x$label) %>% unlist
herr0_seq <- seq(min(dat$herr_OEY_z), max(dat$herr_OEY_z), length.out = 75)
herr1_seq <- seq(min(dat$herr_OEY1_z), max(dat$herr_OEY1_z), length.out = 75)
rkw_seq <- seq(min(dat$rkw_z), max(dat$rkw_z), length.out = 75)

# stock-specific averages for each covariate
stock_means <- dat %>% 
  group_by(stock) %>% 
  summarize(herr_OEY_z = mean(herr_OEY_z),
            herr_OEY1_z = mean(herr_OEY1_z),
            rkw_z = mean(rkw_z),
            .groups = "drop") %>% 
  ungroup()

# make mixed effects predictions using stock mean values
# note that order of covariates in DF varies and is important
mixed_preds_herr0 <- expand.grid(herr_OEY_z = herr0_seq,
                                 stock = unique(dat$stock),
                                 var = "herr0") %>% 
  left_join(., stock_means %>% select(-herr_OEY_z), by = "stock") %>% 
  select(herr_OEY_z, herr_OEY1_z, rkw_z, stock, var)
mixed_preds_herr1 <- expand.grid(herr_OEY1_z = herr1_seq,
                                 stock = unique(dat$stock),
                                 var = "herr1") %>% 
  left_join(., stock_means %>% select(-herr_OEY1_z), by = "stock") %>% 
  select(herr_OEY1_z, herr_OEY_z, rkw_z, stock, var)
mixed_preds_rkw <- expand.grid(rkw_z = rkw_seq,
                               stock = unique(dat$stock),
                               var = "rkw") %>% 
  left_join(., stock_means %>% select(-rkw_z), by = "stock") %>% 
  select(rkw_z, herr_OEY_z, herr_OEY1_z, stock, var)
mixed_pred_list1 <- list(herr0 = mixed_preds_herr0,
                         herr1 = mixed_preds_herr1,
                         rkw = mixed_preds_rkw)
fixed_pred_list1 <- map(mixed_pred_list1, function(x) {
  #replace stock-specific covariates with 0s 
  x[, 2] <- 0
  x[, 3] <- 0 
  x %>% select(-stock) %>% distinct()
})

mixed_pred_dat <- map(mixed_pred_list1, function (x) {
  gen_pred(mod, x, random = TRUE) 
}) %>% 
  bind_rows() %>%
  left_join(., dat %>% select(stock, a_group2) %>% distinct(), by = "stock")
fixed_pred_dat <- map(fixed_pred_list1, function (x) {
  gen_pred(mod, x, random = FALSE) 
}) %>% 
  bind_rows()

# plot fixed effects
ggplot(fixed_pred_dat, aes(x = cov1)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  labs(x = "Standardized Abundance", y = "Predicted Generation Length") +
  facet_wrap(~var)
  
# plot mixed effects w/ observations
obs_dat <- dat %>% 
  select(stock, a_group2, gen_length, herr_OEY_z, herr_OEY1_z, rkw_z) %>% 
  pivot_longer(., cols = c(herr_OEY_z, herr_OEY1_z, rkw_z), 
               names_to = "var", values_to = "z_value") %>% 
  mutate(var = fct_recode(as.factor(var), 
                          herr0 = "herr_OEY_z",
                          herr1 = "herr_OEY1_z",
                          rkw = "rkw_z"))

ggplot(mixed_pred_dat, aes(x = cov1)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up, fill = a_group2), 
              alpha = 0.3) +
  geom_point(data = obs_dat, aes(x = z_value, y = gen_length, fill = a_group2),
             shape = 21, alpha = 0.5) +
  ggsidekick::theme_sleek() +
  labs(x = "Standardized Abundance", y = "Generation Length") +
  facet_grid(var~fct_reorder(stock, as.numeric(a_group2)))



gen_pred2 <- function(mod_in, dat_in, random = FALSE) {
  preds <- predict(mod_in, dat_in, se.fit = TRUE, 
                   newdata.guaranteed = TRUE)
  
  
  #remove zero columns 
  dat_out <- dat_in[c(-2, -3)] %>% 
    mutate(link_fit = as.numeric(preds$fit),
           link_se = as.numeric(preds$se.fit),
           pred_gen = exp(link_fit),
           pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
           pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
    ) 
  colnames(dat_out)[1] <- "cov1"
  
  return(dat_out)
}




# herring only fixed effects
# zero_seal <- new_dat2$seal_anom[which.min(abs(new_dat2$seal_anom - 0))]
new_dat2 %>% 
  ggplot(., aes(x = herr_OEY_z)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  labs(x = "Herring Anomaly", y = "Predicted Generation Length") +
  facet_wrap(~a_group2)



heat_plot <- ggplot(new_dat2, aes(x = seal_anom, y = herr_anom)) +
  geom_raster(aes(fill = pred_surv)) +
  scale_fill_viridis_c(name = "Predicted\nSurvival") +
  labs(x = "SoG Seal Anomaly", y = "SoG Herring Anomaly") +
  ggsidekick::theme_sleek()




