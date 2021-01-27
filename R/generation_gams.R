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

# SOG data 
sog_ocean <- adult_cov$sog %>%
  mutate(ck_ocean_year0 = year) %>% 
  select(-year)


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
  # ocean entry year oceanography during priming period (April-July)
  left_join(.,
            sog_ocean %>% 
              rename(year = ck_ocean_year0),
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
    rkw_z = as.numeric(scale(rollmean_rkw)),
    sst_z = as.numeric(scale(prime_entrance_sst)),
    sal_z = as.numeric(scale(prime_entrance_salinity))
  ) %>% 
  select(stock, stock_name, brood_year, a_group2, gen_length, gen_z, 
         herr_abund_OEY, herr_abund_OEY1, rollmean_rkw, herr_OEY_z, herr_OEY1_z,
         rkw_z, sst_z, sal_z) %>% 
  # remove missing data
  filter(!is.na(rkw_z)#,
         # !is.na(sst_z)
         )


# look at correlations among covariates (minimal)
dat %>%
  select(brood_year, herr_OEY_z, herr_OEY1_z, rkw_z, sst_z, sal_z) %>% 
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


dat %>% 
  select(a_group2, gen_length, herr_OEY_z:rkw_z) %>% 
  pivot_longer(., cols = -a_group2, names_to = "var", values_to = "val") %>% 
  ggplot(., aes(x = a_group2, y = val)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~var)

## NOTE: RKW abundance highly correlated with adult rearing location probably
# can't include categorical variable as covariate


# GAM COMPARISON ---------------------------------------------------------------

# Test various model structures including single, additive or interaction models
# between herring (various lags) and RKW with the potential for stratification
# by adult rearing location

# Herring at one year lags: model structure supports model 2 based on parsimony
h_mod1 <- gam(gen_length ~ s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = dat, method = "REML", family=Gamma(link="log"))
h_mod2 <- gam(gen_length ~ s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                s(stock, bs = "re"),
              data = dat, method = "REML", family=Gamma(link="log"))
# h_mod2 <- gam(gen_length ~ s(herr_OEY1_z, by = a_group2, m = 2, bs = "tp", k = 3) +
#                 a_group2 +
#                 s(stock, bs = "re"),
#               data = dat, method = "REML", family=Gamma(link="log"))
# h_mod3 <- gam(gen_length ~ s(herr_OEY1_z, by = a_group2, m = 2, bs = "tp", k = 3) +
#                 s(herr_OEY1_z, by = stock, m = 1, bs = "tp", k = 3) +
#                 s(stock, bs = "re"),
#               data = dat, method = "REML", family=Gamma(link="log"))
# h_mod4 <- gam(gen_length ~ s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
#                 a_group2 +
#                 s(stock, bs = "re"),
#               data = dat, method = "REML", family=Gamma(link="log"))
AIC(h_mod1, h_mod2)


# Herring at entry - global smooth for juveniles strongly supported
# No region effect because all common
hj_mod1 <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 4) +
                  s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 4) +
                  s(stock, bs = "re"),
               data = dat, method = "REML", family=Gamma(link="log"))
hj_mod2 <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 4) +
                  s(stock, bs = "re"),
               data = dat, method = "REML", family=Gamma(link="log"))
AIC(hj_mod1, hj_mod2)


# SST at entry - global smooth for juveniles strongly supported
# No region effect because all common
sst_mod1 <- gam(gen_length ~ s(sst_z, m = 2, bs = "tp", k = 4) +
                 s(sst_z, by = stock, m = 1, bs = "tp", k = 4) +
                 s(stock, bs = "re"),
               data = dat, method = "REML", family=Gamma(link="log"))
sst_mod2 <- gam(gen_length ~ s(sst_z, m = 2, bs = "tp", k = 4) +
                 s(stock, bs = "re"),
               data = dat, method = "REML", family=Gamma(link="log"))
AIC(sst_mod1, sst_mod2)


# Salinity at entry - global smooth for juveniles strongly supported
# No region effect because all common
sal_mod1 <- gam(gen_length ~ s(sal_z, m = 2, bs = "tp", k = 4) +
                 s(sal_z, by = stock, m = 1, bs = "tp", k = 4) +
                 s(stock, bs = "re"),
               data = dat, method = "REML", family=Gamma(link="log"))
sal_mod2 <- gam(gen_length ~ s(sal_z, m = 2, bs = "tp", k = 4) +
                 s(stock, bs = "re"),
               data = dat, method = "REML", family=Gamma(link="log"))
AIC(sal_mod1, sal_mod2)


# RKW - greatest support for model 1 (overall smooth plus random stock 
# effects)
# k_mod1 <- gam(gen_length ~  s(rkw_z, by = a_group2, m = 2, bs = "tp", k = 4) +
#                s(rkw_z, by = stock, m = 1, bs = "tp", k = 4) +
#                s(stock, bs = "re") +
#                a_group2,
#               data = dat, method = "REML", family=Gamma(link="log"))
k_mod1 <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 4) +
                s(rkw_z, by = stock, m = 1, bs = "tp", k = 4) +
                 s(stock, bs = "re"),
              data = dat, method = "REML", family=Gamma(link="log"))
# k_mod3 <- gam(gen_length ~ s(rkw_z, by = a_group2, m = 2, bs = "tp", k = 4) +
#                 s(stock, bs = "re"),
#               data = dat, method = "REML", family=Gamma(link="log"))
k_mod2 <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 4) +
                s(stock, bs = "re"),
              data = dat, method = "REML", family=Gamma(link="log"))
# k_mod5 <- gam(gen_length ~  s(rkw_z, by = a_group2, m = 2, bs = "tp", k = 4) +
#                 s(rkw_z, by = stock, m = 1, bs = "tp", k = 4) +
#                 s(stock, bs = "re"),
              # data = dat, method = "REML", family=Gamma(link="log"))
AIC(k_mod1, k_mod2)


# Combined models to resolve whether relationships are additive
n_stks <- length(unique(dat$stock))
cat_mod <- gam(gen_length ~ a_group2 ,
                data = dat, method = "REML", family=Gamma(link="log"))
h_hj_mod <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
                  s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
                  s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                  s(stock, bs = "re") +
                  a_group2 ,
                data = dat, method = "REML", family=Gamma(link="log"))
k_hj_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                   s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                   s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
                   s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
                   s(stock, bs = "re"),
              data = dat, method = "REML", family=Gamma(link="log"))
k_hj_int_mod <- gam(gen_length ~ te(herr_OEY_z, rkw_z, 
                                   bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                     t2(herr_OEY_z, rkw_z, stock, m = 2,
                        bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)),
              data = dat, method = "REML", family=Gamma(link="log"))
k_sst_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                 s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                 s(sst_z, m = 2, bs = "tp", k = 3) +
                 s(stock, bs = "re") ,
               data = dat, method = "REML", family=Gamma(link="log"))
k_sal_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                 s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                 s(sal_z, m = 2, bs = "tp", k = 3) +
                 s(stock, bs = "re") ,
               data = dat, method = "REML", family=Gamma(link="log"))
k_h_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
                  s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
                  s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
                  s(stock, bs = "re") ,
              data = dat, method = "REML", family=Gamma(link="log"))
k_h_int_mod <- gam(gen_length ~ te(herr_OEY1_z, rkw_z, 
                                  bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                     t2(herr_OEY1_z, rkw_z, stock, m = 2,
                       bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)),
              data = dat, method = "REML", family=Gamma(link="log"))
k_h_hj_mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
             s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(stock, bs = "re"),
           data = dat,
           method = "REML",
           family=Gamma(link="log"))
k_h_hj_int_mod <- gam(gen_length ~ te(herr_OEY1_z, rkw_z, 
                                      bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                        t2(herr_OEY1_z, rkw_z, stock, m = 2,
                           bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)) +
                        te(herr_OEY_z, rkw_z, 
                           bs = c("tp", "tp"), k = c(3, 3), m = 2) +
                        t2(herr_OEY_z, rkw_z, stock, m = 2,
                           bs = c("tp", "tp", "re"), k = c(3, 3, n_stks)),
                      data = dat,
                      method = "REML", family=Gamma(link="log"))
AIC(cat_mod, h_mod2, hj_mod1, k_mod1, h_hj_mod, k_hj_mod, k_h_mod, k_h_int_mod, 
    k_h_hj_mod, k_hj_int_mod, k_h_hj_int_mod, k_sst_mod, k_sal_mod)
# top supported model is most complex additive (i.e. no support for 
# interactions)


# GAM PREDICTIONS --------------------------------------------------------------

# equivalent model to k_h_hj_mod
mod <- gam(gen_length ~ s(rkw_z, m = 2, bs = "tp", k = 3) +
             s(rkw_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(herr_OEY1_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY_z, m = 2, bs = "tp", k = 3) +
             s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 3) +
             s(stock, bs = "re"),
           data = dat,
           method = "REML",
           family=Gamma(link="log"))


## Check fit (note that fit is on scale of observations so unlikely to violate
# assumption in link space)
fit_preds <- predict(mod, newdata = dat, se.fit = TRUE, 
                     newdata.guaranteed = TRUE)
fit_dat <- dat %>% 
  mutate(link_fit = as.numeric(fit_preds$fit),
         link_se = as.numeric(fit_preds$se.fit),
         pred_gen = exp(link_fit),
         pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
         pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
  )

ggplot(data = fit_dat, aes(x = gen_length, y = pred_gen)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(aes(color = brood_year), size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = pred_gen_lo, 
                     ymax = pred_gen_up),
                 size = 1/2, color = "firebrick4") +
  labs(x = "Observed Generation Length", 
       y = "Predicted Generation Length") +
  theme_bw()

# minimal evidence of temporal autocorrelation
acf(resid(mod))


## Generate predictions

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
    mutate(link_fit = as.numeric(preds$fit),
           link_se = as.numeric(preds$se.fit),
           pred_gen = exp(link_fit),
           pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
           pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
    ) 
  colnames(dat_out)[1] <- "cov1"
  
  return(dat_out)
}


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
  left_join(., dat %>% select(stock, a_group2), by = "stock") %>% 
  mutate(a_group2 = fct_relevel(a_group2, "offshore", "north", "broad", 
                                "south"))
fixed_pred_dat <- map(fixed_pred_list1, function (x) {
  gen_pred(mod, x, random = FALSE) 
}) %>% 
  bind_rows()

# plot fixed effects
fe_splines <- ggplot(fixed_pred_dat, aes(x = cov1)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up), 
              alpha = 0.3) +
  ggsidekick::theme_sleek() +
  labs(x = "Standardized Abundance", y = "Predicted Generation Length") +
  facet_wrap(~var)
  

# heat plots
new_heat_dat <- expand.grid(herr_OEY_z = herr0_seq,
                            herr_OEY1_z = herr1_seq,
                            rkw_z = rkw_seq)
heat_preds <- predict(mod, new_heat_dat, se.fit = TRUE, 
                   exclude = excl_pars[grepl("stock", excl_pars)],
                   newdata.guaranteed = TRUE)
heat_dat <- new_heat_dat %>% 
  mutate(link_fit = as.numeric(heat_preds$fit),
         link_se = as.numeric(heat_preds$se.fit),
         pred_gen = exp(link_fit),
         pred_gen_lo = exp(link_fit + (qnorm(0.025) * link_se)),
         pred_gen_up = exp(link_fit + (qnorm(0.975) * link_se))
  )

heat_plot <- heat_dat %>%
  filter(herr_OEY1_z == median(herr1_seq)) %>% 
  distinct() %>% 
  ggplot(., aes(x = rkw_z, y = herr_OEY_z)) +
  geom_raster(aes(fill = pred_gen)) +
  scale_fill_viridis_c(name = "Predicted\nGeneration\nLength") +
  labs(x = "RKW Anomaly", y = "SoG Herring Anomaly") +
  ggsidekick::theme_sleek()



# plot mixed effects w/ observations
obs_dat <- dat %>% 
  select(stock, a_group2, gen_length, herr_OEY_z, herr_OEY1_z, rkw_z) %>% 
  pivot_longer(., cols = c(herr_OEY_z, herr_OEY1_z, rkw_z), 
               names_to = "var", values_to = "z_value") %>% 
  mutate(var = fct_recode(as.factor(var), 
                          herr0 = "herr_OEY_z",
                          herr1 = "herr_OEY1_z",
                          rkw = "rkw_z"))

mixed_splines <- ggplot(mixed_pred_dat, aes(x = cov1)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up, fill = a_group2), 
              alpha = 0.3) +
  geom_point(data = obs_dat, aes(x = z_value, y = gen_length, fill = a_group2),
             shape = 21, alpha = 0.5) +
  ggsidekick::theme_sleek() +
  labs(x = "Standardized Abundance", y = "Generation Length") +
  facet_grid(var~fct_reorder(stock, as.numeric(a_group2))) +
  coord_cartesian(ylim = c(2, 6)) +
  theme(legend.position = "top")


png(here::here("figs", "gam", "gen_gam", "herr_rkw_fe_heatplot.png"), 
    height = 5, width = 6, res = 300, units = "in")
heat_plot
dev.off()

png(here::here("figs", "gam", "gen_gam", "herr_rkw_fe_splines.png"), 
    height = 2.5, width = 6, res = 300, units = "in")
fe_splines
dev.off()

png(here::here("figs", "gam", "gen_gam", "herr_rkw_mixed_splines.png"), 
    height = 3.5, width = 9, res = 300, units = "in")
mixed_splines 
dev.off()


