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


# FIT GAMS ---------------------------------------------------------------------

# Test various model structures including single, additive or interaction models
# between herring (various lags) and RKW with the potential for stratification
# by adult rearing location

h_mod1 <- gam(gen_length ~ s(herr_OEY_z, m = 2, bs = "tp", k = 4) + 
                s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 4) + 
                s(stock, bs = "re"), 
              data = dat)
h_mod2a <- gam(gen_length ~  s(herr_OEY_z, m = 2, bs = "tp", k = 4) +
                 s(herr_OEY_z, by = a_group2, m = 1, bs = "tp", k = 4) + 
                 a_group2 +
                 s(stock, bs = "re"), 
               data = dat)
h_mod2b <- gam(gen_length ~ s(herr_OEY_z, by = a_group2, m = 2, bs = "tp", k = 4) + 
                 a_group2 +
                 s(stock, bs = "re"), 
               data = dat)
h_mod3a <- gam(gen_length ~ s(herr_OEY_z, by = a_group2, m = 2, bs = "tp", k = 4) + 
                s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 4) + 
                a_group2 +
                s(stock, bs = "re"), 
              data = dat)
h_mod3b <- gam(gen_length ~ s(herr_OEY_z, by = a_group2, m = 2, bs = "tp", k = 4) + 
               s(herr_OEY_z, by = stock, m = 1, bs = "tp", k = 4) + 
                s(a_group2, bs = "re") +
                s(stock, bs = "re"), 
             data = dat)
AIC(h_mod1, h_mod2a, h_mod2b, h_mod3)
AIC(h_mod2a, h_mod2b, h_mod2c)



# generate predictions for interaction model
excl_pars <- map(h_mod1$smooth, function(x) x$label) %>% unlist
herr_seq <- seq(min(dat$herr_OEY_z), max(dat$herr_OEY_z), length.out = 75)
new_dat <- expand.grid(#seal_anom = seal_seq, 
                       herr_OEY_z = herr_seq,
                       a_group2 = unique(dat$a_group2))

# fixed effects predictions
preds <- predict(h_mod2c, new_dat, se.fit = TRUE, 
                 exclude = excl_pars[grepl("stock", excl_pars)],
                 newdata.guaranteed = TRUE)
new_dat2 <- new_dat %>% 
  mutate(link_fit = as.numeric(preds$fit),
         link_se = as.numeric(preds$se.fit),
         pred_gen = plogis(link_fit),
         pred_gen_lo = plogis(link_fit + (qnorm(0.025) * link_se)),
         pred_gen_up = plogis(link_fit + (qnorm(0.975) * link_se))
  )

# herring only fixed effects
# zero_seal <- new_dat2$seal_anom[which.min(abs(new_dat2$seal_anom - 0))]
new_dat2 %>% 
  # filter(seal_anom == zero_seal) %>% 
  ggplot(., aes(x = herr_OEY_z)) +
  geom_line(aes(y = pred_gen)) +
  geom_ribbon(aes(ymin = pred_gen_lo, ymax = pred_gen_up), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  # lims(y = c(0.01, 0.05)) +
  labs(x = "Herring Anomaly", y = "Predicted Generation Length") +
  facet_wrap(~a_group2)



heat_plot <- ggplot(new_dat2, aes(x = seal_anom, y = herr_anom)) +
  geom_raster(aes(fill = pred_surv)) +
  scale_fill_viridis_c(name = "Predicted\nSurvival") +
  labs(x = "SoG Seal Anomaly", y = "SoG Herring Anomaly") +
  ggsidekick::theme_sleek()

