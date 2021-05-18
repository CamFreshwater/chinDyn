## Generation length bayesian
# Hierarchical analysis of size at release on juvenile survival and mean age
# May 17, 2021
## ABANDONED FOR NOW

library(tidyverse)
library(brms)
library(tidybayes)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 2)


chin_dat <- readRDS(here::here("data", "salmon_data", 
                               "cwt_indicator_surv_clean.RDS")) %>% 
  filter(!is.na(avg_weight)) 

gen_dat <- chin_dat %>%
  filter(!is.na(gen_length)) %>% 
  droplevels()
gen_dat <- chin_dat %>%
  filter(!is.na(gen_length)) %>% 
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
  ggsidekick::theme_sleek() +
  geom_point(aes(x = avg_weight, y = gen_length, fill = smolt), shape = 21)


# AGE MODELS -------------------------------------------------------------------

age1 <- brm(data = gen_dat, family = Gamma(link="log"),
            gen_length ~ 1 + avg_weight, 
            prior = c(prior(normal(0, 2), class = "Intercept"),
                      prior(normal(0, 2), class = "b"),
                      prior(gamma(0.01, 0.01), class = "shape")),
            iter = 1000, warmup = 500, chains = 4, cores = 4,
            seed = 13,
            control = list(adapt_delta = 0.9))
posterior_summary(age1)
acf(residuals(age1)[, "Estimate"])

# correlation structure dramatically increases runtime, doesn't seem to eliminate
# autocorr in residuals and cannot be combined with random intercepts in a model
# that converges--remove for now
# age2 <- brm(data = gen_dat, family = Gamma(link = "log"),
#             gen_length ~ 1 + avg_weight + 
#               ar(time = brood_year, gr = stock, p = 1),
#             prior = c(prior(normal(0, 2), class = "Intercept"),
#                       prior(normal(0, 2), class = "b"),
#                       prior(gamma(0.01, 0.01), class = "shape")),
#             iter = 1250, warmup = 500, chains = 4, cores = 4, thin = 10,
#             seed = 13,
#             control = list(adapt_delta = 0.9))
# posterior_summary(age2)
# acf(residuals(age2)[, "Estimate"])


age3 <- brm(data = gen_dat, family = Gamma(link = "log"),
            gen_length ~ 1 + avg_weight + (1 | stock),
            prior = c(prior(normal(0, 2), class = "Intercept"),
                      prior(normal(0, 2), class = "b"),
                      prior(cauchy(0, 1), class = "sd"),
                      prior(gamma(0.01, 0.01), class = "shape")
                      ),
            iter = 1250, warmup = 500, chains = 4, cores = 4,
            seed = 13,
            control = list(adapt_delta = 0.95))
posterior_summary(age3)
acf(residuals(age3)[, "Estimate"])

age4 <- brm(data = gen_dat, family = Gamma(link = "log"),
            gen_length ~ 1 + avg_weight + smolt + avg_weight:smolt + (1 | stock),
            prior = c(prior(normal(0, 2), class = "Intercept"),
                      prior(normal(0, 2), class = "b"),
                      prior(cauchy(0, 2), class = "sd"),
                      prior(gamma(0.01, 0.01), class = "shape"),
                      prior(lkj(2), class = cor)
            ),
            iter = 1000, warmup = 500, chains = 4, cores = 4,
            seed = 13,
            control = list(adapt_delta = 0.9))

temp <- lme4::glmer(data = gen_dat, family = Gamma(link = "log"),
           gen_length ~ 1 + avg_weight + (1 | stock))
