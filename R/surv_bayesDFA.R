## Fit Bayesian DFA models to instantaneous mortality data
# July 20, 2020
# Data fed to chinBayesDFA_new.Rmd for distribution to co-authors; updated
# version of fitsBayesDFA_old.R

library(tidyverse)
library(bayesdfa)
library(rstan)

#helper functions to fit and post-process DFA
source(here("R/functions/dfaFunctions.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

scale_f <- function(x) {
  mn <- mean(x, na.rm = T)
  sd <- sd(x, na.rm = T)
  (x - mn) / sd
}

surv <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                 stringsAsFactors = FALSE) %>% 
  mutate(
    M = -log(surv),
    agg_reg = case_when(
      (is.na(lat)) ~ "north",
      (lat > 52 & !region == "UFR") ~ "north",
      (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                     "WCVI")) ~ "south",
      TRUE ~ "SS"
    ),
    group = paste(agg_reg, smoltType, sep = "_")
  ) %>% 
  arrange(group) %>% 
  mutate_at(vars(stock), list(~ factor(., levels = unique(.)))) %>% 
  rename(year = OEY) %>% 
  group_by(stock) %>% 
  mutate(Mz = scale_f(M)) %>% 
  ungroup()
# saveRDS(surv, here::here("data", "salmonData", "clean_dfa_inputs.rds"))

# plots 
ggplot(surv %>% filter(!is.na(M))) +
  geom_point(aes(x = year, y = Mz, fill = group), shape = 21) + 
  facet_wrap(~stock) +
  theme(legend.position = "top") +
  labs(y = "Scaled M")
ggplot(surv %>% filter(!is.na(M))) +
  geom_histogram(aes(x = M, fill = group)) + 
  facet_wrap(~stock) +
  theme(legend.position = "top")


#focus on subset of BC pops for initial analyses
make_mat <- function(x) {
  x %>%
    select(year, stock, Mz) %>%
    spread(key = stock, value = Mz) %>%
    select(-year) %>% 
    as.matrix() %>% 
    t()
}
m_mat <- make_mat(surv)


## FIT GLOBAL MODEL ------------------------------------------------------------

# n trends = to n groups
global_fit6 <- fit_dfa(y = m_mat, num_trends = length(unique(surv$group)), 
                       zscore = FALSE, iter = 3000, chains = 3, thin = 1,
                       control = list(adapt_delta = 0.95, max_treedepth = 20))
r_global <- rotate_trends(global_fit6)

pdf(here::here("figs", "dfa", "bayes", "global_6.pdf"))
plot_trends(r_global)
plot_fitted(global_fit6)
plot_loadings(r_global) +
  ylim(-3, 3)
dev.off()


## GROUP MODEL SELCTION --------------------------------------------------------

surv_trim2 <- surv %>%
  filter(!group %in% c("north_oceantype", "south_streamtype")) 
surv_tbl <- tibble(group = unique(surv_trim2$group)) %>% 
  mutate(
    m_mat = surv_trim2 %>% 
      group_split(group) %>% 
      map(., make_mat)
  )
tt <- surv_tbl[1:2,]

tt$trend_select_list <- lapply(tt$m_mat, function(x) {
  find_dfa_trends(x, iter = 2500, kmin = 1, kmax = 3, chains = 4,
                  variance =  "equal",
                  control = list(adapt_delta = 0.95, max_treedepth = 20))
})

xx <- find_dfa_trends(tt$m_mat[[2]], iter = 1000, kmin = 1, kmax = 2, chains = 1,
                variance =  "equal",
                control = list(adapt_delta = 0.95, max_treedepth = 20))
