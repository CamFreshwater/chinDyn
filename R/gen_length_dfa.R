## Fit MARSS then Bayes DFA models to mean generation time data
# Jan 8 2021
# Update on surv_bayesDFA; instead of deciding groupings a priori, first perform
# model selection for groupings with MARSS (faster to converge than bayesDFA), 
# then fit group specific models with bayesdfa

library(MARSS)
library(tidyverse)

stk_data <- readRDS(here::here("data/salmonData/cwt_indicator_surv_clean.RDS")) %>% 
  select(brood_year, stock, stock_name, smolt, run, region, a_group, a_group2, 
         year)

gen <- read.csv(here::here("data/salmonData/cwt_indicator_generation_time.csv")) %>% 
  mutate(stock = as.factor(Stock)) %>% 
  select(stock, brood_year = BY, 
         gen_length = GenTim.fishing.mortality.represented.in.calcuations.) %>% 
  # add relevant metadata 
  left_join(., stk_data, by = c("brood_year", "stock")) %>% 
  mutate(a_group =  as.factor(a_group),
         a_group2 = as.factor(a_group2))

gen %>% 
  # select(a_group, a_group2) %>% 
  filter(stock %in% c("STI", "TAK")) %>%
  # select(stock_name) %>%
  distinct()
# TST combines STI/TAK and AKS combines SSA and NSA


## EXPLORATORY -----------------------------------------------------------------

#plot raw generation length data
gen %>%
  filter(!is.na(gen_length),
         !is.na(stock_name)) %>%
  ggplot(.) +
  geom_point(aes(x = year, y = gen_length, fill = a_group2), shape = 21) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group2))) +
  theme(legend.position = "top") +
  labs(y = "Mean Generation Length") +
  ggsidekick::theme_sleek()
