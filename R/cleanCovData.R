## cleanCovData.R
# April 24, 2019
# Script to clean data used as covariates when modeling CK survival 
# 1. Harbor seal population abundance from Strahan Tucker (rawData/sealdata); 
# input data are estimated abundances from population growth models using the 
# "standard" correction factor
# 2. Zooplankton data from Ian Perry (rawData/zpData)
# 3. Juvenile stomach contents data
# -----

library(tidyverse); library(here); library(ggplot2); library(viridis)
  
seals <- read.csv(here("data/salmonData/sealPopEst.csv"), 
                  stringsAsFactors = FALSE) %>% 
  rename(year = "Estimate", mean = "Mean", low = "X2.5th", up = "X97.5th", reg = "Region")
sogZP <- read.csv(here("data/salmonData/totalPreyAnomalies_SOG.csv"),
                  stringsAsFactors = FALSE)
viZP <- read.csv(here("data/salmonData/totalPreyAnomalies_sVI.csv"),
                 stringsAsFactors = FALSE)
diet <- read.csv(here("data/salmonData/ckSummerDiet.csv"), 
                 stringsAsFactors = FALSE)
#-----
### Plot each TS

## Diet data
plotDiet <- diet %>% 
  gather(key = "group", value = "ppn", -year, -season) %>% 
  mutate(ppn = case_when(
    is.na(ppn) ~ 0,
    TRUE ~ ppn
  ))

ggplot(plotDiet, aes(x = year, y = ppn, fill = group)) +
  geom_area(position = "stack") +
  scale_fill_viridis_d() +
  samSim::theme_sleekX()

trimDiet <- diet %>% 
  mutate(fishDietAnom = scale(teleostei)) %>% 
  select(year, fishDietAnom)
  
## ZP data
trimSoG <- sogZP %>% 
  mutate(zpEnvAnom = scale(TotalZooplBiomAnom),
         fishEnvAnom = scale(TotalFishBiomAnom)) %>% 
  select(year = Year, zpEnvAnom, fishEnvAnom)
trimEnv <- viZP %>% 
  mutate(borealAnom = scale(boreal)) %>% 
  select(year = Year, borealAnom) %>% 
  left_join(trimSoG, by = "year")

plotEnv <- trimEnv %>% 
  gather(key = "var", value = "anomaly", -Year)

ggplot(plotEnv, aes(x = Year, y = anomaly, colour = var)) +
  geom_line() +
  scale_colour_viridis_d() +
  samSim::theme_sleekX()


## seal data
ggplot(seals, aes(x = year, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.2) +
  samSim::theme_sleekX() +
  facet_wrap(~reg, scales = "free_y")

trimSeal <- seals %>%
  filter(reg == "SOG") %>%
  mutate(sealAnom = scale(mean)) %>% 
  select(year, sealAnom)

## consolidate all
covOut <- trimDiet %>% 
  left_join(trimEnv, by = "year") %>% 
  left_join(trimSeal, by = "year")

write.csv(covOut, here("data/salmonData/survCovariateAnom.csv"), row.names = F)
