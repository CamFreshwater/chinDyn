## cleanCovData.R
# April 24, 2019
# Script to clean data used as covariates when modeling CK survival 
# 1. Harbor seal population abundance from Strahan Tucker (rawData/sealdata); 
# input data are estimated abundances from population growth models using the 
# "standard" correction factor
# 2. Zooplankton data from Ian Perry (IOS) and Cheryl Morgan (NOAA Newport) 
# (rawData/zpData)
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
newportCope <- read.csv(here("data/salmonData/copeIndices_newport.csv"),
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
trimNewport <- newportCope %>% 
  filter(Month %in% seq(5, 9, by = 1)) %>% 
  group_by(Year) %>% 
  mutate(meanCCI = mean(CCI)) %>% 
  ungroup() %>% 
  select(year = Year, meanCCI) %>% 
  distinct() %>% 
  mutate(cciAnom = scale(meanCCI)[,1])

trimSoG <- sogZP %>% 
  mutate(zpEnvAnom = scale(TotalZooplBiomAnom),
         fishEnvAnom = scale(TotalFishBiomAnom)) %>% 
  select(year = Year, zpEnvAnom, fishEnvAnom)

trimEnv <- viZP %>% 
  mutate(borealAnom = scale(boreal)) %>% 
  select(year = Year, borealAnom) %>% 
  full_join(trimSoG, by = "year") %>% 
  full_join(trimNewport %>% select(-meanCCI), by = "year")

plotEnv <- trimEnv %>% 
  gather(key = "var", value = "anomaly", -year)

ggplot(plotEnv, aes(x = year, y = anomaly, colour = var)) +
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
covOut <- trimSeal %>% 
  full_join(trimDiet, by = "year") %>% 
  full_join(trimEnv, by = "year")

write.csv(covOut, here("data/salmonData/survCovariateAnom.csv"), row.names = F)


# ------------------------------------------------------------------------------

# raw comparison of seal abundance data with survival

eyDat <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)

dat1 <- eyDat %>% 
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         aggReg = case_when(
           (is.na(lat)) ~ "north",
           (lat > 52 & !region == "UFR") ~ "north",
           (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                          "WCVI")) ~ "south",
           TRUE ~ "SS"
         ),
         group = paste(aggReg, smoltType, sep = "_")) %>% 
  rename(year = OEY) %>% 
  left_join(., 
            seals %>% 
              filter(reg == "SOG") %>% 
              select(mean, year),
            by = "year") %>% 
  filter(!is.na(mean),
         !is.na(surv))

seal_plot <- dat1 %>% 
  filter(aggReg == "SS") %>% 
  ggplot() +
  geom_point(aes(x = mean, y = surv, fill = group), shape = 21) +
  ggsidekick::theme_sleek() +
  facet_wrap(~stock)

pdf(here::here("figs", "seal_corr.pdf"), height = 7, width = 8)
seal_plot
dev.off()