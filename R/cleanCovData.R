## cleanCovData.R
# April 24, 2019
# Script to clean data used as covariates when modeling CK survival 
# 1. Harbor seal population abundance from Strahan Tucker (rawData/sealdata); 
# input data are estimated abundances from population growth models using the 
# "standard" correction factor
# 2. Zooplankton data from Ian Perry (IOS) and Cheryl Morgan (NOAA Newport) 
# (rawData/zpData)
# 3. Juvenile stomach contents data
# 4. Various basin-scale oceanographic drivers (PDO, SSTarc, NPGO)
# -----

library(tidyverse); library(here); library(ggplot2); library(viridis)
  
# load and initial clean
seals <- read.csv(here("data/salmonData/sealPopEst.csv"), 
                  stringsAsFactors = FALSE) %>% 
  rename(year = "Estimate", mean = "Mean", low = "X2.5th", up = "X97.5th", 
         reg = "Region")

sogZP <- read.csv(here("data/salmonData/totalPreyAnomalies_SOG.csv"),
                  stringsAsFactors = FALSE)

viZP <- read.csv(here("data/salmonData/totalPreyAnomalies_sVI.csv"),
                 stringsAsFactors = FALSE)

newportCope <- read.csv(here("data/salmonData/copeIndices_newport.csv"),
                        stringsAsFactors = FALSE)

diet <- read.csv(here("data/salmonData/ckSummerDiet.csv"), 
                 stringsAsFactors = FALSE)

npgo <- read.csv(here("data/salmonData/npgo_jul2020.csv"), 
                 stringsAsFactors = FALSE) %>% 
  rename(year = YEAR, month = MONTH, npgo = NPGO.index)

pdo1 <- read.csv(here("data/salmonData/pdo_aug2020.csv"), 
                 stringsAsFactors = FALSE) 
pdo_date_list <- strsplit(as.character(pdo1$Date), "") 
pdo <- pdo1 %>% 
  mutate(
    year = map(pdo_date_list, function (x) paste0(x[1:4], sep = "", 
                                                  collapse = "")) %>% 
      unlist() %>% 
      as.numeric(),
    month = map(pdo_date_list, function (x) paste0(x[5:6], sep = "", 
                                                   collapse = "")) %>% 
      unlist() %>% 
      as.numeric()
  ) %>% 
  select(year, month, pdo = Value)

sstArc <- read.csv(here("data/salmonData/johnstone_indicators.csv"), 
                 stringsAsFactors = FALSE) %>% 
  select(year = Time, sst_arc = SSTarc)

biIndex <- read.csv(here("data/salmonData/bifurcation-index.csv"), 
                   stringsAsFactors = FALSE) 

#-------------------------------------------------------------------------------

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
  mutate(fish_diet_anom = scale(teleostei)[ , 1]) %>% 
  select(year, fish_diet_anom) 

## ZP data
trimNewport <- newportCope %>% 
  filter(Month %in% seq(5, 9, by = 1)) %>% 
  group_by(Year) %>% 
  mutate(meanCCI = mean(CCI)) %>% 
  ungroup() %>% 
  select(year = Year, meanCCI) %>% 
  distinct() %>% 
  mutate(cci_anom = scale(meanCCI)[ , 1])

trimSoG <- sogZP %>% 
  mutate(zp_env_anom = scale(TotalZooplBiomAnom)[ , 1],
         fish_env_anom = scale(TotalFishBiomAnom)[ , 1]) %>% 
  select(year = Year, zp_env_anom, fish_env_anom)

trimEnv <- viZP %>% 
  mutate(boreal_anom = scale(boreal)[ , 1],
         south_anom = scale(southern)[ , 1]) %>% 
  select(year = Year, boreal_anom, south_anom) %>% 
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
  mutate(seal_anom = scale(mean)[ , 1]) %>% 
  select(year, seal_anom)

# for oceanographic indices which have monthly values calculate annual means and
# means for April-July (priming period)
sstArc2 <- sstArc %>%
  mutate(year_n = floor(year),
         year_c = as.character(year_n),
         month_cont = (year - year_n) * 12,
         month = ceiling(month_cont)
         )
phys_ocean <- left_join(pdo, npgo, by = c("month", "year")) %>% 
  left_join(., sstArc2 %>% select(year = year_n, month, sst_arc), 
            by = c("month", "year")) %>% 
  pivot_longer(., cols = c("pdo", "npgo", "sst_arc"), names_to = "index",
               values_to = "value") %>% 
  #identify months making up priming period
  mutate(prime_period = ifelse(month > 3 & month < 8, "y", "n")) %>% 
  group_by(year, prime_period, index) %>% 
  mutate(prime_mean = mean(value)) %>% 
  group_by(year, index) %>% 
  mutate(annual_mean = mean(value)) %>% 
  ungroup() %>% 
  mutate(prime_anom = scale(prime_mean)[ , 1],
         annual_anom = scale(annual_mean)[ , 1]) %>% 
  filter(prime_period == "y") %>% 
  select(year, index, prime_mean:annual_anom) %>% 
  distinct()
  
phys_ocean_anom <- phys_ocean %>% 
  select(-prime_mean, -annual_mean) %>% 
  pivot_wider(names_from = index, values_from = c(prime_anom, annual_anom)) %>% 
  glimpse()

## consolidate all anomalies
covOut <- phys_ocean_anom %>% 
  full_join(trimSeal, by = "year") %>% 
  full_join(bi_index %>% select(year, bi_anom = bi_stnd), by = "year") %>% 
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