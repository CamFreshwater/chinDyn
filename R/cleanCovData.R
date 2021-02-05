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
# 5. Sea surface temperature and salinity data from Entrance Island lighthouse
# 6. Herring age-2 recruit index from Jaclyn Cleary
# 7. SoG bloom timing from Susan Allen
# -----

library(tidyverse); library(here); library(ggplot2); library(viridis)
  
## Load and initial clean of each dataset
seals <- read.csv(here("data/salmonData/sealPopEst.csv")) %>% 
  rename(year = "Estimate", mean = "Mean", low = "X2.5th", up = "X97.5th", 
         reg = "Region")

sogZP <- read.csv(here("data/salmonData/totalPreyAnomalies_SOG.csv"))

viZP <- read.csv(here("data/salmonData/totalPreyAnomalies_sVI.csv"))

newportCope <- read.csv(here("data/salmonData/copeIndices_newport.csv"))

diet <- read.csv(here("data/salmonData/ckSummerDiet.csv"))

whales <- read.csv(here("data/salmonData/rkw_salmon_jun2019.csv"))

npgo <- read.csv(here("data/salmonData/npgo_jul2020.csv"))
colnames(npgo) <- c("year", "month", "npgo")


pdo1 <- read.csv(here("data/salmonData/pdo_aug2020.csv")) 
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


sstArc <- read.csv(here("data/salmonData/johnstone_indicators.csv")) 
names(sstArc)[1] <- "year"
sstArc <- sstArc %>%  select(year, sst_arc = SSTarc)


biIndex <- read.csv(here("data/salmonData/bifurcation-index.csv")) %>% 
  mutate(bi_stnd = scale(bifurcation_index)[, 1])


# herring data
# helper function to rename herring data
split_bind <- function(x) {
  strsplit(x, "_") %>% 
    map(., .f = function(y) y[1]) %>% 
    unlist(.)
}

herring <- read.csv(here::here("data/salmonData/herring_r.csv")) %>% 
  pivot_wider(., names_from = "Indicator", values_from = "Value") %>% 
  pivot_longer(cols = contains("_Med"), names_to = "stock", 
               values_to = "median") %>%
  pivot_longer(cols = contains("_LowerCI"), names_to = "stock2", 
               values_to = "lo_ci") %>%
  pivot_longer(cols = contains("_UpperCI"), names_to = "stock3", 
               values_to = "up_ci") %>% 
  mutate(
    stock = split_bind(stock),
    stock2 = split_bind(stock2), 
    stock3 = split_bind(stock3)
  ) %>% 
  #retain only rows with stocks matching
  filter(
    (stock == stock2 & stock2 == stock3)
  ) %>%
  select(year = Year, stock, median, lo_ci, up_ci) 


# Entrance Island data
sog_sst_wide <- read.csv(here::here("data/salmonData/entrance_island_sst.csv"))
colnames(sog_sst_wide) <- tolower(colnames(sog_sst_wide))
sog_salinity_wide <- read.csv(here::here("data/salmonData/entrance_island_salinity.csv"))
colnames(sog_salinity_wide) <- tolower(colnames(sog_salinity_wide))

sog_sst <- pivot_longer(sog_sst_wide, cols = -year, names_to = "month", 
                        values_to = "entrance_sst") 
sog_salinity <- pivot_longer(sog_salinity_wide, cols = -year, names_to = "month", 
                             values_to = "entrance_salinity") 
sog_ocean <- left_join(sog_sst, sog_salinity, by = c("year", "month")) %>% 
  mutate(month_f = fct_relevel(as.factor(month), 
                             "jan", "feb", "mar", "apr", "may", "jun", "jul", 
                             "aug", "sep", "oct", "nov", "dec"),
         month = as.numeric(month_f),
         #replace filler values
         entrance_sst = ifelse(entrance_sst == "99.99",
                               NA,
                               entrance_sst),
         entrance_salinity = ifelse(entrance_salinity == "99.99",
                               NA,
                               entrance_salinity)) %>% 
  filter(!is.na(entrance_sst),
         !is.na(entrance_salinity))


# Bloom timing data
bloom_path <- here::here("data/salmonData/bloomdates_29mar2019.dat")
# skip first three lines which are metadata
bloom <- read.table(bloom_path, header = TRUE, skip = 3) %>% 
  rename(year = Year, bloom_day = YearDay)
  

#-------------------------------------------------------------------------------

### Plot each TS

## Diet data
# plotDiet <- diet %>% 
#   gather(key = "group", value = "ppn", -year, -season) %>% 
#   mutate(ppn = case_when(
#     is.na(ppn) ~ 0,
#     TRUE ~ ppn
#   )) %>% 
#   ggplot(., aes(x = year, y = ppn, fill = group)) +
#   geom_area(position = "stack") +
#   scale_fill_viridis_d() +
#   samSim::theme_sleekX()

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

# trimEnv %>% 
#   gather(key = "var", value = "anomaly", -year) %>% 
#   ggplot(., aes(x = year, y = anomaly, colour = var)) +
#   geom_line() +
#   scale_colour_viridis_d() +
#   samSim::theme_sleekX()


# seal data
# ggplot(seals, aes(x = year, y = mean)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.2) +
#   samSim::theme_sleekX() +
#   facet_wrap(~reg, scales = "free_y")

trimSeal <- seals %>%
  filter(reg == "SOG") %>%
  mutate(seal_anom = scale(mean)[ , 1]) %>% 
  select(year, seal_anom)


# herring data
# ggplot(herring, aes(x = year, y = median)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lo_ci, ymax = up_ci), alpha = 0.4) +
#   samSim::theme_sleekX() +
#   facet_wrap(~stock, scales = "free_y")

trimHerring <- herring %>% 
  group_by(stock) %>% 
  mutate(herr_anom = scale(median)[ , 1]) %>% 
  ungroup() %>% 
  filter(stock %in% c("WCVIHerringR", "SoGHerringR")) %>% 
  pivot_wider(names_from = stock, values_from = herr_anom) %>% 
  select(year, sog_herr_anom = SoGHerringR, wcvi_herr_anom = WCVIHerringR) 

# resident killer whale data
# whales %>% 
#   pivot_longer(., cols = c(NRKW_N, SRKW_N), names_to = "population", 
#                values_to = "abundance") %>% 
#   ggplot(.) +
#   geom_line(aes(x = Year, y = abundance, color = population)) 

trimWhales <- whales %>% 
  #scale and subtract two years to estimate effects at maturity
  mutate(lagged_year = Year - 2,
         nrkw_anom = scale(NRKW_N)[ , 1],
         srkw_anom = scale(SRKW_N)[ , 1],
         rkw_anom = scale(NRKW_N + SRKW_N)[ , 1]) %>% 
  select(year = lagged_year, nrkw_anom, srkw_anom, rkw_anom)

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
  left_join(., sog_ocean, by = c("month", "year")) %>% 
  pivot_longer(., cols = c("pdo", "npgo", "sst_arc", "entrance_sst",
                           "entrance_salinity"), 
               names_to = "index",
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
covWide <- phys_ocean_anom %>% 
  full_join(trimSeal, by = "year") %>% 
  full_join(trimWhales, by = "year") %>% 
  full_join(trimHerring, by = "year") %>% 
  full_join(biIndex %>% select(year, bi_anom = bi_stnd), by = "year") %>% 
  full_join(trimDiet, by = "year") %>% 
  full_join(trimEnv, by = "year") 

covOut <- covWide %>% 
  pivot_longer(cols = 2:ncol(.), names_to = "metric", values_to = "anomaly") %>%
  mutate(
    time_step = case_when(
      grepl("prime", metric) ~ "seasonal",
      metric == "cci_anom" ~ "seasonal",
      TRUE ~ "annual"
    ),
    class = case_when(
      grepl("annual", metric) ~ "physical",
      grepl("prime", metric) ~ "physical",
      metric == "bi_anom" ~ "physical",
      grepl("seal", metric) ~ "predator",
      grepl("rkw", metric) ~ "predator",
      TRUE ~ "diet"
    ),
    region = case_when(
      metric == "srkw_anom" ~ "sog",
      grepl("rkw", metric) ~ "basin",
      class == "physical" ~ "basin",
      metric == "cci_anom" ~ "ca_current",
      metric %in% c("zp_env_anom", "fish_env_anom", "seal_anom",
                    "sog_herr_anom") ~ "sog",
      metric %in% c("boreal_anom", "south_anom", "fish_diet_anom",
                    "wcvi_herr_anom") ~ "wcvi"
    )
  )
  
write.csv(covOut, here("data/salmonData/survCovariateAnom.csv"), row.names = F)


## Export whale, seal and herring data as list for juvenile analyses (i.e. SoG
# only, herring lagged by two years)

# Eventually incorporate uncertainty using CIs and assuming:
# SD = sqrt(n-1) x ((up - lo) / 3.92) where n is length of time series for seals
# (should be DF of model) and number of MCMC iterations for herring estimates
# herr_clean <- herring %>% 
#   mutate(sd = sqrt(1000) * (abs(lo_ci - up_ci) / 3.92)) %>% 
#   glimpse()
# seal_clean <- seals %>% 
#   filter(reg == "SOG") %>% 
#   mutate(n = length(unique(year)),
#          sd = sqrt(n - 1) * (abs(low - up) / 3.92))

juv_cov <- herring %>% 
  filter(stock == "SoGHerringR") %>%
  #subtract by two since index is for age-2 recruits and age-0 likely most 
  #important
  mutate(year2 = year - 2) %>%
  select(herring_model_year = year, year = year2, herr_abund = median) %>% 
  left_join(., 
            seals %>%
              filter(reg == "SOG") %>% 
              select(year, seal_abund = mean),
            by = "year") %>% 
  left_join(.,
            phys_ocean %>% 
              filter(index %in% c("entrance_sst", "entrance_salinity"),
                     !is.na(prime_mean)) %>% 
              pivot_wider(., -c(prime_anom, annual_anom, annual_mean), 
                          names_from = "index", values_from = "prime_mean",
                          names_prefix = "prime_"),
            by = "year") %>% 
  left_join(., bloom, by = "year") 

saveRDS(juv_cov, here::here("data/salmonData/cov_subset_juv.rds"))


# adult cleaning requires more complex stock-specific considerations 
# incorporated in generation_gams
adult_cov <- list(herring = herring %>%
                    mutate(herring_age0_year = year - 2) %>% 
                    select(herring_model_year = year, 
                           herring_age0_year,
                           herr_abund = median, 
                           stock),
                  rkw = whales %>% 
                    mutate(total_n = SRKW_N + NRKW_N) %>% 
                    select(year = Year, srkw_n = SRKW_N, nrkw_n = NRKW_N, 
                           total_n),
                  sog = phys_ocean %>% 
                    filter(index %in% c("entrance_sst", "entrance_salinity"),
                           !is.na(prime_mean)) %>% 
                    pivot_wider(., -c(prime_anom, annual_anom, annual_mean), 
                                names_from = "index", values_from = "prime_mean",
                                names_prefix = "prime_") %>% 
                    left_join(., bloom, by = "year") 
)

saveRDS(adult_cov, here::here("data/salmonData/cov_subset_adult.rds"))


## EXPLORATORY -----------------------------------------------------------------

## Looks at covariates
cov_ts <- covOut %>% 
  ggplot(.) +
  geom_line(aes(x = year, y = anomaly, colour = class)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~metric)

png(here::here("figs", "dfa", "covEffects", "cov_anomalies.png"), height = 5,
    width = 7, units = "in", res = 300)
cov_ts
dev.off()


## Correlelogram 
corCov <- covWide %>% 
  select(-year) %>% 
  na.omit() %>% 
  as.matrix() %>% 
  cor()
corrplot::corrplot(corCov, method = "number", "upper")


## Length of time series
covOut %>% 
  filter(!is.na(anomaly)) %>% 
  group_by(metric) %>% 
  summarize(min = min(year),
            max = max(year))
