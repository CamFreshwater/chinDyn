## cleanChinData.R
# April 24, 2019
# Script to clean and convert to long format data from following data
# 1. Age-2 survival data from coastwide CWT stocks by brood year (copied from 
#    SurvivalRateStats_AllCWTindicators_ERA2017 - Copy for Strahan.xlsx)
# 2. Age-2 survival data from coastwide CWT stocks by entry year (copied from 
#    SurvivalRateStats_AllCWTindicators_ERA2017 - Copy for Strahan.xlsx)
# 3. CWT stock information from PST (copied from 
#    SurvivalRateStats_AllCWTindicators_ERA2017 - Copy for Strahan.xlsx)
# 4. CTC raw escapement data from all stocks (mix of hatchery and wild) (copied
#    from EscGraphs2018_Master_4.5.19.xlsx)
# Updated Sep 28 - switch survival data to updated version provided by C. Parken
# (original file surv rate stats raw_2020 analysis.xlsx)
# -----

require(tidyverse); require(here)


# new survival data
by_raw <- read.csv(here::here("data","salmonData", "cwt_indicator_surv_sep2020.csv"), 
                      stringsAsFactors = FALSE)
colnames(by_raw)[1] <- "year"
stock_key <- data.frame(
  stock = colnames(by_raw)[2:58],
  stock_name = t(by_raw[1, 2:58])
) %>% 
  rename(stock_name = X1)
by_dat1 <- by_raw[-1, ] %>% 
  pivot_longer(., 2:ncol(.), names_to = "stock", values_to = "survival") %>% 
  left_join(., stock_key, by = "stock") %>% 
  mutate(survival = as.numeric(survival)) %>% 
  select(brood_year = year, stock, stock_name, survival) %>% 
  arrange(stock) %>% 
  #remove stocks that are aggregates of others on CP's advice
  filter(!stock %in% c("TST", "AKS"))

# import old survival data to add some features
# old_surv <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
#                      stringsAsFactors = FALSE)
# export temporary .csv to fill in metadata
# temp_out <- by_dat1 %>%
#   left_join(., old_surv %>% select(stock, jurisdiction:long) %>% distinct(), 
#             by = "stock") %>%   
#   # filter(is.na(jurisdiction)) %>% 
#   select(stock_name, jurisdiction:long) %>% 
#   distinct()
# write.csv(temp_out, here::here("data", "salmonData", "metadata.csv"))

metadata <- read.csv(here::here("data", "salmonData", "metadata_clean.csv")) %>% 
  select(-X)

by_dat <- by_dat1 %>% 
  left_join(., metadata, by = "stock_name") %>%
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         M = -log(survival),
         j_group = case_when(
           (lat > 52 & !region == "UFR") ~ "north",
           region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                          "WCVI") ~ "south",
           TRUE ~ "salish"
         ),
         a_group = case_when(
           smoltType == "streamtype" ~ "offshore",
           #subset of ECVI stocks are north-migrating
           stock_name %in% c("Puntledge River Summer", "Quinsam River Fall") ~
             "north",
           region %in% c("ECVI", "LFR", "HOODC", "SPGSD", "NPGSD") ~ "south",
           grepl("COLR", region) ~ "columbia",
           TRUE ~ "north"
         ),
         run = tolower(adultRunTiming),
         j_group2 = paste(j_group, smoltType, sep = "_")
         ) %>% 
  select(brood_year:stock_name, smolt = smoltType, run, 
         region:long, j_group, j_group2, a_group, 
         survival, M)

write.csv(by_dat, 
          here::here("data", "salmonData", "cwt_indicator_surv_clean.csv"),
          row.names = FALSE)

a_palette <- disco::disco("bright", n = length(unique(by_dat$a_group)))
names(a_palette) <- unique(by_dat$a_group)
j_palette <- disco::disco("muted", n = length(unique(by_dat$j_group)))
names(j_palette) <- unique(by_dat$j_group)
pals <- list(a_palette, j_palette)
saveRDS(pals, here::here("data", "color_pals.RDS"))

#How many stocks per region?
by_dat %>% 
  group_by(j_group) %>% 
  summarize(nStocks = length(unique(stock)))
by_dat %>% 
  group_by(a_group) %>% 
  summarize(nStocks = length(unique(stock)))

#----- 

## Prep escapement data
# Focus only on Salish Sea stocks
escDatWide <- read.csv(here("data", "salmonData", "rawData", "escapementData",
                            "escDataWideTrim.csv"), 
                       stringsAsFactors = FALSE)

escDat <- escDatWide %>% 
  gather(key = stock, value = esc, -Year) %>% 
  #restrict to Salish Sea only stocks
  filter(stock %in% c("Upper.Geo.", "Cowichan", "Nanaimo", "Fraser..Sp.1.3",
                      "Fraser..Sp.1.2", "Nicola.Sp.1.2", "Fraser..Sum.1.3", 
                      "Fraser..Sum.0.3", "L..Shuswap", "Harrison", "Skagit.Spr",
                      "Skagit.Sum", "Stillaguamish", "Snohomish", "Green", 
                      "Nooksack.Spr", "Lake.Washington")) %>% 
  mutate(year = as.numeric(
    case_when(
      Year > 74 & Year < 100 ~ paste("19", Year, sep = ""),
      Year < 10 ~ paste("200", Year, sep =  ""),
      TRUE ~ paste("20", Year, sep = ""))
    )
  ) %>% 
  select(year, stock, esc)

write.csv(escDat, here::here("data", "salmonData", "CLEANsalishSea_escData.csv"),
          row.names = FALSE)
