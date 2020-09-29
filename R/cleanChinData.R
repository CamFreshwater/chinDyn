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


# import old survival data to add some features
old_surv <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                     stringsAsFactors = FALSE)
# new survival data
by_raw <- read.csv(here("data","salmonData","rawData",
                           "cwt_indicator_surv_sep2020.csv"), 
                      stringsAsFactors = FALSE)
stock_key <- data.frame(
  stock = colnames(by_raw)[2:58],
  stock_name = t(by_raw[1, 2:58])
) %>% 
  rename(stock_name = X1)
by_dat1 <- by_raw[-1, ] %>% 
  pivot_longer(., 2:ncol(.), names_to = "stock", values_to = "survival") %>% 
  left_join(., stock_key, by = "stock") %>% 
  mutate(survival = as.numeric(survival)) %>% 
  select(year = X, stock, stock_name, survival) %>% 
  arrange(stock)

# export temporary .csv to fill in metadata
# temp_out <- by_dat1 %>%
#   left_join(., old_surv %>% select(stock, jurisdiction:long) %>% distinct(), 
#             by = "stock") %>%   
#   # filter(is.na(jurisdiction)) %>% 
#   select(stock_name, jurisdiction:long) %>% 
#   distinct()
# write.csv(temp_out, here::here("data", "salmonData", "metadata.csv"))

metadata <- read.csv(here("data", "salmonData", "metadata_clean.csv")) %>% 
  select(-X)

by_dat <- by_dat1 %>% 
  left_join(., metadata, by = "stock_name") %>%
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         M = -log(survival),
         agg_reg = case_when(
           (is.na(lat)) ~ "north",
           (lat > 52 & !region == "UFR") ~ "north",
           (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                          "WCVI")) ~ "south",
           TRUE ~ "SS"
         ),
         run = tolower(adultRunTiming),
         group = paste(agg_reg, smoltType, sep = "_")) %>% 
  select(year:stock_name, smolt = smoltType, run, 
         region:long, agg_reg, group, survival, M)

write.csv(by_dat, here("data", "salmonData", "cwt_indicator_surv_clean.csv"))


#How many stocks per region?
by_dat %>% 
  group_by(region) %>% 
  summarize(nStocks = length(unique(stock)))
by_dat %>% 
  group_by(group) %>% 
  summarize(nStocks = length(unique(stock)))


## OLD ##
# byDatWide <- read.csv(here("data","salmonData","rawData","cwt2SRData",
#                            "cwtInd_age2SR_BY.csv"), 
#                   stringsAsFactors = FALSE)
# eyDatWide <- read.csv(here("data","salmonData","rawData","cwt2SRData",
#                            "cwtInd_age2SR_OEY.csv"), 
#                   stringsAsFactors = FALSE)
# stockInfo <- read.csv(here("data","salmonData","rawData", "cwt2SRData",
#                            "pstID.csv"), 
#                       stringsAsFactors = FALSE) %>% 
#   #remove some extra columns to keep things clean
#   select(-oceanStartAge, -terminalNetAge, -maxAge, -smoltAgeCode,
#          -adultRunTimingCode, -country, -comments) %>% 
#   mutate(lat = as.numeric(lat),
#          long = as.numeric(long),
#          aggReg = case_when(
#            (is.na(lat)) ~ "north",
#            (lat > 52 & !region == "UFR") ~ "north",
#            (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
#                           "WCVI")) ~ "south",
#            TRUE ~ "SS"
#          )) %>% 
#   filter(aggReg == "SS")
# 
# 
# 
# #Lengthen each dataset
# byDat <- byDatWide %>% 
#   gather(key = "stock", value = "surv", -BY) 
# eyDat <- eyDatWide %>% 
#   gather(key = "stock", value = "surv", -OEY)
# 
# # Function to merge and filter each dataset w/ less than 10 years of survival
# # data
# cleanPST <- function(dat, stockInfo) {
#   #ID stocks w/ long TS
#   longStks <- dat %>% 
#     filter(!is.na(surv)) %>% 
#     group_by(stock) %>% 
#     summarize(tsLength = length(surv)) %>% 
#     filter(!tsLength <= 10) 
#   #join and subset
#   dat %>% 
#     left_join(stockInfo) %>% 
#     filter(stock %in% longStks$stock)
# }
# byDat <- cleanPST(byDat, stockInfo)
# write.csv(byDat, here("data", "salmonData", "CLEANcwtInd_age2SR_BY.csv"), 
#           row.names = FALSE)
# eyDat <- cleanPST(eyDat, stockInfo)
# write.csv(eyDat, here("data", "salmonData", "CLEANcwtInd_age2SR_OEY.csv"), 
#           row.names = FALSE)


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
