## cleanChinData.R
# April 24, 2019
# Script to clean and convert to long format data from following data
# 1. Age-2 survival data from coastwide CWT stocks by brood year (copied from 
#    SurvivalRateStats_AllCWTindicators_ERA2017 - Copy for Strahan.xlsx)
# 2. Age-2 survival data from coastwide CWT stocks by entry year (copied from 
#    SurvivalRateStats_AllCWTindicators_ERA2017 - Copy for Strahan.xlsx)
# 3. CWT stock information from PST (copied from 
#    SurvivalRateStats_AllCWTindicators_ERA2017 - Copy for Strahan.xlsx)
# -----
require(tidyverse); require(here)
byDatWide <- read.csv(here("data/salmonData/cwtInd_age2SR_BY.csv"), 
                  stringsAsFactors = FALSE)
eyDatWide <- read.csv(here("data/salmonData/cwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)
stockInfo <- read.csv(here("data/salmonData/pstID.csv"), 
                      stringsAsFactors = FALSE) %>% 
  #remove some extra columns to keep things clean
  select(-jurisdiction, -oceanStartAge, -terminalNetAge, -maxAge, -smoltAgeCode,
         -adultRunTimingCode, -country, -comments)

#Lengthen each dataset
byDat <- byDatWide %>% 
  gather(key = "stock", value = "surv", -BY) 
eyDat <- eyDatWide %>% 
  gather(key = "stock", value = "surv", -OEY)

# Function to merge and filter each dataset w/ less than 10 years of survival
# data
cleanPST <- function(dat, stockInfo) {
  #ID stocks w/ long TS
  longStks <- dat %>% 
    filter(!is.na(surv)) %>% 
    group_by(stock) %>% 
    summarize(tsLength = length(surv)) %>% 
    filter(!tsLength <= 10) 
  #join and subset
  dat %>% 
    left_join(stockInfo) %>% 
    filter(stock %in% longStks$stock)
}
byDat <- cleanPST(byDat, stockInfo)
eyDat <- cleanPST(eyDat, stockInfo)
