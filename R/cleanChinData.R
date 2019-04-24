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
byDatWide <- read.csv(here("data/salmonData/rawData/cwtInd_age2SR_BY.csv"), 
                  stringsAsFactors = FALSE)
eyDatWide <- read.csv(here("data/salmonData/rawData/cwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)
stockInfo <- read.csv(here("data/salmonData/rawData/pstID.csv"), 
                      stringsAsFactors = FALSE) %>% 
  #remove some extra columns to keep things clean
  select(-oceanStartAge, -terminalNetAge, -maxAge, -smoltAgeCode,
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
write.csv(byDat, here("data", "salmonData", "CLEANcwtInd_age2SR_BY.csv"), 
          row.names = FALSE)
eyDat <- cleanPST(eyDat, stockInfo)
write.csv(eyDat, here("data", "salmonData", "CLEANcwtInd_age2SR_OEY.csv"), 
          row.names = FALSE)

#-----

## Random exploration

#How many stocks per region?
byDat %>% 
  group_by(region) %>% 
  summarize(nStocks = length(unique(stock)))
byDat %>% 
  group_by(jurisdiction) %>% 
  summarize(nStocks = length(unique(stock)))
