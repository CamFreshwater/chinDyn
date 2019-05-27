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
# -----


require(tidyverse); require(here)
byDatWide <- read.csv(here("data","salmonData","rawData","cwt2SRData",
                           "cwtInd_age2SR_BY.csv"), 
                  stringsAsFactors = FALSE)
eyDatWide <- read.csv(here("data","salmonData","rawData","cwt2SRData",
                           "cwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)
stockInfo <- read.csv(here("data","salmonData","rawData", "cwt2SRData",
                           "pstID.csv"), 
                      stringsAsFactors = FALSE) %>% 
  #remove some extra columns to keep things clean
  select(-oceanStartAge, -terminalNetAge, -maxAge, -smoltAgeCode,
         -adultRunTimingCode, -country, -comments) %>% 
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         aggReg = case_when(
           (is.na(lat)) ~ "north",
           (lat > 52 & !region == "UFR") ~ "north",
           (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                          "WCVI")) ~ "south",
           TRUE ~ "SS"
         )) %>% 
  filter(aggReg == "SS")



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
      Year > 74 & Year < 99 ~ paste("19", Year, sep = ""),
      Year < 10 ~ paste("200", Year, sep =  ""),
      TRUE ~ paste("20", Year, sep = ""))
    )
  )

write.csv(escDat, here::here("data", "salmonData", "CLEANsalishSea_escData.csv"),
          row.names = FALSE)
