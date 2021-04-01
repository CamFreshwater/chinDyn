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
by_raw <- read.csv(here::here("data","salmonData", 
                              "cwt_indicator_surv_sep2020.csv"), 
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
  mutate(survival = as.numeric(survival),
         brood_year = as.numeric(year)) %>% 
  select(brood_year, stock, stock_name, survival) %>% 
  arrange(stock) 

# mean generation length data
gen1 <- read.csv(here::here("data/salmonData/cwt_indicator_generation_time.csv")) %>% 
  mutate(stock = as.factor(Stock)) %>% 
  select(stock, brood_year = BY, 
         gen_length = GenTim.fishing.mortality.represented.in.calcuations.)


## Generate metadata using old survival as a template than modify in excel
# import old survival data to add some features
# also includes some stocks added manually
# old_surv <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"),
#                      stringsAsFactors = FALSE) %>% 
#   mutate(lat = as.numeric(lat),
#          long = as.numeric(long)) %>% 
#   select(stock, jurisdiction:long) %>%
#   distinct()

# export temporary .csv to fill in metadata
# temp_out <- by_dat1 %>%
#   left_join(., old_surv, by = "stock") %>%
#   # filter(is.na(jurisdiction)) %>%
#   select(stock, stock_name, jurisdiction:long) %>%
#   distinct() %>% 
#   arrange(stock)
# write.csv(temp_out, here::here("data", "salmonData", "metadata.csv"))


# import version cleaned by hand (added lat/longs and two systems HOK and SMK)
metadata <- read.csv(here::here("data", "salmonData", "metadata_clean.csv")) 

by_dat <- metadata %>% 
  left_join(.,
            expand.grid(brood_year = unique(by_dat1$brood_year),
                        stock = unique(metadata$stock)),
            by = "stock") %>% 
  left_join(., by_dat1, by = c("stock", "stock_name", "brood_year")) %>% 
  full_join(., gen1, by = c("stock", "brood_year")) %>%
  mutate(
    lat = as.numeric(hatch_lat),
    long = as.numeric(hatch_long),
    M = -log(survival),
    # change Elwha's region given catch dist similar to Puget
    region = ifelse(stock == "ELW", "NPGSD", region),
    j_group4 = case_when(
      grepl("COLR", region) ~ "col",
      region %in% c("HOODC", "SPGSD", "NPGSD") ~ "puget",
      (lat > 52 & !region == "UFR") ~ "north",
      stock_name == "Transboundary Rivers" ~ "north",
      region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                    "WCVI") ~ "south",
      TRUE ~ "sog"
    ),
    j_group3 = case_when(
      j_group4 %in% c("south", "col") ~ "south",
      TRUE ~ j_group4
    ),
    j_group2 = case_when(
      j_group3 %in% c("puget", "sog") ~ "salish",
      TRUE ~ j_group3
    ),
    j_group1 = case_when(
      j_group2 %in% c("south", "north") ~ "shelf",
      TRUE ~ j_group2
    ),
    a_group4 = case_when(
      #subset of ECVI stocks are north-migrating
      stock_name %in% c("Puntledge River Summer", "Quinsam River Fall",
                        "Lyons Ferry Yearling", "Willamette Spring", 
                        "Atnarko Yearling", "Kitsumkalum Yearling") ~ "north",
      stock_name %in% c("Big Qualicum River Fall", "Chilliwack River Fall",
                        "Cowichan River Fall", "Nanaimo River Fall", 
                        "Harrison River") ~ "sog",
      j_group4 == "puget" ~ "puget",
      region == "LCOLR" ~ "broad",
      smoltType == "streamtype" ~ "offshore",
      TRUE ~ "north"
    ),
    a_group3 = case_when(
      a_group4 %in% c("sog", "puget") ~ "south",
      TRUE ~ a_group4
    ),
    a_group2 = case_when(
      a_group3 %in% c("south", "broad") ~ "south",
      TRUE ~ a_group3
    ),
    a_group1 = case_when(
      a_group2 %in% c("south", "north") ~ "shelf",
      TRUE ~ a_group2
    ),
    run = tolower(adultRunTiming)
  ) %>% 
  select(stock, stock_name, brood_year, survival, M, gen_length,
         jurisdiction, smolt = smoltType, run, 
         region:long, j_group1, j_group2, j_group3, j_group4, a_group1, 
         a_group2, a_group3, a_group4) %>% 
  mutate_at(vars(stock), list(~ factor(., levels = unique(.)))) %>% 
  mutate(year = ifelse(smolt == "streamtype", brood_year + 2, 
                       brood_year + 1),
         smolt = as.factor(smolt),
         j_group4 = as.factor(j_group4),
         j_group3 = as.factor(j_group3),
         j_group2 = as.factor(j_group2),
         j_group1 = as.factor(j_group1),
         j_group4b  = as.factor(paste(j_group4, smolt, sep = "_")),
         j_group3b  = as.factor(paste(j_group3, smolt, sep = "_")),
         j_group2b  = as.factor(paste(j_group2, smolt, sep = "_")),
         j_group1b = as.factor(paste(j_group1, smolt, sep = "_"))
         ) %>% 
  #add final category that separates SoG/PS for subyearlings, but not 
  #yearlings
  mutate(j_group5b = case_when(
    j_group2 == "salish" & smolt == "streamtype" ~ "salish_streamtype",
    TRUE ~ as.character(j_group3b)
  ) %>% 
    as.factor()
  )

saveRDS(by_dat,
        here::here("data", "salmonData", "cwt_indicator_surv_clean.RDS"))

groupings_table <- by_dat %>% 
  select(stock, stock_name, smolt, run, region, j_group1:j_group4, 
         a_group1:a_group4) %>% 
  distinct() %>% 
  arrange(juvenile_grouping)
write.csv(groupings_table, here::here("data", "salmonData", 
                                      "groupings_table.csv"),
          row.names = FALSE)

# How many stocks per region?
by_dat %>% 
  group_by(j_group) %>% 
  summarize(nStocks = length(unique(stock)))
by_dat %>% 
  group_by(a_group) %>% 
  summarize(nStocks = length(unique(stock)))


#-------------------------------------------------------------------------------

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

