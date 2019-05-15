library(ggplot2)
library(PBSmapping)
library(maps)
library(ggmap)
library(mapdata)


stockDat <- read.csv("C:/github/chinDyn/data/salmonData/CLEANcwtInd_age2SR_OEY.csv", 
         stringsAsFactors = FALSE) %>% 
  select(-OEY, - surv) %>% 
  distinct() %>% 
  filter(!lat == "na") %>% 
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         aggReg = case_when(
           (lat > 52 & !region == "UFR") ~ "north",
           (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                         "WCVI")) ~ "south",
           TRUE ~ "SS"
         )) %>% 
  mutate(grp = paste(smoltType, aggReg, sep = "_"))

na <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))

ggplot(data = na, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-136, -119), ylim = c(42.5, 59.5), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray") +
  geom_jitter(data = stockDat, aes(x = long, y = lat, color = as.factor(grp)), 
             inherit.aes = FALSE, width = 0.25, height = 0.5)

