library(ggplot2)
library(PBSmapping)
library(maps)
library(ggmap)
library(mapdata)


stockDat  <- read.csv(here::here("data", "salmonData", "metadata_clean.csv")) %>% 
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         aggReg = case_when(
           (lat > 52 & !region == "UFR") ~ "north",
           (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                         "WCVI")) ~ "south",
           TRUE ~ "SS"
         )) %>% 
  mutate(grp = as.factor(paste(smoltType, aggReg, sep = "_"))) %>% 
  arrange(desc(lat)) %>% 
  mutate(grp = factor(grp, unique(grp))) %>% 
  filter(!grp %in% c("oceantype_north", "streamtype_south")) %>% 
  mutate(grp = fct_recode(grp, "Yearling\nNorth" = "streamtype_north", 
                          "Yearling\nSalish Sea" = "streamtype_SS", 
                          "Subyearling\nSalish Sea" = "oceantype_SS",
                          "Subyearling\nSouth" = "oceantype_south"))

na <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))

mapP <- ggplot(data = na, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-136, -119), ylim = c(42.5, 59.5), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_jitter(data = stockDat, aes(x = long, y = lat, color = as.factor(grp)), 
             inherit.aes = FALSE, width = 0.1, height = 0.1) + 
  labs(color = "Life History Grouping", x = "Longitude", y = "Latitude") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.2, 0.15))

png(here::here("figs", "hatcheryMap.png"), height = 5.5, width = 6, units = "in",
    res = 300)
print(mapP)
dev.off()
