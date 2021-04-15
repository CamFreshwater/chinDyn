## Map of approximate hatchery locations
# Updated April 15, 2021

library(tidyverse)
library(maptools)
library(rmapshaper)
library(mapdata)


# import sf dataframe (made in chin_tagging/R/prep_spatial_data.R)
coast_sf <- readRDS(
  here::here("data", "salmon_data", "coast_major_river_sf_plotting_large.RDS"))

# generate base plot of study area
base_map <- ggplot() +
  geom_sf(data = coast_sf, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(fill = "darkgrey")) +
  #removes extra border
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 


# salmon data
gen_raw <- readRDS(here::here("data", "salmon_data", 
                              "cwt_indicator_surv_clean.RDS")) 
dat <- read.csv(here::here("data", "salmon_data", "metadata_clean.csv")) %>% 
  left_join(., gen_raw %>% select(stock, j_group3) %>% distinct(), 
            by = "stock") %>% 
  mutate(
    life_history = fct_recode(smoltType, 
                              Yearling = "streamtype", 
                              Subyearling = "oceantype"),
    juv_region = fct_relevel(j_group3, "north", "puget", "sog", "south"),
    juv_region = fct_recode(juv_region, 
                            "North" = "north", "Puget Sound" = "puget",
                            "Strait of Georgia" = "sog", "South" = "south")
    # grouping = fct_recode(j_group3b,
    #                       "North Subyearling" = "north_oceantype", 
    #                       "North Yearling" = "north_streamtype", 
    #                       "Puget Subyearling" = "puget_oceantype", 
    #                       "Puget Yearling" = "puget_streamtype", 
    #                       "SoG Subyearling" = "sog_oceantype", 
    #                       "SoG Yearling" = "sog_streamtype",
    #                       "South Subyearling" = "south_oceantype",
    #                       "South Yearling" = "south_streamtype"
    # ),
    # grouping = fct_relevel(grouping,
    #                        "North Subyearling", 
    #                        "North Yearling", 
    #                        "Puget Subyearling", 
    #                        "Puget Yearling", 
    #                        "SoG Subyearling", 
    #                        "SoG Yearling",
    #                        "South Subyearling",
    #                        "South Yearling")
  ) %>% 
  filter(!stock %in% c("TST", "AKS"))

# plot with salmon data
fill_pal <- RColorBrewer::brewer.pal(name = "Spectral", 
                                     n = length(unique(dat$juv_region)))
names(fill_pal) <- unique(dat$juv_region)
shape_pal <- c(21, 23)
names(shape_pal) <- c("Subyearling", "Yearling")

oey_map <- base_map +
  geom_point(data = dat,
             aes(x = hatch_long, y = hatch_lat, 
                 fill = juv_region, shape = life_history),
             size = 3,
             inherit.aes = FALSE) +
  scale_fill_manual(name = "Ocean Entry Region",
                    values = fill_pal) +
  scale_shape_manual(name = "Juvenile Life History",
                     values = shape_pal) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


# inset map
w_can <- map_data("world", region = c("usa", "canada")) %>%
  filter(lat > 25 & lat < 80,
         long >-175 & long < -60) %>% 
  fortify(.)
inset_map <-  ggplot() +
  geom_polygon(data = w_can, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "darkgrey") + 
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  geom_rect(aes(xmin = -136.5, xmax = -118, ymin = 42, ymax = 60),
            fill = NA, lty = 2, col = "black") +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "top",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  coord_fixed(ratio = 1.3)


png(here::here("figs", "ms_figs", "stock_locs.png"), 
    height = 7, width = 6, units = "in", res = 300)
cowplot::ggdraw() +
  cowplot::draw_plot(oey_map) +
  cowplot::draw_plot(inset_map, x = 0.13, y = 0.01, width = 0.3, height = 0.35)
dev.off()
