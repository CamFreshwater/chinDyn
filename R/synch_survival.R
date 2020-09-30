## Explore synchrony in CWT survival to age 2 data
# Sep 29, 2020

library(synchrony)
library(tidyverse)

surv <- read.csv(here::here("data/salmonData/cwt_indicator_surv_clean.csv"), 
                 stringsAsFactors = FALSE)

stks <- surv %>% 
  filter(!is.na(M)) %>% 
  group_by(stock) %>% 
  summarize(min_yr = min(year),
            max_yr = max(year),
            .groups = "drop") %>% 
  ungroup() %>% 
  filter(min_yr < 1981,
         max_yr > 2014) %>% 
  pull(stock)

# function to take dataframe, subset and convert to matrix for community.synch() 
make_mat <- function(x) {
  out <- x %>%
    select(year, stock, M) %>%
    spread(key = stock, value = M) %>%
    select(-year) %>% 
    as.matrix() 
  rownames(out) <- unique(x$year)
  return(out)
}


# make matrices consisting of various subgroupings (limited by necessity of
# not having any NAs)
gappy_dat <- surv %>% 
  filter(stock %in% stks,
         !stock %in% c("WRY", "TAK", "STL", "PPS", "LRW"), 
         !year < 1981,
         !year > 2014) 
# matrix with maximum number of stocks 
m_stocks <- gappy_dat %>%
  filter(!year < 1985) %>% 
  make_mat()
# matrix with maximum number of continuous years
m_length <- make_mat(gappy_dat) %>% 
  t(.) %>% 
  .[complete.cases(.), ] %>% 
  t(.)
# matrix with southern ocean types
m_south <- gappy_dat %>% 
  filter(!year < 1985,
         group == "south_oceantype") %>% 
  make_mat()
# matrix with salish sea ocean types
m_salish <- gappy_dat %>% 
  filter(!year < 1985,
         group == "SS_oceantype") %>% 
  make_mat()

# check regional coverage of remaining stocks
life_hist <- surv %>% 
  filter(stock %in% row.names(m_stocks)) %>% 
  select(stock_name, smolt, region, group) %>% 
  distinct() %>% 
  arrange(group)

# combine into tibble and estimate synchrony
synch_tbl <- tibble(group = c("coastwide", "coastwide_long", "south_subyearling",
                              "salish_subyearling"),
                    m_mat = list(m_stocks, m_length, m_south, m_salish)) %>% 
  mutate(
    # calculate synchrony for each matrix
    synch_est = map(m_mat, .f = function(y) {
      data.frame(
        year = rownames(y),
        synch = zoo::rollapplyr(y, width = window_size, function(x) community.sync(x)$obs,
                                fill = NA, by.column = FALSE)
        )
  }))

# plot synchrony trends
synch_m_trends <- synch_tbl %>% 
  select(-m_mat) %>% 
  unnest(synch_est) %>% 
  filter(!is.na(synch)) %>% 
  mutate(year = as.numeric(year)) %>% 
  ggplot(., aes(x = year, y = synch, color = group)) +
  geom_line()

pdf(here::here("figs", "synch_m.pdf"))
synch_m_trends
dev.off()