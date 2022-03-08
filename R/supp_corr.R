## Supplementary analyses
# All data originates with PST Chinook technical committee and was provided by
# C. Parken

library(tidyverse)

dat <- readRDS(here::here("data/salmon_data/cwt_indicator_surv_clean.RDS")) %>% 
  mutate(logit_surv = qlogis(survival))


# RAW BOX PLOTS ----------------------------------------------------------------

raw_data <- dat %>% 
  mutate(
    group = factor(
      j_group3b, 
      levels = c("north_oceantype", "north_streamtype", "sog_oceantype",
                 "sog_streamtype", "puget_oceantype", "puget_streamtype",    
                 "south_oceantype", "south_streamtype"),
      labels = c("North\nSubyearling", "North\nYearling", "SoG\nSubyearling", 
                 "SoG\nYearling", "Puget\nSubyearling", "Puget\nYearling", 
                 "South\nSubyearling",  "South\nYearling")
    ),
    smolt = factor(smolt, labels = c("Subyearling", "Yearling")),
    region_juv = factor(j_group3, levels = c("north", "sog", "puget", "south"),
                        labels = c("North", "SoG", "Puget", "South"))
  ) %>% 
  filter(!stock %in% c("TST", "AKS")) 

scale_1 <- c("#8073AC", "#E08214")
names(scale_1) <- levels(raw_data$smolt)


png(here::here("figs", "supp_figs", "survival_box.png"), height = 7, 
    width = 7.5, res = 300, units = "in")
ggplot() +
  geom_boxplot(data = raw_data, aes(x = as.factor(year), y = survival,
                                    fill = smolt)) +
  facet_grid(smolt ~ region_juv,
             scales = "free_y") +
  labs(x = "Ocean Entry Year", y = "Juvenile Marine Survival Rate") +
  scale_x_discrete(
    breaks = seq(1970, 2020, by = 10)
  ) +
  scale_fill_manual(values = scale_1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
dev.off()

png(here::here("figs", "supp_figs", "age_box.png"), height = 7, 
    width = 7.5, res = 300, units = "in")
ggplot() +
  geom_boxplot(data = raw_data, aes(x = as.factor(year), y = gen_length,
                                    fill = smolt)) +
  facet_grid(smolt ~ region_juv,
             scales = "free_y") +
  labs(x = "Ocean Entry Year", y = "Mean Age-At-Maturity") +
  scale_x_discrete(
    breaks = seq(1970, 2020, by = 10)
  ) +
  scale_fill_manual(values = scale_1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
dev.off()

png(here::here("figs", "supp_figs", "size_box.png"), height = 7, 
    width = 7.5, res = 300, units = "in")
ggplot(raw_data %>% filter(!is.na(avg_weight))) +
  geom_boxplot(aes(x = as.factor(year), y = avg_weight,
                                    fill = smolt)) +
  facet_grid(smolt ~ region_juv,
             scales = "free_y") +
  labs(x = "Ocean Entry Year", y = "Release Size") +
  scale_x_discrete(
    breaks = seq(1970, 2020, by = 10)
  ) +
  scale_fill_manual(values = scale_1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
dev.off()



## COVARIANCE ------------------------------------------------------------------

# function to visualize raw estimates
plot_foo <- function(x_in = "gen_length", y_in = "survival",
                     xlab = "Mean Age-At-Maturity", ylab = "Survival") {
  dat %>% 
    filter(!is.na(.data[[x_in]]),
           !is.na(.data[[y_in]])) %>% 
    droplevels() %>% 
    mutate(stock = fct_reorder(stock, as.numeric(j_group3b))) %>% 
    ggplot(.) +
    geom_point(aes_string(x = x_in, y = y_in, fill = "j_group3b"), shape = 21) +
    scale_fill_brewer(type = "qual", palette = 2, name = "Juvenile\nGrouping") +
    facet_wrap(~stock, scales = "free") +
    labs(x = xlab, y = ylab) +
    ggsidekick::theme_sleek() 
}

pdf(here::here("figs", "supp_figs", "covs.pdf"))
plot_foo(x_in = "gen_length", y_in = "survival")
plot_foo(x_in = "avg_weight", y_in = "survival",
         xlab = "Mean Mass", ylab = "Survival")
plot_foo(x_in = "avg_weight", y_in = "gen_length",
         xlab = "Mean Mass", ylab = "Mean Age-At-Maturity")
dev.off()


# fit gams to estimate effects
age_surv_mod <- mgcv::gam(gen_length ~ s(logit_surv, m = 2, bs = "tp") + 
                            s(logit_surv, m = 1, bs = "tp", by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat)
wt_surv_mod <- mgcv::gam(logit_surv ~ s(avg_weight, m = 2, bs = "tp") + 
                            s(avg_weight, m = 1, bs = "tp", by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat#,
                          # family = mgcv::betar(link="logit")
                         )
age_wt_mod <- mgcv::gam(gen_length ~ s(avg_weight, m = 2, bs = "tp") + 
                            s(avg_weight, m = 1, bs = "tp", by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat#,
                          # family = Gamma(link="log")
                          )

# generate predictive dataframes



## SYNCHRONY -------------------------------------------------------------------

dat %>% 
  filter(!is.na(survival)) %>% 
  group_by(j_group3b, year) %>% 
  summarize(n_stocks = length(unique(stock)), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_line(aes(x = year, y = n_stocks, colour = j_group3b))
  glimpse()
