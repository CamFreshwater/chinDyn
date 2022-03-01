## Supplementary analyses
# All data originates with PST Chinook technical committee and was provided by
# C. Parken

library(tidyverse)

dat <- readRDS(here::here("data/salmon_data/cwt_indicator_surv_clean.RDS")) %>% 
  mutate(logit_surv = qlogis(survival))


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
age_surv_mod <- mgcv::gam(survival ~ s(gen_length, m = 2, bs = "tp") + 
                            s(gen_length, m = 1, bs = "tp", by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat,
                          family = mgcv::betar(link="logit"))
wt_surv_mod <- mgcv::gam(survival ~ s(avg_weight, m = 2, bs = "tp") + 
                            s(avg_weight, m = 1, bs = "tp", by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat,
                          family = mgcv::betar(link="logit"))
age_surv_mod <- mgcv::gam(gen_length ~ s(avg_weight, m = 2, bs = "tp") + 
                            s(avg_weight, m = 1, bs = "tp", by = j_group3b) +
                            s(stock, bs = "re"),
                          data = dat,
                          family = Gamma(link="log"))


## SYNCHRONY -------------------------------------------------------------------

dat %>% 
  filter(!is.na(survival)) %>% 
  group_by(j_group3b, year) %>% 
  summarize(n_stocks = length(unique(stock)), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_line(aes(x = year, y = n_stocks, colour = j_group3b))
  glimpse()
