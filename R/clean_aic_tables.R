## Maximum Likelihood Model Selection Tables 
# June 28, 2021
# Clean AIC tables generated in surv_dfa.R and gen_length_dfa.R, export as .csv

library(tidyverse)

surv_tab <- readRDS(here::here("data", "survival_fits", "marss_aic_tab.RDS"))
gen_tab <- readRDS(here::here("data", "generation_fits",
                              "marss_aic_tab_scalingA_raw.RDS"))

foo <- function(x) {
  x %>% 
    mutate(
     H2 = case_when(
        grepl("j_", H) ~ "juvenile distribution",
        grepl("a_", H) ~ "adult distribution",
        TRUE ~ NA_character_
      ),
      H3 = case_when(
        grepl("group1", H) ~ "very coarse",
        grepl("group2", H) ~ "coarse",
        grepl("group3", H) ~ "fine",
        grepl("group4", H) ~ "very fine",
        TRUE ~ NA_character_
      ),
      H_new = case_when(
        grepl("smolt", H) ~ "life history",
        grepl("b", H) ~ paste(H2, H3, "+ life history", sep = " "),
        grepl("group", H) ~ paste(H2, H3, sep = " "),
        TRUE ~ as.character(H)
      )
      ) %>% 
    select("Trait Model" = H_new, "Number of Trends" = m, 
           "n Parameters" = num.param,  
           "Log-Likelihood" = logLik, AICc, deltaAICc, 
           "AIC Weight" = aic_weight)
}

write.csv(foo(surv_tab), 
          here::here("data", "manuscript_tables", "cleaned_marss_aic_surv.csv"),
          row.names = FALSE)

write.csv(foo(gen_tab), 
          here::here("data", "manuscript_tables", "cleaned_marss_aic_age.csv"),
          row.names = FALSE)
