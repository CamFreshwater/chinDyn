## Fit MARSS then Bayes DFA models to mean generation time data
# Jan 8 2021
# Update on surv_bayesDFA; instead of deciding groupings a priori, first perform
# model selection for groupings with MARSS (faster to converge than bayesDFA), 
# then fit group specific models with bayesdfa

library(MARSS)
library(tidyverse)

gen_raw <- readRDS(here::here("data", "salmonData", "cwt_indicator_surv_clean.RDS")) 

#remove stocks with no gen_length data
gen <- gen_raw %>% 
  filter(!is.na(gen_length)) %>% 
  mutate(gen_z = as.numeric(scale(gen_length)))
  
# dataframe of only stocks and adult groupings
stk_tbl <- gen %>% 
  select(stock, stock_name, run, a_group:a_group3) %>% 
  distinct()


## EXPLORATORY -----------------------------------------------------------------

#plot raw generation length data
gen  %>% 
  ggplot(.) +
  geom_point(aes(x = year, y = gen_length, fill = a_group2), shape = 21) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group2))) +
  theme(legend.position = "top") +
  labs(y = "Mean Generation Length") +
  ggsidekick::theme_sleek()


#distribution of generation length 
gen %>%
  # filter(!is.na(gen_length)) %>% 
  group_by(stock) %>% 
  ggplot(.) +
  geom_histogram(aes(x = gen_z, fill = a_group)) +
  facet_wrap(~ fct_reorder(stock, as.numeric(a_group))) +
  theme(legend.position = "top") +
  ggsidekick::theme_sleek()


## MARSS MODEL RUNS ------------------------------------------------------------

# make matrix of natural mortality rates
gen_mat1 <- gen %>% 
  select(year, stock, gen_length) %>% 
  pivot_wider(names_from = stock, values_from = gen_length) %>% 
  arrange(year) %>% 
  as.matrix() %>% 
  t()
gen_mat <- gen_mat1[2:nrow(gen_mat1), ]
colnames(gen_mat) <- seq(min(gen$year), max(gen$year), by = 1)

n_ts <- nrow(gen_mat)
tt <- ncol(gen_mat)


## Generic MARSS approach

# specify the z models based on different groups
z1 <- factor(stk_tbl$run)
z2 <-  factor(stk_tbl$a_group)
z3 <- factor(stk_tbl$a_group2)
z4 <- factor(stk_tbl$a_group3)
z_models <- list(z1, z2, z3, z4)
names(z_models) <- c("run", "off-shelf", "region", "region2")

q_models <- c("diagonal and equal", 
              "diagonal and unequal", 
              "equalvarcov",
              "unconstrained")

U <- "unequal"
R <- "diagonal and equal"
A <- "scaling"
B <- "identity"
x0 <- "unequal"
V0 <- "zero"
model_constants <- list(U = U, R = R, A = A, B = B, x0 = x0, V0 = V0)


# function to fit models
fit_marss <- function(z_name, z_in, q_in) {
  fit_model <- c(list(Z = z_in, Q = q_in), model_constants)
  fit <- MARSS(gen_mat, model = fit_model,
               silent = FALSE, control = list(minit = 100, maxit = 500))
  #use BFGS except when equalvarcov (can't fit)
  if (fit$convergence != 0 & q_in != "equalvarcov"){
    fit <- MARSS(gen_mat, 
                 model = fit_model, 
                 control = list(maxit=4000, trace=1),
                 inits = as.matrix(coef(fit)[1]),
                 method = "BFGS")
  }
  out <- data.frame(
    H = z_name, Q = q_in, U = U,
    logLik = fit$logLik, AICc = fit$AICc, num.param = fit$num.params,
    m = length(unique(z_in)),
    num.iter = fit$numIter, 
    converged = !fit$convergence,
    stringsAsFactors = FALSE)
  list(fit = fit, out = out)
}

# tibble containing model combinations
mod_names = expand.grid(q = q_models, z = names(z_models)) %>% 
  mutate(name = paste(z, q, sep = "-"))
mod_tbl <- tibble(
  mod_name = mod_names$name,
  z_name = mod_names$z,
  q_name = mod_names$q,
  z_models = rep(z_models, each = 4),
  q_models = rep(q_models, times = 4)
)

# MARSS(gen_mat, model =  c(list(Z = z_models[[1]], 
#                                Q = "unconstrained"), 
#                           model_constants),
#       silent = FALSE, control = list(maxit = 100),
#       method = "BFGS")

# fit generic MARSS models
marss_list <- furrr::future_pmap(list(z_name = mod_tbl$z_name,
                                      z_in = mod_tbl$z_models,
                                      q_in = mod_tbl$q_models),
                                 .f = fit_marss,
                                 .progress = TRUE)

marss_aic_tab <- purrr::map(marss_list, "out") %>% 
  bind_rows() %>% 
  arrange(AICc) %>% 
  mutate(deltaAICc = AICc - min(AICc),
         rel_like = exp(-1 * deltaAICc / 2),
         aic_weight = rel_like / sum(rel_like))

saveRDS(marss_list, here::here("data", "generation_fits",
                               "marss_selection_fits.RDS"))
saveRDS(marss_aic_tab, here::here("data", "generation_fits",
                                  "marss_aic_tab.RDS"))
