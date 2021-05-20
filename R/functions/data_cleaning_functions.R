## Helper functions for data cleaning

#Function to spread and label input matrices for bayesdfa
make_mat <- function(x) {
  mat1 <- x %>%
    select(year, stock, gen_length) %>%
    spread(key = stock, value = gen_length) %>%
    as.matrix() 
  out_mat <- t(mat1[, 2:ncol(mat1)])
  colnames(out_mat) <- mat1[, "year"]
  return(out_mat)
}

# Function to identify regimes
regime_f <- function(rots_in) {
  dum <- vector(nrow(rots_in$trends_mean), mode = "list")
  for(i in 1:nrow(rots_in$trends_mean)) {
    dum[[i]] <- find_regimes(
      rots_in$trends_mean[i, ], 
      sds = (rots_in$trends_upper - rots_in$trends_mean)[i, ] / 1.96,
      max_regimes = 3,
      iter = 3000,
      control = list(adapt_delta = 0.99, max_treedepth = 20)
    )
  }
  return(dum)
}

# Function to make covariate dataframe to pass to fit_dfa
make_cov_df <- function(raw_dat, group, mat, 
                        scale = c(NULL, "scale", "center")) {
  wt <- raw_dat %>% 
    filter(j_group3b == group) %>% 
    select(brood_year, stock, avg_weight)
  
  cov_dat <- expand.grid(stock = dimnames(mat)[[1]],
              brood_year = as.numeric(dimnames(mat)[[2]]),
              covariate = 1) %>% 
    mutate(time = as.numeric(as.factor(brood_year)),
           timeseries = as.numeric(stock)) %>% 
    left_join(., wt, by = c("stock", "brood_year")) %>% 
    dplyr::select(time, timeseries, covariate, value = avg_weight) %>% 
    filter(!is.na(value))
  
  if (scale == "scale") {
    cov_dat <- cov_dat %>% 
      group_by(timeseries) %>% 
      mutate(value = as.numeric(scale(value, center = TRUE, scale = TRUE))) %>% 
      ungroup()
  }
  if (scale == "center") {
    cov_dat <- cov_dat %>% 
      group_by(timeseries) %>% 
      mutate(value = as.numeric(scale(value, center = TRUE, scale = FALSE))) %>% 
      ungroup()
  }
  
  return(as.data.frame(cov_dat))
}