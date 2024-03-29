## function to prepare predicted fits for plot_fitted_bayes
# modelfit <- surv_dfa[[2]]; names <- surv_tbl$names[[2]];
# years <- surv_tbl$years[[2]]
# descend_order = FALSE
fitted_preds <- function(modelfit, names = NULL, years = NULL,
                         descend_order = FALSE, subset = NULL,
                         year1_last_mean = 2011) {
  n_ts <- dim(modelfit$data)[1]
  n_years <- dim(modelfit$data)[2]
  if (is.null(years)) {
    years <- seq_len(n_years)
  }
  pred <- predicted(modelfit)
  
  df_pred <- data.frame(ID = rep(seq_len(n_ts), n_years),
                        Time = sort(rep(years, n_ts)),
                        mean = c(t(apply(pred, c(3, 4), mean))),
                        lo = c(t(apply(pred, c(3, 4), quantile, 0.05))),
                        hi = c(t(apply(pred, c(3, 4), quantile, 0.95)))
                        ) %>% 
    mutate(stock = names$stock[ID])
  
  df_obs <- data.frame(ID = rep(seq_len(n_ts), n_years),
                       Time = sort(rep(years,
                                       n_ts)),
                       obs_y = c(modelfit$data)) %>%
    filter(!is.na(obs_y)) 
  
  # new categorical version
  last_gen_mean <-  final_prob(modelfit = modelfit, names = names, 
                               years = years, year1_last_mean = year1_last_mean
                               ) %>% 
    mutate(prob = ifelse(last_mean > 0, prob_above_0, prob_below_0)) %>% 
    select(-c(prob_above_0, prob_below_0))
  
  out <- df_pred %>%
    left_join(.,
              last_gen_mean, 
              by = "stock") %>%
    left_join(., df_obs, by = c("ID", "Time")) 
  
  if (!is.null(subset)) {
    samp_seq <- sample(unique(df_pred$ID), size = subset)
    out <- out %>% 
      filter(ID %in% samp_seq)
  }
  
  out %>% 
    mutate(ID = names$stock_name[ID] %>% 
             as.factor(.) %>% 
             fct_reorder(., last_mean, .desc = descend_order)) 
}

## function to plot fits in link space (based on bayesdfa::plot_fitted)
# df_pred <- surv_pred_list[[1]]
plot_fitted_pred <- function(df_pred, #ylab = NULL, 
                             print_x = TRUE, 
                             col_ramp = c(-1, 1),
                             col_ramp_direction = -1,
                             facet_row = NULL, facet_col = NULL,
                             leg_name = NULL,
                             year1_last_mean = 2011,
                             drop = TRUE) {  
  #limits for y axis
  y_lims <- max(abs(df_pred$obs_y), na.rm = T) * c(-1, 1)
  x_int <- year1_last_mean
  
  #make palette for last five year mean based on bins and col_ramp values 
  breaks <- seq(min(col_ramp), max(col_ramp), length.out = 9)
  df_pred$color_ids <- cut(df_pred$last_mean, 
                           breaks=breaks, 
                           include.lowest=TRUE, 
                           right=FALSE)
  col_pal <- c("#a50f15", "#de2d26", "#fb6a4a", "#fc9272", "#9ecae1", "#6baed6",
               "#3182bd", "#08519c", "grey60")
  names(col_pal) <- c(levels(df_pred$color_ids), "historic")
  # replace color ID label so that low probabilities are historic (i.e. grey)
  df_pred$color_ids2 <- ifelse(df_pred$prob < 0.90, 
                               "historic", 
                               as.character(df_pred$color_ids))
  
  dum <- df_pred %>% 
    group_by(stock) %>% 
    #calculate SD of ts for horizontal line
    mutate(ts_mean = mean(mean),
           ts_mean_sd = sd(mean)) %>% 
    ungroup()
  
  labs <- df_pred %>% 
    filter(!is.na(obs_y)) %>% 
    group_by(ID) %>% 
    tally()
  
  p <- ggplot(dum %>% filter(Time >= x_int),
              aes_string(x = "Time", y = "mean")) + 
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "color_ids2",
                           fill = "color_ids2"), alpha = 0.6) + 
    geom_line(aes_string(colour = "color_ids2"), 
              size = 1.25) +
    geom_ribbon(data = dum %>% filter(Time <= x_int),
                aes_string(ymin = "lo", ymax = "hi"),
                fill = "grey60", colour = "grey60", alpha = 0.6) +
    geom_line(data = dum %>% filter(Time <= x_int),
              size = 1) +
    geom_hline(aes(yintercept = ts_mean), lty = 2) +
    # geom_hline(aes(yintercept = ts_mean + ts_mean_sd), lty = 3) +
    # geom_hline(aes(yintercept = ts_mean - ts_mean_sd), lty = 3) +
    geom_vline(xintercept = x_int, lty = 1, alpha = 0.6) +
    scale_fill_manual(values = col_pal) +
    scale_colour_manual(values = col_pal) +
    scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
    geom_point(data = dum %>% filter(!is.na(obs_y)), 
               aes_string(x = "Time", y = "obs_y"),  
               size = 1, alpha = 0.6, shape = 21, fill = "black") + 
    facet_wrap(~ID, nrow = facet_row, ncol = facet_col, drop = drop) +
    ggsidekick::theme_sleek() +
    coord_cartesian(y = y_lims) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank(),
          legend.position = "none",
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank()) +
    geom_text(
      data = labs, aes(x = -Inf, y = -Inf, label = n),
      hjust = -0.2, vjust = -0.4
    )

  if (print_x == FALSE) {
    p <- p + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  return(p)
}

## function to plot fits in real space (based on bayesdfa::plot_fitted)
# df_pred <- real_surv_pred_list[[1]]
plot_fitted_pred_real <- function(df_pred, #ylab = NULL, 
                             y_lims = NULL, 
                             print_x = TRUE, 
                             facet_row = NULL, facet_col = NULL,
                             year1_last_mean = 2011
                             ) {  
  x_int <- year1_last_mean
  y_int <- df_pred %>% 
    group_by(ID) %>% 
    summarize(ts_mean_logit = mean(uncent_mean_logit),
              sd_mean_logit = sd(uncent_mean_logit), 
              ts_mean_sd_lo = plogis(ts_mean_logit + 
                                       (qnorm(0.025) * sd_mean_logit)),
              ts_mean_sd_hi = plogis(ts_mean_logit + 
                                       (qnorm(0.975) * sd_mean_logit)),
              ts_uncent_mean = mean(uncent_mean), 
              .groups = "drop") %>% 
    distinct()
  
  # specify that color greyed out if relatively uncertain (in logit space)
  df_pred2 <- df_pred %>% 
    left_join(., y_int, by = c("ID")) %>% 
    mutate(
      color_id = case_when(
        prob < 0.9 ~ "historic",
        last_mean < ts_mean_sd_lo ~ "very low",
        ts_mean_sd_lo < last_mean & last_mean < ts_uncent_mean ~ "low",
        ts_mean_sd_hi > last_mean & last_mean  > ts_uncent_mean ~ "high",
        last_mean > ts_mean_sd_hi ~ "very high"
      ),
      color_id = fct_reorder(as.factor(color_id), 
                             last_mean - ts_uncent_mean),
      # necessary to order correctly
      ID_key = fct_reorder(as.factor(ID), as.numeric(color_id))
      ) %>% 
    droplevels()
  y_int2 <- y_int %>% 
    left_join(., df_pred2 %>% select(ID, ID_key) %>% distinct(), by = "ID")
  
  #make palette for last five year mean based on bins and col_ramp values 
  col_pal <- c("#a50f15",  "#fc9272", "#9ecae1",  "#08519c", "grey60")
  names(col_pal) <- c("very low", "low", "high", "very high", "historic")
   
  labs <- df_pred2 %>% 
    filter(!is.na(obs_y)) %>% 
    group_by(ID_key) %>% 
    tally()
  
  p <- ggplot(df_pred2 %>% filter(Time >= x_int),
              aes_string(x = "Time", y = "uncent_mean")) + 
    geom_ribbon(aes_string(ymin = "uncent_lo", ymax = "uncent_hi", 
                           colour = "color_id",
                           fill = "color_id"), alpha = 0.6) + 
    geom_line(aes_string(colour = "color_id"), 
              size = 1.25) +
    geom_ribbon(data = df_pred2 %>% filter(Time <= x_int),
                aes_string(ymin = "uncent_lo", ymax = "uncent_hi"),
                fill = "grey60", colour = "grey60", alpha = 0.6) +
    geom_line(data = df_pred2 %>% filter(Time <= x_int), 
              size = 1) +
    geom_hline(data = y_int2, aes(yintercept = ts_uncent_mean), lty = 2) +
    # geom_hline(data = y_int2, aes(yintercept = ts_mean_sd_hi), lty = 3) +
    # geom_hline(data = y_int2, aes(yintercept = ts_mean_sd_lo), lty = 3) +
    geom_vline(xintercept = x_int, lty = 1, alpha = 0.6) +
    scale_fill_manual(values = col_pal) +
    scale_colour_manual(values = col_pal) +
    geom_point(data = df_pred2 %>% filter(!is.na(obs_y)), 
               aes_string(x = "Time", y = "survival"),  
               size = 1, alpha = 0.6, shape = 21, fill = "black") + 
    facet_wrap(~ID_key,  
               nrow = facet_row, ncol = facet_col) +
    ggsidekick::theme_sleek() +
    coord_cartesian(y = c(0, 0.2), expand = 0) +
    scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank(),
          legend.position = "none",
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank()) +
    geom_text(
      data = labs, aes(x = -Inf, y = Inf, label = n),
      hjust = -0.2, vjust = 1.1
    )

  if (print_x == FALSE) {
    p <- p + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  return(p)
}


# as above but for uncentered data
plot_fitted_pred_uncent <- function(df_pred, #ylab = NULL, 
                                    print_x = TRUE, 
                                    col_ramp = c(-1, 1),
                                    col_ramp_direction = -1,
                                    facet_row = NULL, facet_col = NULL,
                                    leg_name = NULL,
                                    year1_last_mean = 2011,
                                    drop = TRUE) {  
  #limits for y axis
  y_lims <- c(min(abs(df_pred$obs_y), na.rm = T), 
              max(abs(df_pred$obs_y), na.rm = T))
  x_int <- year1_last_mean
  
  #make palette for last five year mean based on bins and col_ramp values 
  breaks <- seq(min(col_ramp), max(col_ramp), length.out = 9)
  df_pred$color_ids <- cut(df_pred$last_mean, 
                           breaks=breaks, 
                           include.lowest=TRUE, 
                           right=FALSE)
  col_pal <- c("#a50f15", "#de2d26", "#fb6a4a", "#fc9272", "#9ecae1", "#6baed6",
               "#3182bd", "#08519c", "grey60")
  names(col_pal) <- c(levels(df_pred$color_ids), "historic")
  # replace color ID label so that low probabilities are historic (i.e. grey)
  df_pred$color_ids2 <- ifelse(df_pred$prob < 0.90, 
                               "historic", 
                               as.character(df_pred$color_ids))
  
  dum <- df_pred %>% 
    group_by(stock) %>% 
    #calculate SD of ts for horizontal line
    mutate(ts_mean_sd = sd(mean)) %>% 
    ungroup()
  
  labs <- dum %>% 
    filter(!is.na(obs_y)) %>% 
    group_by(ID) %>% 
    tally()
  
  p <- ggplot(dum %>% filter(Time >= x_int),
              aes_string(x = "Time", y = "mean")) + 
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "color_ids2",
                           fill = "color_ids2"), alpha = 0.6) + 
    geom_line(aes_string(colour = "color_ids2"), 
              size = 1.25) +
    geom_ribbon(data = dum %>% filter(Time <= x_int),
                aes_string(ymin = "lo", ymax = "hi"),
                fill = "grey60", colour = "grey60", alpha = 0.6) +
    geom_line(data = dum %>% filter(Time <= x_int),
              size = 1) +
    geom_hline(aes(yintercept = obs_mean_age), lty = 2) +
    # geom_hline(aes(yintercept = obs_mean_age + ts_mean_sd), lty = 3) +
    # geom_hline(aes(yintercept = obs_mean_age - ts_mean_sd), lty = 3) +
    geom_vline(xintercept = x_int, lty = 1, alpha = 0.6) +
    scale_fill_manual(values = col_pal) +
    scale_colour_manual(values = col_pal) +
    scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
    geom_point(data = dum %>% filter(!is.na(obs_y)), 
               aes_string(x = "Time", y = "obs_y"),  
               size = 1, alpha = 0.6, shape = 21, fill = "black") + 
    facet_wrap(~ID, nrow = facet_row, ncol = facet_col, drop = drop) +
    ggsidekick::theme_sleek() +
    coord_cartesian(y = y_lims) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank(),
          legend.position = "none",
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank()) +
    geom_text(
      data = labs, aes(x = -Inf, y = -Inf, label = n),
      hjust = -0.2, vjust = -0.4
    )
  
  if (print_x == FALSE) {
    p <- p + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  return(p)
}


## function to calculate probability that estimates below average in last 
# n_years
# modelfit = surv_dfa[[2]]; names = surv_tbl$names[[2]]; years = surv_tbl$years[[2]]

final_prob <- function(modelfit, names, 
                       years = years, year1_last_mean = 2010,
                       year2_last_mean = NULL
                       ) {
  tt <- reshape2::melt(predicted(modelfit), 
                       varnames = c("iter", "chain", "time", "stock")) %>% 
    left_join(., 
              data.frame(year = years,
                         time = unique(.$time)),
              by = "time") 
  tt$stock <- as.factor(names$stock[tt$stock])
  
  # yr_range <- seq(max(tt$year) - (n_years - 1), max(tt$year), by = 1)
  if (is.null(year2_last_mean)) {
    year2_last_mean <- max(tt$year)
  }
  yr_range <- seq(year1_last_mean, year2_last_mean, by = 1)
  tt %>% 
    group_by(stock) %>% 
    filter(!year > year2_last_mean) %>%
    mutate(overall_mean = mean(value)) %>% 
    filter(year %in% yr_range) %>%
    group_by(stock, iter) %>%
    mutate(mean_value = mean(value)) %>% 
    group_by(stock) %>% 
      summarize( 
        last_mean = mean(mean_value),
        prob_below_0 = sum(mean_value < overall_mean) / length(mean_value),
        prob_above_0 = sum(mean_value > overall_mean) / length(mean_value)
      ) 
}


## function to prepare rotated model fit for plotting trends (based on 
# bayesdfa::plot_trends)
prep_trends <- function (rotated_modelfit, years, group) { 
  rotated <- rotated_modelfit
  n_ts <- dim(rotated$Z_rot)[2]
  n_trends <- dim(rotated$Z_rot)[3]
  n_years <- dim(rotated$trends_mean)[2]
  
  data.frame(x = c(t(rotated$trends_mean)), 
             lo = c(t(rotated$trends_lower)), 
             hi = c(t(rotated$trends_upper)), 
             trend = paste0("Trend ",
                            sort(rep(seq_len(n_trends), n_years))), 
             time = rep(years, n_trends),
             group = group)
}


## function to plot trends (based on bayesdfa::plot_trends)
plot_one_trend <- function(trend_dat, facet_var = FALSE) {
  p <- ggplot(trend_dat, 
         aes_string(x = "time", y = "x")) + 
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "life_history",
                           fill = "life_history"), 
                alpha = 0.4) + 
    geom_line(aes_string(colour = "life_history"), size = 1.2) + 
    # scale_colour_brewer(type = "qual", name = "") +
    # scale_fill_brewer(type = "qual", name = "") +
    geom_hline(yintercept = 0, lty = 2) +
    # xlab("Brood Year") + 
    ylab("Estimated Trend") +
    scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
    facet_wrap(~group, nrow = 1) +
    ggsidekick::theme_sleek() + 
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.title.x = element_blank())
  
  if (facet_var == TRUE) {
    p <- p + 
      facet_grid(group~var)
  }
  
  return(p)
}


## function to prep regime model fit for plotting (based on 
# bayesdfa::plot_regime_model)
prep_regime <- function(regime_model, probs = c(0.05, 0.95),
                        regime_prob_threshold = 0.9, flip_regimes = FALSE,
                        years, group) {
  gamma_tk <- rstan::extract(regime_model$model, pars = "gamma_tk")[[1]]
  mu_k <- rstan::extract(regime_model$model, pars = "mu_k")[[1]]
  l <- apply(gamma_tk, 2:3, quantile, probs = probs[[1]])
  u <- apply(gamma_tk, 2:3, quantile, probs = probs[[2]])
  med <- apply(gamma_tk, 2:3, quantile, probs = 0.5)
  range01 <- function(x) (x - min(x))/(max(x) - min(x))
  mu_k_low <- apply(mu_k, 2, quantile, probs = probs[[1]])
  mu_k_high <- apply(mu_k, 2, quantile, probs = probs[[2]])
  mu_k <- apply(mu_k, 2, median)
  confident_regimes <- apply(
    gamma_tk, 2:3, function(x) mean(x > 0.5) > regime_prob_threshold
  )
  regime_indexes <- apply(confident_regimes, 1, function(x) {
    w <- which(x)
    if (length(w) == 0) 
      NA
    else w
  })
  
  #should regimes be flipped for plotting
  if (flip_regimes) {
    mu_k <- 1 - mu_k
    u <- 1 - u
    l <- 1 - l
    med <- 1 - med
  }
  
  plot_prob_indices <- seq_len(ncol(med))
  df_l <- reshape2::melt(l, varnames = c("Time", "State"), 
                         value.name = "lwr")
  df_u <- reshape2::melt(u, varnames = c("Time", "State"), 
                         value.name = "upr")
  df_m <- reshape2::melt(med, varnames = c("Time", "State"), 
                         value.name = "median")
  df_y <- data.frame(y = range01(regime_model$y), 
                     Time = seq_along(regime_model$y))
  dplyr::inner_join(df_l, df_u, by = c("Time", "State")) %>% 
    dplyr::inner_join(df_m, by = c("Time", "State")) %>% 
    dplyr::filter(.data$State %in% plot_prob_indices) %>% 
    dplyr::mutate(State = paste("State", .data$State),
                  time = rep(years, length(unique(State))),
                  group = group)
}


## function to plot regimes (based on bayesdfa::plot_trends/plot_regime_model)
plot_one_regime <- function(regime_dat, facet_var = FALSE, y_lab = NULL) {
  p <- ggplot(regime_dat, 
              aes_string(x = "time", y = "median")) + 
    geom_ribbon(aes_string(ymin = "lwr", ymax = "upr", colour = "life_history",
                           fill = "life_history"), 
                alpha = 0.4, lty = 6) + 
    geom_line(aes_string(colour = "life_history"), size = 1.2, lty = 6) + 
    # scale_colour_brewer(type = "qual", name = "") +
    # scale_fill_brewer(type = "qual", name = "") +
    # xlab("Brood Year") + 
    ylab(y_lab) +
    scale_x_continuous(limits = c(1972, 2018), expand = c(0, 0)) +
    facet_wrap(~group, nrow = 1) +
    ggsidekick::theme_sleek() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.title.x = element_blank()
      )
  
  if (facet_var == TRUE) {
    p <- p +
      facet_grid(group ~ var)
  }
  
  return(p)
}



## function to prepare rotated model fit for plotting loadings (based on 
# bayesdfa::plot_loadings)
prep_loadings <- function (rotated_modelfit, names, group, conf_level = 0.95) { 
  v <- reshape2::melt(rotated_modelfit$Z_rot, 
                      varnames = c("iter", "name", "trend")) 
  v$name <- as.factor(names$stock[v$name])
  v %>% 
    mutate(trend = as.factor(paste0("Trend ", trend))) %>% 
    group_by(name, trend) %>% 
    mutate(q_lower = sum(value < 0) / length(value), 
           q_upper = 1 - q_lower, 
           prob_diff0 = max(q_lower, q_upper),
           group = group) 
}


##function to plot loadings
plot_load <- function(x, group = NULL, guides = FALSE, y_lims = c(-0.5, 0.5)) {
  p <- ggplot(x, aes_string(x = "name", y = "value", fill = "trend", 
                            alpha = "prob_diff0")) + 
    scale_alpha_continuous(name = "Probability\nDifferent") +
    scale_fill_brewer(name = "", palette = "Paired") +
    geom_violin(color = NA, position = position_dodge(0.3)) + 
    geom_hline(yintercept = 0, lty = 2) + 
    coord_flip() + 
    xlab("Time Series") + 
    ylab("Loading") +
    scale_y_continuous(limits = y_lims, expand = c(0, 0)) +
    ggsidekick::theme_sleek() +
    guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
    theme(#axis.text.y = element_text(angle = 45, vjust = -1, size = 7),
          axis.title = element_blank()) +
    annotate("text", x = Inf, y = -Inf, label = group, hjust = -0.05, 
             vjust = 1.1, size = 3.5)
  
  if (guides == FALSE) {
    p <- p +
      theme(legend.position = "none")
  }
  
  return(p)
}

