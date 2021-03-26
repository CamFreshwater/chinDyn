## function to prepare predicted fits for plot_fitted_bayes
fitted_preds <- function(modelfit, names = NULL, years = NULL,
                         descend_order = FALSE, subset = NULL) {
  n_ts <- dim(modelfit$data)[1]
  n_years <- dim(modelfit$data)[2]
  if (is.null(years)) {
    years <- seq_len(n_years)
  }
  pred <- predicted(modelfit)
  df_pred <- data.frame(ID = rep(seq_len(n_ts), n_years),
                        Time = sort(rep(years,
                                        n_ts)),
                        mean = c(t(apply(pred, c(3, 4), mean))),
                        lo = c(t(apply(pred, c(3, 4), quantile, 0.025))),
                        hi = c(t(apply(pred, c(3, 4), quantile, 0.975)))
                        )
  df_obs <- data.frame(ID = rep(seq_len(n_ts), n_years),
                       Time = sort(rep(years,
                                       n_ts)),
                       obs_y = c(modelfit$data)) %>%
    filter(!is.na(obs_y))
  
  # calculate mean over last generation to color code
  gen_seq <- seq(max(df_pred$Time - 5), max(df_pred$Time), by = 1)
  last_gen_mean <- df_pred %>%
    filter(Time %in% gen_seq) %>%
    group_by(ID) %>% 
    summarize(last_mean = mean(mean))
  
  out <- df_pred %>%
    left_join(.,
              last_gen_mean, 
              by = "ID") %>%
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


## function to plot fits (based on bayesdfa::plot_fitted)
plot_fitted_pred <- function(df_pred, #ylab = NULL, 
                             print_x = TRUE, 
                             col_ramp = c(-1, 1), 
                             col_ramp_direction = -1, facet_col = NULL,
                             leg_name = NULL) {  
  #limits for y axis
  y_lims <- max(df_pred$obs_y, na.rm = T) * c(-1, 1)
  x_int <- max(df_pred$Time, na.rm = T) - 5
  
  p <- ggplot(df_pred, aes_string(x = "Time", y = "mean")) + 
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "last_mean",
                           fill = "last_mean"), 
                alpha = 0.6) + 
    geom_line(aes_string(colour = "last_mean"), size = 1.25) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_vline(xintercept = x_int, lty = 1, alpha = 0.6) +
    scale_fill_distiller(type = "div", limit = col_ramp, 
                         direction = col_ramp_direction, 
                         palette = "RdYlBu", name = leg_name) +
    scale_colour_distiller(type = "div", limit = col_ramp, 
                           direction = col_ramp_direction, 
                           palette = "RdYlBu", name = leg_name) +
    # scale_y_continuous(name = ylab, position = 'right', sec.axis = dup_axis(),
    #                    labels = scales::number_format(accuracy = 1)) + 
    scale_x_continuous(limits = c(1972, 2016), expand = c(0, 0)) +
    geom_point(aes_string(x = "Time", y = "obs_y"),  
               size = 1, alpha = 0.6, shape = 21, fill = "black") + 
    # facet_wrap(~ID, ncol = facet_col) +
    facet_wrap(~ID, nrow = 1) +
    ggsidekick::theme_sleek() +
    coord_cartesian(y = y_lims) +
    # labs(y = ylab) +
    theme(axis.title.x = element_blank(),
          axis.title.y.left = element_blank(),
          legend.position = "none",
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())
  
  if (print_x == FALSE) {
    p <- p + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  return(p)
}



## function to calculate probability that estimates below average in last 
# n_years
final_prob <- function (modelfit, names, n_years = 5) {
  tt <- reshape2::melt(predicted(modelfit), 
                       varnames = c("iter", "chain", "year", "stock")) 
  tt$stock <- as.factor(names$stock_name[tt$stock])
  yr_range <- seq(max(tt$year) - (n_years - 1), max(tt$year), by = 1)
  tt %>% 
    filter(year %in% yr_range) %>% 
    group_by(stock) %>% 
    summarize( 
      prob_below_0 = sum(value < 0) / length(value),
      prob_above_0 = sum(value > 0) / length(value)
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
plot_one_trend <- function(trend_dat) {
  ggplot(trend_dat, 
         aes_string(x = "time", y = "x")) + 
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi", colour = "life_history",
                           fill = "life_history"), 
                alpha = 0.4) + 
    geom_line(aes_string(colour = "life_history"), size = 1.2) + 
    scale_colour_brewer(type = "qual", name = "") +
    scale_fill_brewer(type = "qual", name = "") +
    geom_hline(yintercept = 0, lty = 2) +
    xlab("Brood Year") + 
    ylab("Estimated Trend") +
    scale_x_continuous(limits = c(1972, 2016), expand = c(0, 0)) +
    facet_grid(group ~ var) +
    ggsidekick::theme_sleek() + 
    theme(legend.position = "top")
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
    scale_fill_brewer(name = "", palette = "Set2") +
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

