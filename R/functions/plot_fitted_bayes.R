## corrected version of bayesdfa::plot_fitted

plot_fitted_bayes <- function (modelfit, names = NULL, years = NULL) 
{
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
                       y = c(modelfit$data))
  if (!is.null(names)) {
    df_pred$ID <- names[df_pred$ID]
    df_obs$ID <- names[df_obs$ID]
  }
  p1 <- ggplot(df_pred, aes_string(x = "Time", y = "mean")) + 
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), 
                alpha = 0.4) + geom_line() + 
    geom_point(data = df_obs, 
               aes_string(x = "Time", y = "y"), col = "red", 
               size = 0.5, alpha = 0.4) + 
    facet_wrap("ID", scales = "free_y") + 
    xlab("Time") + ylab("") +
    ggsidekick::theme_sleek()
  p1
}