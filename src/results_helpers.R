
library(ggplot2)
library(reshape2)

plot_genealogy_theta <- function(theta_hist, anc, edge_id = 1,
                                 base_col = "black", alive_col = "red",
                                 pch = 19, cex = .25, lwd = .25)
{
  M <- dim(theta_hist)[1L];  T <- dim(theta_hist)[2L]
  
  ## (1) rebuild *all* trajectories for that coordinate
  traj_all <- matrix(0, M, T)
  for (tt in 1:(T - 1))
    traj_all[, tt] <- theta_hist[ anc[, tt + 1L], tt, edge_id ]
  traj_all[, T] <- theta_hist[, T, edge_id]
  
  ## (2) rebuild the ones that survive to the end
  traj_surv <- matrix(0, M, T)
  traj_surv[, T] <- theta_hist[, T, edge_id]
  idx <- anc[, T]
  for (tt in (T - 1L):1L) {
    traj_surv[, tt] <- theta_hist[idx, tt, edge_id]
    idx <- anc[idx, tt]
  }
  
  ## (3) plot
  matplot(t(traj_all), col = base_col, pch = pch, cex = cex,
          type = "p", lty = 1, lwd = lwd,
          xlab = "time step", ylab = bquote(theta[.(edge_id)]),
          main = paste0("Genealogy of θ[", edge_id, "]"))
  # matlines(t(traj_surv), col = alive_col, lty = 1, lwd = lwd)
}

plot_theta_histograms <- function(theta_mean, k_set = 1:3, colors = "lightblue") {
  par(mfrow = c(1, length(k_set)), mar = c(4, 4, 2, 1))
  
  for (k in k_set) {
    hist(theta_mean[, k],
         main = bquote("Histogram of " ~ theta[.(k)]),
         xlab = bquote(theta[.(k)]),
         col = colors,
         border = "white",
         breaks = 100)
  }
}


plot_theta_paths <- function(theta_mean, theta_se, theta_true = NULL, k_values = 1:3) {
  n_time <- nrow(theta_mean)
  df_all <- do.call(rbind, lapply(k_values, function(k) {
    df <- data.frame(
      t = 1:n_time,
      mu = theta_mean[, k],
      lo = theta_mean[, k] - theta_se[, k],
      hi = theta_mean[, k] + theta_se[, k],
      k = factor(k)
    )
    
    # Handle theta_true formats
    if (!is.null(theta_true)) {
      if (length(theta_true) == 1) {
        df$true <- rep(theta_true, n_time)
      } else if (is.vector(theta_true) && length(theta_true) == n_time) {
        df$true <- theta_true
      } else if (is.matrix(theta_true) && all(dim(theta_true) == dim(theta_mean))) {
        df$true <- theta_true[, k]
      } else {
        warning("theta_true must be either a scalar, a vector of length T, or a matrix with dim = [T, K]")
        df$true <- NA
      }
    }
    
    df
  }))
  
  p <- ggplot(df_all, aes(x = t, y = mu, color = k, fill = k)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    labs(title = expression("Posterior mean ±1 s.e."),
         y = expression(theta[k]), x = "time t") +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Pastel2")
  
  if (!is.null(theta_true)) {
    p <- p + geom_line(aes(y = true, group = k), color = "black", linetype = "dashed") +
      labs(title = expression("Posterior mean ±1 s.e. with true " ~ theta[k]))
  }
  
  print(p)
}

