
library(ggplot2)
library(reshape2)

plot_genealogy_theta <- function(theta_hist, anc, edge_id = 1,
                                 base_col = "black", alive_col = "red",
                                 pch = 19, cex = .25, lwd = .25,
                                 ylim = NULL)  # ← added ylim input
{
  M <- dim(theta_hist)[1L]
  T <- dim(theta_hist)[2L]
  
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
  matplot(t(traj_all),
          col = base_col, pch = pch, cex = cex,
          type = "p", lty = 1, lwd = lwd,
          xlab = "time step", ylab = bquote(theta[.(edge_id)]),
          main = paste0("Genealogy of θ[", edge_id, "]"),
          ylim = ylim)  # ← passed through
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


library(ggplot2)
library(tidyr)
library(dplyr)

plot_theta_paths_df <- function(theta_mean, theta_sd = NULL, theta_true = NULL, legend_labels = NULL, plot_title = NULL) {
  n_time <- nrow(theta_mean)
  vars <- colnames(theta_mean)
  df <- data.frame(t = 1:n_time)
  
  # Melt theta_mean
  df_long <- df %>%
    bind_cols(theta_mean) %>%
    pivot_longer(-t, names_to = "var", values_to = "mu")
  
  # Add standard error if provided
  if (!is.null(theta_sd)) {
    stopifnot(nrow(theta_mean) == nrow(theta_sd))
    df_se <- df %>%
      bind_cols(theta_sd) %>%
      pivot_longer(-t, names_to = "var", values_to = "se")
    df_long <- df_long %>%
      left_join(df_se, by = c("t", "var")) %>%
      mutate(
        lo = mu - se,
        hi = mu + se
      )
  }
  
  # Handle theta_true
  if (!is.null(theta_true)) {
    if (length(theta_true) == 1) {
      df_true <- df_long %>%
        select(t, var) %>%
        mutate(true = theta_true)
    } else if (is.vector(theta_true) && length(theta_true) == n_time) {
      df_true <- df_long %>%
        select(t, var) %>%
        group_by(var) %>%
        mutate(true = theta_true) %>%
        ungroup()
    } else if (is.data.frame(theta_true) && all(dim(theta_true) == dim(theta_mean))) {
      df_true <- df %>%
        bind_cols(theta_true) %>%
        pivot_longer(-t, names_to = "var", values_to = "true")
    } else {
      warning("theta_true must be a scalar, a vector of length T, or a data frame with same shape as theta_mean")
      df_true <- NULL
    }
  } else {
    df_true <- NULL
  }
  
  # Apply custom labels
  if (!is.null(legend_labels)) {
    stopifnot(length(legend_labels) == length(vars))
    names(legend_labels) <- vars
    df_long$var <- factor(df_long$var, levels = vars, labels = legend_labels)
    if (!is.null(df_true)) {
      df_true$var <- factor(df_true$var, levels = vars, labels = legend_labels)
    }
  }
  
  # Build the plot
  p <- ggplot(df_long, aes(x = t))
  
  # Add ribbon if se is available
  if (!is.null(theta_sd)) {
    p <- p + geom_ribbon(aes(ymin = lo, ymax = hi, fill = var), alpha = 0.15, color = NA)
  }
  
  # Add posterior mean line
  p <- p + geom_line(aes(y = mu, color = var), linewidth = 0.9)
  
  # Add true value line
  if (!is.null(df_true)) {
    p <- p + geom_line(data = df_true,
                       aes(x = t, y = true, group = var, linetype = "True value"),
                       color = "black", linewidth = 0.8)
  }
  
  # Final styling
  p <- p +
    labs(
      title = plot_title %||% expression("Posterior mean ±1 s.e."),
      y = expression(theta), x = "time t",
      color = NULL, fill = NULL, linetype = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.15),
      legend.justification = c("right", "bottom"),
      legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.text = element_text(size = 13)
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Pastel2") +
    scale_linetype_manual(values = c("True value" = "dashed")) +
    guides(
      color = guide_legend(order = 1),
      fill = guide_legend(order = 1),
      linetype = guide_legend(order = 2, override.aes = list(color = "black"))
    )
  
  print(p)
  return(p)
}







