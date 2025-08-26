
library(here)
source(here::here("src","config.R"))


######################################################################

res_stand <- readRDS("simul_results/static_dgp/standard_3_test.rds")
res_block <- readRDS("simul_results/static_dgp/block_stat_3.rds")

res_stand$log_model_evidence
res_block$log_model_evidence

res_stand <- readRDS("simul_results/static_dgp/standard_3_pi05_nmh1_N250_KM")
res_SS <- readRDS("simul_results/static_dgp/standard_3_pi05_nmh1_N250_SSVS")

res_stand$theta_mean


plot_theta_paths(res_stand$theta_mean, res_stand$theta_se, k=3, theta_true = 0.2)
plot_theta_paths(res_SS$theta_mean, res_SS$theta_se, k=3, theta_true = 0.2)




k=1
df_mean <- data.frame(
  block = res_block$theta_mean[, k],
  stand = res_stand$theta_mean[, k]
)

df_se <- data.frame(
  block = res_block$theta_se[, k],
  stand = res_stand$theta_se[, k]
)
true_values_list = list(0.66,  0.70,  0.53,  0.38,  0.70,  0.70,  0.58)
p <- plot_theta_paths_df(df_mean, df_se, theta_true=0.7, legend_labels = c("Tempering", "Standard"), 
                         plot_title = expression("Posterior path of " * theta[12]))

library(patchwork)

true_values_list <- rev(list(0.03,  0.15,  0.20,  0.00,  0.00))
true_values_list <- rev(list(0.2,  0.15,  0.03,  0.00,  0.00))
true_values_list <- rev(list(0,0,0,0))


p <- 18
plots <- lapply(seq_along(true_values_list), function(k) {
  df_mean <- data.frame(
    block = res_block$gamma_mean[, k+p],
    stand = res_stand$gamma_mean[, k+p]
  )
  
  df_se <- data.frame(
    block = res_block$gamma_se[, k+p],
    stand = res_stand$gamma_se[, k+p]
  )
  
  plot_theta_paths_df(
    theta_mean = df_mean,
    theta_sd   = df_se,
    theta_true = true_values_list[[k]],
    legend_labels = c("Tempering", "Standard"),
    plot_title = bquote("Posterior path of " * gamma[.(k)])
  )
})

# Add empty placeholders to fill 3x3 layout
while (length(plots) < 4) {
  plots[[length(plots) + 1]] <- patchwork::plot_spacer()
}

# Combine and display in a 3x3 grid
combined_plot <- wrap_plots(plots, ncol = 2)
combined_plot


ggsave(filename = "figures/static_8d_tr4_gamma.png", plot = combined_plot, width = 12, height = 12) 


df_ESS <- data.frame(
  block = res_block$diag_log$ESS,
  #stand = res_stand$diag_log$ESS,
  threshold = rep(1000, length( res_stand$diag_log$ESS))
)

df_ESS <- df_ESS[1:3000,]

sum(df_ESS$block<1000)

p <- plot_theta_paths_df(df_ESS,  plot_title = "Effective Sample Size - 169 hits")
ggsave(filename = "figures/static_8d_ESS_block.png", plot = p, width = 8, height = 6)


df_unique <- data.frame(
  block = res_block$diag_log$unique,
  #stand = res_stand$diag_log$unique,
  initial_value = rep(2000, length( res_stand$diag_log$unique))
)

df_unique <- df_unique[1:3000,]

p <- plot_theta_paths_df(df_unique,  plot_title = "Unique number of particles")
ggsave(filename = "figures/static_3d_unique_block.png", plot = p, width = 8, height = 6)







#U  <- sim_static_cop_3(N = 1000)
# U  <- sim_ar1_copula_corr_3(N = 1000)
# N  <- nrow(U)
# d  <- ncol(U)
# cfg <- build_cfg(d)


#results <- run_standard_smc(U, cfg, type="standard")
# results <- run_block_smc(U, cfg, type="block")
# 
# 
# 
# 
# cat("\n\n===== FINAL MODEL EVALUATION =====\n")
# cat(sprintf("Log Model Evidence: %.4f\n", results$log_model_evidence))

#================================================================================
# Save results
#================================================================================

# results[["cfg"]] <- cfg
# saveRDS(results, file = "simul_results/block_dynamic_3.rds")



#================================================================================
# Explore results
#================================================================================

plot(res_stand$diag_log$ESS, type="l")
#  
plot_genealogy_theta(res_SS$theta_hist, res_SS$ancestorIndices, edge_id = 1)   
#  
# plot_theta_histograms(results$theta_mean, k_set = c(1))
# 
# plot_theta_paths(results$theta_mean, results$theta_se, k=1, theta_true = 0.6)
#  
#  


#results <- readRDS("simul_results/block_stat_3.rds")


sum(res_SS$diag_log$ESS<500)



###----------------------------------------------------------------------------------------------------------------------------------------------------------
# NEW



res <- readRDS("simul_results/static_dgp/standard_3_test")


## fam_hist[ particle , time , edge ]
fh <- res$fam_hist            # 1000 × 300 × K   ➊
T  <- dim(fh)[2]              # 300
K  <- dim(fh)[3]              # number of edges (3 here)

## --- utility: returns a data.frame for ONE edge -------------------
edge_family_ts <- function(e, fh) {
  prop <- function(code) colMeans(fh[ , , e] == code)   # over particles
  data.frame(
    t        = seq_len(dim(fh)[2]),
    indep    = prop(0L),
    gaussian = prop(1L),
    bb1      = prop(2L)
  )
}

## --- compute all edges --------------------------------------------
edge_series <- lapply(seq_len(K), edge_family_ts, fh = fh)

## --- plot each edge ------------------------------------------------
old_par <- par(mfrow = c(K, 1), mar = c(3, 4, 2, 1))     # one panel per edge
on.exit(par(old_par), add = TRUE)

cols <- c("black", "blue", "red")       # indep, gaussian, bb1
for (e in seq_len(K)) {
  df <- edge_series[[e]]
  matplot(df$t, df[ , 2:4], type = "l", lty = 1, col = cols,
          xlab = "time t", ylab = "proportion",
          main = paste("Edge", e, ": family probabilities over time"))
  legend("topright", legend = c("indep", "gaussian", "bb1"),
         col = cols, lty = 1, bty = "n")
}







####################################################################################################################################
####################################################################################################################################

library(here)
source(here("src/R", "config.R"))  


#### VAR FORECASTING

out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/standard_20250807_174624.rds")


out$risk$VaR[4000,,]
out2$risk$VaR[4000,,]


length(out$log_pred)

dim(out$port$VaR)



actual = readRDS("data/returns_actual.rds")[,2:4]

dim(actual)

port_actual <- rowMeans(actual)#[757:1501] 

length(port_actual)


out$port$VaR[,4]



sum(port_actual<out$port$VaR[,4])/length(port_actual)



sum(actual[757:1501,3]<out$risk$VaR[,3,4])/(1501-757)

actual[757:1501,stock]  

# Extract series
actual_series <- port_actual[2757:3500]
var_series <- out$port$VaR[, 1]  # make sure dimensions match

# Plot actual portfolio values
plot(actual_series, type = "l", col = "black", lwd = 2,
     ylab = "Value", xlab = "Time", main = "Actual vs VaR",
     ylim = range(c(actual_series, var_series)))

# Add the VaR line
lines(var_series, col = "red", lwd = 2, lty = 2)

# Optional: Add legend
legend("topright", legend = c("Actual", "VaR"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)



##### UNCERTAINTY

out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/standard_20250807_174624.rds")


plot_fam_lines_edges <- function(fam_hist, edges = c(1L, 2L),
                                 dates = NULL, use_prop = TRUE,
                                 smooth_k = NULL, tick_by = "years") {
  # global family codes across selected edges -> consistent colors/legend
  lev_all <- sort(unique(as.vector(fam_hist[ , , edges])))
  lev_all <- lev_all[!is.na(lev_all)]
  K <- length(lev_all)
  cols <- seq_len(K)
  
  op <- par(mfrow = c(length(edges), 1),
            mar = c(3.2, 4.6, 1.5, 0.5),               # <- wider left margin
            mgp = c(2.0, 0.6, 0))                      # <- title/ticks closer
  on.exit(par(op), add = TRUE)
  on.exit(par(op), add = TRUE)
  
  for (edge in edges) {
    fam <- fam_hist[ , , edge, drop = FALSE][ , , 1]   # [N x T]
    N   <- nrow(fam); Tt <- ncol(fam)
    
    # X-axis: dates or indices
    if (!is.null(dates)) {
      x <- as.Date(dates)
      if (length(x) != Tt) {
        warning("Length of 'dates' (", length(x), ") != T (", Tt, "); falling back to index.")
        x <- seq_len(Tt)
        use_dates <- FALSE
      } else {
        use_dates <- TRUE
      }
    } else {
      x <- seq_len(Tt)
      use_dates <- FALSE
    }
    
    # counts per family per time using global levels
    idx <- matrix(match(as.vector(fam), lev_all), nrow = N, ncol = Tt)
    counts <- sapply(seq_len(Tt), function(t) {
      tabulate(idx[, t][!is.na(idx[, t])], nbins = K)
    })
    
    denom <- colSums(!is.na(idx))
    props <- sweep(counts, 2, pmax(denom, 1L), "/")
    ymat  <- if (use_prop) props else counts  # K x T
    
    # optional cheap smoothing (one-sided moving average)
    if (!is.null(smooth_k) && smooth_k > 1) {
      filt <- rep(1 / smooth_k, smooth_k)
      ymat <- t(apply(ymat, 1, function(v) as.numeric(stats::filter(v, filt, sides = 1))))
    }
    
    # plot
    matplot(x, t(ymat), type = "l", lty = 1, lwd = 1, col = cols,
            xlab = if (use_dates) "Date" else "time",
            ylab = if (use_prop) "proportion" else "count",
            ylim = if (use_prop) c(0, 1) else range(ymat, na.rm = TRUE),
            xaxt = if (use_dates) "n" else "s",
            main = paste("Edge", edge))
    
    if (use_dates) {
      axis.Date(1, at = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), by = tick_by),
                format = "%Y")
    }
    legend("topright", bty = "n", legend = paste0("fam=", lev_all),
           col = cols, lty = 1, lwd = 1)
  }
}

# Example:
dates <- readRDS("data/dates.rds")
png("figures/fam_mix_edges.png", width = 9, height = 5, units = "in", res = 300)
plot_fam_lines_edges(out$fam_hist, edges = c(1L, 2L), dates = dates, use_prop = TRUE)
dev.off()




plot_edge1_conditional_params <- function(fam_edge1, p1_edge1, p2_edge1,
                                          dates = NULL, min_n = 20, tick_by = "years",
                                          ylims = NULL, clip = c(0.02, 0.98), pad = 0.05,
                                          legend_pos = "topleft") {  # <- set to NULL to hide legend
  fam <- fam_edge1; p1 <- p1_edge1; p2 <- p2_edge1
  N <- nrow(fam); Tt <- ncol(fam)
  
  use_dates <- !is.null(dates) && length(dates) == Tt
  x <- if (use_dates) as.Date(dates) else seq_len(Tt)
  
  stat_t <- function(param_col, mask, min_n) {
    v <- param_col[mask & is.finite(param_col)]
    n <- length(v); if (n < min_n) return(c(NA_real_, NA_real_, NA_real_))
    qs <- stats::quantile(v, c(0.025, 0.975), names = FALSE, type = 8)
    c(mean(v), qs[1], qs[2])
  }
  series_stats <- function(code, par_mat) {
    res <- sapply(seq_len(Tt), function(t) stat_t(par_mat[, t], fam[, t] == code, min_n))
    list(mean = res[1, ], lo = res[2, ], hi = res[3, ])
  }
  
  gauss  <- series_stats(1L, p1)
  bb1_p1 <- series_stats(2L, p1)
  bb1_p2 <- series_stats(2L, p2)
  
  choose_ylim <- function(m, lo, hi, override = NULL) {
    if (!is.null(override)) return(override)
    v <- c(m, lo, hi); v <- v[is.finite(v)]
    if (!length(v)) return(c(-1, 1))
    rng <- if (length(clip) == 2) stats::quantile(v, clip, names = FALSE, type = 8) else range(v, na.rm = TRUE)
    rng <- range(rng); d <- diff(rng); rng + c(-pad, pad) * d
  }
  get_override <- function(name) if (is.list(ylims)) ylims[[name]] else if (is.numeric(ylims)) ylims else NULL
  yg  <- choose_ylim(gauss$mean,  gauss$lo,  gauss$hi,  get_override("gauss"))
  yb1 <- choose_ylim(bb1_p1$mean, bb1_p1$lo, bb1_p1$hi, get_override("bb1_p1"))
  yb2 <- choose_ylim(bb1_p2$mean, bb1_p2$lo, bb1_p2$hi, get_override("bb1_p2"))
  
  op <- par(mfrow = c(3, 1), mar = c(3.2, 4.6, 1.5, 0.5), mgp = c(2, 0.6, 0)); on.exit(par(op), add = TRUE)
  
  plot_one <- function(xx, m, lo, hi, ylim, ylab, main) {
    plot(xx, m, type = "l", lwd = 1.2,
         xlab = if (use_dates) "Date" else "time",
         ylab = ylab, xaxt = if (use_dates) "n" else "s",
         ylim = ylim, main = main)
    lines(xx, lo, lty = 2); lines(xx, hi, lty = 2)
    if (use_dates) axis.Date(1, at = seq(min(xx, na.rm = TRUE), max(xx, na.rm = TRUE), by = tick_by), format = "%Y")
    if (!is.null(legend_pos)) legend(legend_pos, bty = "n",
                                     lty = c(1, 2), lwd = c(1.2, 1),
                                     legend = c("mean", "95% CI"))
  }
  
  plot_one(x, gauss$mean,  gauss$lo,  gauss$hi,  yg,  "Gaussian param1",  "Edge 1 · Family 1 (Gaussian)")
  plot_one(x, bb1_p1$mean, bb1_p1$lo, bb1_p1$hi, yb1, "BB1 param1 (θ)",   "Edge 1 · Family 2 (BB1) · param1")
  plot_one(x, bb1_p2$mean, bb1_p2$lo, bb1_p2$hi, yb2, "BB1 param2 (δ)",   "Edge 1 · Family 2 (BB1) · param2")
}






# Example usage (edge 1 slices):
fam1 <- out$fam_hist[ , , 1]
p1e1 <- out$par1_hist[ , , 1]   # param1 for edge 1
p2e1 <- out$par2_hist[ , , 1]   # param2 for edge 1  (as you described)
dates <- readRDS("data/dates.rds")

png("figures/param_mix_edges.png", width = 9, height = 5, units = "in", res = 300)
plot_edge1_conditional_params(
  fam1, p1e1, p2e1, dates = dates,
  ylims = list(gauss = c(0, 1), bb1_p1 = c(0, 0.7), bb1_p2 = c(1, 3))
)
dev.off()

