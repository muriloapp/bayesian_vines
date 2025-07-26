
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

























