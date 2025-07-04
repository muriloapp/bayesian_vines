
library(here)
source(here::here("src","config.R"))


######################################################################

res_stand <- readRDS("simul_results/static_dgp/standard_stat_3.rds")
res_block <- readRDS("simul_results/static_dgp/block_stat_3.rds")

res_stand$log_model_evidence
res_block$log_model_evidence




res_stand <- readRDS("simul_results/static_dgp/standard_8_pi07_M2000_G3_nmh5")
res_block <- readRDS("simul_results/static_dgp/block_8_pi07_M2000_G3_nmh5")

res_stand$


plot_theta_paths(res_block$theta_mean, res_block$theta_se, k=1, theta_true = 0.5)




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

# plot(results$diag_log$unique, type="l")
#  
# plot_genealogy_theta(results$theta_hist, results$ancestorIndices, edge_id = 1, ylim = c(0.5, 1.1))   
#  
# plot_theta_histograms(results$theta_mean, k_set = c(1))
# 
# plot_theta_paths(results$theta_mean, results$theta_se, k=1, theta_true = 0.6)
#  
#  


#results <- readRDS("simul_results/block_stat_3.rds")












