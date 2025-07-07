
library(here)
source(here::here("src","config.R"))


######################################################################

res_stand <- readRDS("simul_results/static_dgp/standard_stat_3.rds")
res_block <- readRDS("simul_results/static_dgp/block_stat_3.rds")

res_stand$log_model_evidence
res_block$log_model_evidence

res_stand <- readRDS("simul_results/static_dgp/standard_6_M200_ivgamma_beta_l_2000_alt4")
res_SS <- readRDS("simul_results/static_dgp/standard_3_pi05_nmh1_N250_SSVS")

res_stand$theta_mean


plot_theta_paths(res_stand$theta_mean, res_stand$theta_se, k=1, theta_true = 0.6)
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


################################################################################################
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


mean(res_stand$diag_log$ESS)


res_stand$mh_acc_pct


res_stand$diag_log$unique



theta_mean <- apply(res_stand$theta_hist, MARGIN = c(2, 3), FUN = mean)
theta_sd <- apply(res_stand$theta_hist, MARGIN = c(2, 3), FUN = sd)



plot_theta_paths(tanh(theta_mean), 2*theta_sd, k=1, theta_true = 0)





theta_matrix  <- matrix(c(0,0.4, 0.6, 0,0,0.4, 0,0,0), 3, 3, byrow = FALSE)
RVineGrad(U, data$RVM, theta_matrix)



data$RVM












# In 2D, we only have one pair of variables (1, 2) and one copula.
# The structure matrix is fixed.
sim_matrix    <- matrix(c(2,1, 0,1), 2, 2, byrow = FALSE)
family_matrix <- matrix(c( 0,1, 0,0), 2, 2, byrow = FALSE)
theta_matrix  <- matrix(c( 0,0.4, 0,0), 2, 2, byrow = FALSE)

# Construct the 2D RVineMatrix object.
RVM <- RVineMatrix(Matrix = sim_matrix, 
                   family = family_matrix,
                   par = theta_matrix, 
                   names = c("V1", "V2"))

part  <- matrix(c( 0, 0.20, 0, 0), 2, 2, byrow = FALSE)

out2 <- RVineGrad(simdata, RVM, par = part)
out2$gradient/nrow(simdata)

