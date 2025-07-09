
library(here)
source(here::here("src","config.R"))


######################################################################

res_stand <- readRDS("simul_results/static_dgp/standard_stat_3.rds")
res_block <- readRDS("simul_results/static_dgp/block_stat_3.rds")

res_stand$log_model_evidence
res_block$log_model_evidence

res_stand <- readRDS("simul_results/static_dgp/rjmcmc_3_RJMCMC")
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





tic()
for (i in 1:1000){j= i+1
print(j)}
toc()




mean_theta <- apply(res_stand$theta_hist, c(2, 3), mean)


plot_theta_paths(results$theta_mean, results$theta_se, k=1, theta_true = 0.6)


dim(res_stand$model_hist)













## --------------------------------------------------
##  Inputs (already in your workspace)
##  res_stand$is_slab_hist : 1000 × 20 × 3   (0 = spike, 1 = slab)
##  res_stand$model_hist   : 1000 × 20 × 3   (1 = Gaussian, 2 = Clayton)
##  res_stand$theta_hist   : 1000 × 20 × 3   (θ values)
## --------------------------------------------------

## 1.  Derive a single “family code” array.
##     0 = independent, 1 = Gaussian, 2 = Clayton
family <- ifelse(res_stand$is_slab_hist == 0L,
                 0L,                               # spike → independent
                 res_stand$model_hist)             # slab  → use model id

N <- dim(family)[1]     # number of particles  (= 1000)
T <- dim(family)[2]     # time points          (= 20)
P <- dim(family)[3]     # parameters           (= 3)

## 2.  Proportion of each family, for every (t, p) pair.
prop_family <- array(NA_real_, c(T, P, 3),
                     dimnames = list(time   = seq_len(T),
                                     param  = seq_len(P),
                                     family = c("indep","gauss","clayton")))

for (t in seq_len(T)) {
  for (p in seq_len(P)) {
    counts <- tabulate(family[, t, p] + 1L, nbins = 3)  # +1 → 1-based ids
    prop_family[t, p, ] <- counts / N
  }
}

## 3.  Conditional means of θ for each family.
mean_theta <- array(NA_real_, c(T, P, 3),
                    dimnames = dimnames(prop_family))

for (t in seq_len(T)) {
  for (p in seq_len(P)) {
    for (f in 0:2) {                    # 0 = indep, 1 = gauss, 2 = clayton
      idx <- family[, t, p] == f
      if (any(idx)) {
        mean_theta[t, p, f + 1] <- mean(res_stand$theta_hist[idx, t, p])
      }
    }
  }
}

## 4.  (Optional) Split outputs into convenient lists of 20×3 matrices
prop_list  <- lapply(1:3, \(k) prop_family[ , , k])
mean_list  <- lapply(1:3, \(k) mean_theta[ , , k])
names(prop_list) <- names(mean_list) <- c("indep", "gauss", "clayton")

##  prop_list[["gauss"]][t, p]  is the proportion of Gaussian particles
##  mean_list[["clayton"]][t, p] is E[θ | Clayton]  at time t, parameter p



res_stand$model_hist[,,]

mean_theta <- apply(res_stand$model_hist, c(2, 3), mean)



apply(res_stand$gamma_hist, c(2, 3), mean)














