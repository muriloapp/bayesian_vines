# ──────────────────────────────────────────────────────────────────────────────
#  mc_experiment.R  ── full Monte-Carlo for C-vine SMC  ── 2025-07-02
# ──────────────────────────────────────────────────────────────────────────────

library(here)
source(here("src", "config.R"))          # build_cfg(), sim_static_cop_6(), …
# source(here("src", "smc_kernels.R"))   # run_standard_smc(), run_block_smc()

## -------- 1. experiment specs ------------------------------------------------
n_sim     <- 5        # Monte-Carlo replications
N_obs     <- 1000        # sample size per replication
seed_base <- 42   # reproducible but different for each sim

cfg_variants <- list(
  list(
    label      = "M200_ivgamma_beta",
    tau_prior  = "inv_gamma",
    pi_prior   = "beta",
    M          = 2000
  )
  ## add more variants here
)

## -------- 2. utilities -------------------------------------------------------
out_root <- here("simul_results", "mc_static_dgp")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

save_fit <- function(res, label, alg, sim_tag) {
  dir.create(file.path(out_root, label), showWarnings = FALSE)
  saveRDS(res,
          file.path(out_root, label, sprintf("%s_%s.rds", alg, sim_tag)))
  message(sprintf("saved %s | %s | %s", label, alg, sim_tag))
}

## -------- 3. Monte-Carlo loop ------------------------------------------------
for (sim in seq_len(n_sim)) {
  set.seed(seed_base + sim)
  sim_tag <- sprintf("sim%03d", sim)
  
  ## 3a. simulate data and store it once --------------------------------------
  data <- sim_static_cop_6(N = N_obs)
  d    <- ncol(data$U)
  data_dir <- file.path(out_root, "_data")
  dir.create(data_dir, showWarnings = FALSE)
  saveRDS(data, file.path(data_dir,
                          sprintf("data_%s.rds", sim_tag)))
  
  ## 3b. loop over cfg variants and algorithms --------------------------------
  for (v in cfg_variants) {
    label  <- v$label
    tweaks <- v[ setdiff(names(v), "label") ]
    cfg    <- modifyList(build_cfg(d), tweaks)
    cfg$label <- label
    
    for (alg in c("standard","block")) {
      
      t_sec <- system.time(
      res <- switch(
        alg,
        standard = run_standard_smc(data, cfg, type = "standard"),
        block    = run_block_smc(data,    cfg, type = "block")
      )
      )[["elapsed"]]                       # elapsed wall-clock time (s)
      res$cfg      <- cfg
      res$sim_tag  <- sim_tag
      res$alg      <- alg
      res$elapsed <- t_sec 
      
      save_fit(res, label, alg, sim_tag)
    }
  }
  
  ## optional -- keep memory tidy between reps
  rm(data); gc(verbose = FALSE)
}





# library(here)
# 
# # all block-SMC fits for the “M200_ivgamma_beta” variant
# block_files <- list.files(
#   here("simul_results/mc_static_dgp/M200_ivgamma_beta"),
#   pattern = "^standard_.*\\.rds$", full.names = TRUE
# )
# res <- lapply(block_files, readRDS)
# 
# # corresponding data sets
# data_files <- list.files(
#   here("simul_results/mc_static_dgp/_data"),
#   pattern = "^data_.*\\.rds$", full.names = TRUE
# )
# sim_data <- lapply(data_files, readRDS)
# 
# 
# 
# sim_data[[1]]$RVM
# 
# 
# theta_mean <- apply(res[[1]]$theta_hist, MARGIN = c(2, 3), FUN = mean)
# theta_sd <- apply(res[[1]]$theta_hist, MARGIN = c(2, 3), FUN = sd)
# 
# 
# 
# plot_theta_paths(tanh(theta_mean), theta_sd, k=12, theta_true = -0.37)
# 
# 








# 
# pkgs <- unique(renv::dependencies()$Package)
# 
# # create a vector of "package==version"
# lines <- sapply(pkgs, function(pkg) {
#   ver <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "not_installed")
#   paste0(pkg, "==", ver)
# })
# 
# writeLines(lines, "package_versions.txt")
# 






















