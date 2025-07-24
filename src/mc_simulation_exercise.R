# ──────────────────────────────────────────────────────────────────────────────
#  mc_experiment.R  ── full Monte-Carlo for C-vine SMC  ── 2025-07-02
# ──────────────────────────────────────────────────────────────────────────────

library(here)
source(here("src", "config.R"))          # build_cfg(), sim_static_cop_6(), …
# source(here("src", "smc_kernels.R"))   # run_standard_smc(), run_block_smc()

## ---------- experiment specs -------------------------------------------------
n_sim     <- 10
N_obs     <- 1500
seed_base <- 42
dims      <- c(3, 6)                     # <-- asked for both

cfg_variants <- list(
  list(label = "M2000", M = 2000),
  list(label = "M1000", M = 1000)
)

## ---------- helpers ----------------------------------------------------------
out_root <- here("simul_results", "mc_static_dgp")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

save_fit <- function(res, label, alg, sim_tag, d) {
  path <- file.path(out_root, sprintf("d%d", d), label)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file.path(path, sprintf("%s_%s.rds", alg, sim_tag)))
  message(sprintf("saved   d=%d | %s | %s | %s",
                  d, label, alg, sim_tag))
}

## ---------- Monte-Carlo ------------------------------------------------------
for (sim in seq_len(n_sim)) {
  for (d in dims) {
    
    set.seed(seed_base + 100 * d + sim)          # unique per (d, sim)
    sim_tag <- sprintf("sim%03d", sim)
    
    ## 1️⃣  simulate data once per replication
    data <- sim_static_cop(d = d, N = N_obs)
    data_dir <- file.path(out_root, sprintf("d%d", d), "_data")
    dir.create(data_dir, showWarnings = FALSE)
    saveRDS(data, file.path(data_dir,
                            sprintf("data_%s.rds", sim_tag)))
    
    ## 2️⃣  run all cfg variants
    for (v in cfg_variants) {
      label  <- v$label
      cfg    <- modifyList(build_cfg(d), v)      # M overwritten
      cfg$label <- label
      
      ## === choose algorithms here =====================================
      for (alg in "standard") {
        t_sec <- system.time(
          res <- run_standard_smc(data, cfg, type = "standard")
        )[["elapsed"]]
        
        res$cfg      <- cfg
        res$sim_tag  <- sim_tag
        res$alg      <- alg
        res$elapsed  <- t_sec
        res$d        <- d
        
        save_fit(res, label, alg, sim_tag, d)
      }
    }
    
    rm(data); gc(verbose = FALSE)
  }
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






















