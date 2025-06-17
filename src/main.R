# ──────────────────────────────────────────────────────────────────────────────
#  C-vine SMC with spike-and-slab prior
#  2025-06-16
# ──────────────────────────────────────────────────────────────────────────────
set.seed(42)
quiet_assert <- function() {
  assignInNamespace("assert_that", function(...) invisible(TRUE), ns = "assertthat")
  assignInNamespace("see_if",      function(...) invisible(TRUE), ns = "assertthat")
}

load_packages <- function() {
  pkgs <- c(
    "rvinecopulib", "VineCopula", "data.table", "tictoc", "Rcpp",
    "here", "parallel", "RcppThread", "profvis", "ggplot2", "reshape2"
  )
  lapply(pkgs, require, character.only = TRUE)
}

build_cfg <- function(d) {
  list(
    d            = d,
    K            = d * (d - 1) / 2,
    M            = 1000,
    pi0_edge     = 0.30,
    slab_sd      = 0.50,
    ess_thr      = 0.50,
    W            = 1000L,
    k_step       = 1L,                  
    proc_sd      = 0,
    p_dyn_flip   = 0,
    n_mh         = 1L,
    step_sd      = 0.05,
    p_flip_edge  = 0.10,
    indep_copula = bicop_dist("indep"),
    W_predict    = 5L,
    seed         = 42,
    G            = 2L                      # Group in which tree
  )
}



quiet_assert()
load_packages()

source(here("src", "core_functions.R"))
source(here("src", "simulation.R"))
source(here("src", "smc_stand_vine.R"))
source(here("src", "smc_block_vine.R"))
source(here("src", "results_helpers.R"))



U  <- sim_static_cop_3(N = 6)
N  <- nrow(U)
d  <- ncol(U)
cfg <- build_cfg(d)


#results <- run_standard_smc(U, cfg, type="standard")
results <- run_block_smc(U, cfg, type="block")

cat("\n\n===== FINAL MODEL EVALUATION =====\n")
cat(sprintf("Log Model Evidence: %.4f\n", results$log_model_evidence))




#================================================================================
# Explore results
#================================================================================

# plot(results$diag_log$unique)
# 
# plot_genealogy_theta(results$theta_hist, results$ancestorIndices, edge_id = 1)   
# 
# plot_theta_histograms(results$theta_mean, k_set = c(1))
# 
# plot_theta_paths(results$theta_mean, results$theta_se, k=3, theta_true = 0.2)
# 
# 







