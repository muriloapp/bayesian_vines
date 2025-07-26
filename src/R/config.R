# ──────────────────────────────────────────────────────────────────────────────
#  C-vine SMC with spike-and-slab prior
#  2025-06-16
# ──────────────────────────────────────────────────────────────────────────────
set.seed(126)
quiet_assert <- function() {
  assignInNamespace("assert_that", function(...) invisible(TRUE), ns = "assertthat")
  assignInNamespace("see_if",      function(...) invisible(TRUE), ns = "assertthat")
}

load_packages <- function() {
  pkgs <- c(
    "rvinecopulib", "VineCopula", "data.table", "tictoc", #"Rcpp",
    "here", "parallel", "profvis", "ggplot2", "reshape2", #, "RcppThread"
    "data.table", "parallel"
  )
  lapply(pkgs, require, character.only = TRUE)
}

## build_cfg()  ─────────────────────────────────────────────────────────
build_cfg <- function(d,
                      lambda        = 1,
                      step_sd       = 0.05,
                      q_flip   = NULL,
                      K       = d * (d - 1) / 2,
                      families      = c("indep", "gaussian", "bb1"),  # NEW
                      families_first = c("indep", "gaussian", "bb1"),
                      families_deep  = c("indep", "gaussian"),
                      adapt_step_sd = TRUE) {
  
  if (is.null(q_flip))        # default that mimics “stay vs leave equal”
    q_flip <- 1 / (K + 1)
  
  list(
    d       = d,
    K       = K,
    M       = 1000,
    ess_thr = 0.50,
    W       = 1000L,
    k_step  = 1L,
    n_mh    = 3L,
    W_predict    = 2000L,
    q_flip=q_flip,
    step_sd = step_sd,
    lambda  = lambda,
    families = families,          # ← store user choice
    families_first = families_first,
    families_deep  = families_deep,
    adapt_step_sd = adapt_step_sd,
    seed    = 42, G = 2L,
    edge_tree  = edge_tree_map(d),
    nc       = max(parallel::detectCores()-1, 1),
    type     = "standard",
    alphas     = c(.05, .025)
  )
}

source(here("src/R", "core_functions.R"))
source(here("src/simulation", "simulation.R"))
source(here("src/models", "smc_stand_vine.R"))
#source(here("src/models", "smc_block_vine.R"))
source(here("src/R", "results_helpers.R"))
source(here("src/models", "main_empirical.R"))
source(here("src/R", "utils.R"))


quiet_assert()                       
load_packages()                      































