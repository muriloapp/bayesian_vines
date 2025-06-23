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
    "here", "parallel", "profvis", "ggplot2", "reshape2"#, "RcppThread"
  )
  lapply(pkgs, require, character.only = TRUE)
}

build_cfg <- function(d) {
  list(
    d            = d,
    K            = d * (d - 1) / 2,
    M            = 1000,
    pi0_edge     = 0.50, #0.3
    slab_sd      = 0.50,
    ess_thr      = 0.50,
    W            = 1000L,      #1000L,
    k_step       = 1L,                  
    proc_sd      = 0,
    p_dyn_flip   = 0,
    n_mh         = 3L,
    step_sd      = 0.05,
    p_flip_edge  = 0.25,
    indep_copula = bicop_dist("indep"),
    W_predict    = 5L,
    seed         = 126, #42
    G            = 2L                      # Group in which tree
  )
}


source(here("src", "core_functions.R"))
source(here("src", "simulation.R"))
source(here("src", "smc_stand_vine.R"))
source(here("src", "smc_block_vine.R"))
source(here("src", "results_helpers.R"))

quiet_assert()                       # your helper that silences assertthat
load_packages()                      # loads rvinecopulib, here(), … (already defined)
































