
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


build_cfg <- function(d,
                      lambda         = 1,
                      step_sd        = 0.05,
                      q_flip         = NULL, #0.66, # stay or 2 possible families to go 
                      K              = d * (d - 1) / 2,
                      families       = c("indep", "gaussian", "bb1"),  # NEW
                      families_first = c("indep", "gaussian", "bb1"),
                      families_deep  = c("indep", "gaussian"),
                      adapt_step_sd  = TRUE) {
  
  if (is.null(q_flip))        # default that mimics “stay vs leave equal”
    q_flip <- 1 / (K + 1)
  
  list(
    d       = d,
    K       = K,
    M       = 1000,                                   # number of particles
    ess_thr = 0.50,                                   # threshold for resample move
    W       = 252L,                                  # rolling-window
    k_step  = 1L,                                     # print diagnostic every k_step  
    n_mh    = 3L,                                     # number mh moves
    W_predict = 756L,                                 # training period     
    q_flip=q_flip,
    step_sd = step_sd,
    lambda  = lambda,
    families = families,                              # store user choice
    families_first = families_first,
    families_deep  = families_deep,
    adapt_step_sd = adapt_step_sd,
    seed    = 42, G = 2L,
    edge_tree  = edge_tree_map(d),
    nc       = max(parallel::detectCores()-1, 1),
    type     = "standard",
    alphas     = c(0.1, .05, .025, 0.01)
  )
}

source(here("src/R", "utils.R"))
source(here("src/R", "core_functions.R"))
source(here("src/simulation", "simulation.R"))
source(here("src/models", "smc_stand_vine.R"))
#source(here("src/models", "smc_block_vine.R"))
source(here("src/R", "results_helpers.R"))
source(here("src/models", "main_empirical.R"))
source(here("src/R", "metrics.R"))



quiet_assert()                       
load_packages()                      































