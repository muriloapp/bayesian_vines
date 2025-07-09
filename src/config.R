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

build_cfg <- function(d, tau_prior = c("fixed", "inv_gamma"),
                      tau0      = 0.025,      # ≈ old spike_sd
                      c_slab    = 37,      # so c*tau0 ≈ old slab_sd = 0.50
                      a0 = 3, b0 = (a0 - 1) * tau0^2,   # hyper-pars if IG
                      pi_prior = c("fixed","beta"),
                      pi0   = 0.50,                     # centre if fixed
                      a_pi  = 2, b_pi = 2,
                      adapt_step_sd = TRUE,
                      n_fam = 2
  ) {
  tau_prior <- match.arg(tau_prior)
  
  list(
    d            = d,
    K            = d * (d - 1) / 2,
    M            = 1000,
    ess_thr      = 0.50,
    W            = 1000L,     
    k_step       = 1L,                  
    proc_sd      = 0,
    p_dyn_flip   = 0,
    n_mh         = 3L,
    step_sd      = 0.05,
    
    #RJCMCM
    pi0_edge     = 0.50,
    slab_sd      = 0.50,
    p_flip_edge  = 0.25,
    indep_copula = bicop_dist("indep"),
    p_birth      = 0.20,   # prob of proposing a birth vs. death when edge set ≠ {0,K}
    p_within     = 0.60,   # prob of proposing a RW inside the current model
    p_death      = 0.20,  
    q_theta_sd   = 0.50, 
    
    W_predict    = 5L,
    seed         = 42,
    G            = 2L,                      # Group in which tree
    adapt_step_sd = adapt_step_sd,
    pi_prior     = pi_prior,
    pi0          = pi0,                     # centre if fixed
    a_pi         = a_pi, 
    b_pi         = b_pi,                    # Beta hyper-pars

    c_slab    = c_slab,       
    tau_prior = tau_prior,     
    tau0      = tau0,          
    a0        = a0, b0 = b0,
    
    n_fam  = n_fam,                      # Gaussian, Frank, Clayton, Gumbel …
    omega  = rep(1/n_fam, n_fam),   # π(m) – eq. (6)
    p_switch = 0.30,
    
    delta_spike <- 0.05      # |θ| < δ → forced spike
    #eps_overlap <- 1e-3 
  )
}


source(here("src", "core_functions.R"))
source(here("src", "core_functions_rjmcmc.R"))
source(here("src", "simulation.R"))
source(here("src", "smc_stand_vine.R"))
source(here("src", "smc_block_vine.R"))
source(here("src", "smc_vines_rjmcmc.R"))
source(here("src", "results_helpers.R"))

quiet_assert()                       
load_packages()                      































