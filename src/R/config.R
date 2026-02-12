
quiet_assert <- function() {
  assignInNamespace("assert_that", function(...) invisible(TRUE), ns = "assertthat")
  assignInNamespace("see_if",      function(...) invisible(TRUE), ns = "assertthat")
}

load_packages <- function() {
  pkgs <- c(
    "rvinecopulib", "VineCopula", "data.table", "tictoc", #"Rcpp",
    "here", "parallel", "profvis", "ggplot2", "reshape2", #, "RcppThread"
    "data.table", "assertthat"
  )
  lapply(pkgs, require, character.only = TRUE)
}


build_cfg <- function(d,
                      lambda         = 1,
                      step_sd        = 0.025,
                      q_flip         = NULL,                 # <- default now computed later
                      K              = NULL,                 # <- ignored if truncation given
                      families       = c("bb1","bb1r180","bb7","bb7r180","t"),
                      families_first = c("bb1","bb1r180","bb7","bb7r180","t"),
                      #families_deep  = c("bb1","bb1r180","bb7","bb7r180","t"),
                      families_deep  = c("t"),
                      adapt_step_sd  = TRUE,
                      trunc_tree     = NULL,
                      W              = 252L,
                      W_predict      = 252
                      ) {
  
  d <- as.integer(d)
  K_full <- as.integer(d * (d - 1L) / 2L)
  
  trunc_tree <- ifelse(!is.null(trunc_tree), trunc_tree, d - 1L)
  K_trunc <- K_from_trunc(d, trunc_tree)
  
  
  edge_tree <- edge_tree_map(d, trunc_tree)                 # edge-to-tree mapping consistent with truncation (length == K_trunc)
  if (is.null(q_flip)) q_flip <- 1 / (d - 1L)               # default q_flip uses *truncated* K (prob. of proposing a family flip vs. stay)
  
  cfg <- list(
    d            = d,
    K            = K_trunc,                    # truncated K used everywhere downstream
    K_full       = K_full,                     # for reference if needed
    trunc_tree   = trunc_tree,                 # truncation level
    M            = 2000L,                      # number of particles
    ess_thr      = 0.50,                       # ESS threshold
    n_mh         = 3L,                         # number of mh iterations
    W            = W,                          # look-back period to mh_step and to compute empirical tails
    k_step       = 1000,                   # print diagnostic every k_step obs
    
    W_predict    = W_predict,                  # training period -> start to predict after W_predict
    q_flip       = q_flip,                     # prob of proposing a family flip
    step_sd      = step_sd,                    # standard deviation of the proposal
    lambda       = lambda,
    families     = families,
    families_first = families_first,
    families_deep  = families_deep,
    adapt_step_sd  = adapt_step_sd,
    adapt_iter = 0L,
    seed         = 11111L,
    G            = 2L,
    edge_tree    = edge_tree,                  # <- length == K
    nc           = 4,                          # cores
    type         = "standard",
    alphas       = c(0.1, 0.05, 0.025, 0.01),  # alpha for tail risk
    use_tail_informed_prior = FALSE,           # turn on strong tail-centered priors
    tip_method   = "EmpTC",                    # FRAPO::tdc method: "EmpTC" or "EVT"
    tip_k        = NULL,                       # if NULL, FRAPO’s default k=floor(sqrt(n)) -> bias-variance trade-off
    tip_sd_logit = 0.025,                      # variance of prior -> small sd -> STRONG prior around empirical λ

    ## Tree 1 edge map (E1 x 2), in the SAME order as edges where edge_tree==1L
    ## If you already have it, keep using yours.
    edge_pair    = NULL,
    
    
    # currently not in use --> EXCLUDE
    use_weighted_ll = FALSE,
    tauL      = 0.1,                        # per-margin lower-quantile threshold for "tail"
    joint_k   = 2L,                         # how many margins must be in the lower tail in a row
    tail_eps  = 0.30                       # weight for non-tail rows (0<eps<=1)
  )
  
  # sanity checks that prevent the “length mismatch” errors downstream
  stopifnot(length(cfg$edge_tree) == cfg$K)
  cfg
}


source(here("src/R", "utils.R"))
source(here("src/R", "core_functions.R"))
#source(here("src/models", "smc_stand_vine.R"))
#source(here("src/models", "smc_block_vine.R"))
source(here("src/R", "results_helpers.R"))
#source(here("src/models", "main_empirical.R"))
source(here("src/R", "metrics.R"))





quiet_assert()                       
load_packages()                      































