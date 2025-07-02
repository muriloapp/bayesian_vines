# ---------------------------------------------------------------
#  6-D C-vine simulator   (Gaussian pair copulas + random zeros)
# ---------------------------------------------------------------
library(VineCopula)           # RVineMatrix(), RVineSim()

sim_static_cop_6 <- function(N      = 200,
                             p_zero = 0.5,
                             rho_lo = -0.99,
                             rho_hi =  0.99) {
    d <- 6
  
  ## 1. C-vine structure matrix (same pattern as your 3-D prototype)
  sim_matrix <- matrix(0, d, d)
  for (j in 1:d)
    sim_matrix[j:d, j] <- j:d         # column j :  j, j+1, …, d
  
  ## 2. Family- and parameter matrices (lower-triangular part only)
  family_matrix <- matrix(0, d, d)
  theta_matrix  <- matrix(0, d, d)
  
  for (j in 1:(d - 1)) {              # columns
    for (i in (j + 1):d) {            # rows below the diagonal
      if (runif(1) > p_zero) {        # keep the edge?
        family_matrix[i, j] <- 1      # 1 = Gaussian copula
        theta_matrix[i,  j] <- runif(1, rho_lo, rho_hi)
      }
    }
  }
  
  ## 3. Build the R-vine object and simulate data
  RVM <- RVineMatrix(sim_matrix,
                     family = family_matrix,
                     par    = theta_matrix)
  
  U <- RVineSim(N, RVM)
  colnames(U) <- paste0("U", 1:d)
  
  ## 4. Return both the data and the true specification
  list(U      = U,           # simulated observations
       RVM    = RVM,         # full C-vine object
       family = family_matrix,
       theta  = theta_matrix)
}


# out <- sim_static_cop_6(N = 200, seed = 2)
# print(out$RVM)        # concise summary of structure, families, and θ
# out$family            # 6 × 6 family matrix (0 = independence, 1 = Gaussian)
# out$theta             # 6 × 6 correlation matrix (zeros on skipped edges)
# head(out$U)           # your data





# ──────────────────────────────────────────────────────────────────────────────
#  mc_experiment.R  ── full Monte-Carlo for C-vine SMC  ── 2025-07-02
# ──────────────────────────────────────────────────────────────────────────────

library(here)
source(here("src", "config.R"))          # build_cfg(), sim_static_cop_6(), …
# source(here("src", "smc_kernels.R"))   # run_standard_smc(), run_block_smc()

## -------- 1. experiment specs ------------------------------------------------
n_sim     <- 100        # Monte-Carlo replications
N_obs     <- 300        # sample size per replication
seed_base <- 20250618   # reproducible but different for each sim

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
  message(sprintf("✓ saved %s | %s | %s", label, alg, sim_tag))
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
    
    for (alg in c("standard", "block")) {
      res <- switch(
        alg,
        standard = run_standard_smc(data, cfg, type = "standard"),
        block    = run_block_smc(data,    cfg, type = "block")
      )
      res$cfg      <- cfg
      res$sim_tag  <- sim_tag
      res$alg      <- alg
      save_fit(res, label, alg, sim_tag)
    }
  }
  
  ## optional -- keep memory tidy between reps
  rm(data); gc(verbose = FALSE)
}

