
# Block vine Block_vine
# Gene tree

# ================================================================
#  SMC para C-vine com spike-and-slab 
# ================================================================
library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(Rcpp)
library(here)
library(parallel)
library(RcppThread)
library(assertthat)
library(profvis)

assignInNamespace("assert_that", function(...) invisible(TRUE), ns = "assertthat")
assignInNamespace("see_if", function(...) invisible(TRUE), ns = "assertthat")

source(here('src','core_functions.R'))
source(here('src','simulation.R'))

run_standard_smc <- function(data,
                    cfg,
                    type       = c("standard", "block"),
                    n_cores    = max(parallel::detectCores()-1, 1)) {
  
  U <- data$U
  
  type <- match.arg(type)
  skeleton  <- vinecop(U, family_set = "gaussian", structure = data$RVM$Matrix[nrow(data$RVM$Matrix):1, ])

  # ── dimensions ───────────────────────────────────────────────────
  N <- nrow(U); K <- cfg$K; M <- cfg$M

  # ── pre-allocate diagnostics ───────────────────────────────────────────────
  out <- list(
    log_pred   = numeric(N),
    theta_mean = matrix(NA_real_, N, K),
    theta_se   = matrix(NA_real_, N, K),
    gamma_mean = matrix(NA_real_, N, K),
    gamma_se   = matrix(NA_real_, N, K),
    diag_log   = data.table::data.table(
      t      = integer(N),
      tr     = integer(N),
      ESS    = numeric(N),
      unique = integer(N),
      euc    = numeric(N),
      sparsity = numeric(N)
    ),
    mh_acc_pct      = rep(NA_real_, N),
    step_sd_hist      = rep(NA_real_, N),
    theta_hist      = array(NA_real_,    dim = c(M, N, K)),
    gamma_hist      = array(NA_integer_, dim = c(M, N, K)),
    ancestorIndices = matrix(0L, M, N),
    incl_hist = matrix(NA_real_, N, K),
    theta_q025   = matrix(NA_real_, N, K),
    theta_q975   = matrix(NA_real_, N, K)
  )
  
  # ── initial state ──────────────────────────────────────────────────────────
  particles <- replicate(M, new_particle(cfg), simplify = FALSE)
  
  out$ancestorIndices[, 1] <- seq_len(M)
  tr  <- 0L        #  “tree” counter — keep if you resample by tree else 0
  pos <- 1L        #  row pointer for diag_log
  
  # ── parallel backend ───────────────────────────────────────────────────────
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)          # safe cleanup
  parallel::clusterSetRNGStream(cl, cfg$seed)
  
  parallel::clusterExport(
    cl,
    c(
      ## core SMC kernels
      "mh_step", "mh_step_in_tree", #, "propagate_particles",
      "update_weights", "ESS", "systematic_resample", "resample_move",
      ## log-target & proposals
      "log_prior", "fast_vine_from_particle", "rtnorm_vec",
      ## likelihood helpers from rvinecopulib
      "bicop_dist", "vinecop_dist", "dvinecop",
      ## shared data objects
      "skeleton", "cfg",
      ## diagnostics & prediction
      "diagnostic_report", "compute_predictive_metrics", "compute_log_incr",
      ## small utilities needed inside workers
      "w_mean", "w_var", "mc_se", "w_quantile"
    ),
    envir = environment()
  )
  
  
  # ── main SMC loop ──────────────────────────────────────────────────────────
  for (t_idx in seq_len(N)) {
    u_row <- U[t_idx, , drop = FALSE]
    
    # 1. propagate step ──────────────────────────────────────────────────────
    #particles <- propagate_particles(particles, cfg)
    
    # 1. predictive metrics (after burn-in) ──────────────────────────────────
    if (t_idx > cfg$W_predict) {
      w_prev <- vapply(particles, `[[`, numeric(1), "w")
      w_prev <- w_prev / sum(w_prev)
      
      pm <- compute_predictive_metrics(u_row, particles,
                                       skeleton, w_prev, cfg)
      
      out$log_pred[t_idx]    <- pm$log_pred_density
      out$theta_mean[t_idx,] <- pm$theta_mean
      out$theta_se[t_idx,]   <- pm$theta_se
      out$gamma_mean[t_idx,] <- pm$gamma_mean
      out$gamma_se[t_idx,]   <- pm$gamma_se
    }
    
    # 1. weight update ──────────────────────────────────────────────────────
    log_incr  <- compute_log_incr(particles, u_row, skeleton, cfg)
    particles <- update_weights(particles, log_incr)
    w_new     <- vapply(particles, `[[`, numeric(1), "w")
    
    ## 3. diagnostics ------------------------------------------
    dg <- diagnostic_report(t_idx, tr, U, particles, w_new, cfg)
    out$diag_log[pos, `:=`(
      t        = t_idx,
      tr       = tr,
      ESS      = dg$ESS,
      unique   = dg$unique,
      euc      = dg$euc,
      sparsity = dg$sparsity
    )]
    pos <- pos + 1L
    
    # 5. resample + move if ESS below threshold ─────────────────────────────
    if (ESS(w_new) < cfg$ess_thr * M && t_idx < N) {
      data_up_to_t <- U[max(1, t_idx - cfg$W + 1):t_idx, , drop = FALSE]
      newAnc <- stratified_resample(w_new)
      move_out     <- resample_move(particles, newAnc, data_up_to_t,
                                    cl, type, cfg, skeleton=skeleton)
      particles <- move_out$particles
      out$mh_acc_pct[t_idx] <- move_out$acc_pct
      
      if (cfg$adapt_step_sd) {
      {cfg$step_sd <- compute_adapt_step_sd(cfg, move_out$acc_pct)}
      out$step_sd_hist[t_idx] <- cfg$step_sd
      }
    } else {
      step_prev <- t_idx - 1L
      newAnc    <- if (step_prev < 1L) seq_len(M) else out$ancestorIndices[, step_prev]
    }
    out$ancestorIndices[, t_idx] <- newAnc
    
    # 6. save history ───────────────────────────────────────────────────────
    out$theta_hist[, t_idx, ] <- t(vapply(particles, `[[`, numeric(K), "theta"))
    out$gamma_hist[, t_idx, ] <- t(vapply(particles, `[[`, integer(K), "gamma"))
    
    out$theta_q025[t_idx, ] <- dg$edges$q025[[1]]
    out$theta_q975[t_idx, ] <- dg$edges$q975[[1]]
  }
  
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  out$particles_final    <- particles
  return(out)
}
  

