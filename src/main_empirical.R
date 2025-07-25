# ===========================================================================
#   empirical_run.R        tidy “full-history” SMC + risk metrics
# ===========================================================================

library(here)
source(here("src", "config.R"))          # build_cfg(), constants …


save_result <- function(res, dir_out = here("empirical_results")) {
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  fn <- file.path(dir_out,
                  sprintf("%s_%s.rds", res$alg, res$cfg$label))
  saveRDS(res, fn)
  message("✓ saved → ", fn)
}


make_skeleton <- function(U_train) {
  # any structure selection you like – here: Gaussian C-vine
  CVM <- RVineStructureSelect(U_train, familyset = c(1), type=1, indeptest = TRUE, level = 0.1)
  old_M  <- CVM$Matrix          # the VineCopula matrix you showed
  order  <- old_M[, 1]
  skeleton <- vinecop(U_train, family_set = "gaussian", structure = cvine_structure(order))
  skeleton
}

make_cluster <- function(n_cores, seed, exports) {
  cl <- makeCluster(n_cores)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, exports , envir = parent.frame())
  cl
}
# ───────────────────────────────────────────────────────────────────────────
#  0. helper: risk metrics for one forecast date
# ───────────────────────────────────────────────────────────────────────────
risk_stats <- function(R_pred,   # L × d matrix  (future PIT-draws ⇒ returns)
                       mu_fc,    # length-d mean forecast
                       sig_fc,   # length-d scale forecast
                       alphas = c(.05, .025)) {
  
  ## 1.   back-transform to returns
  R_t <- sweep(R_pred, 2, sig_fc, `*`) + rep(mu_fc, each = nrow(R_pred))
  
  mu_hat  <- colMeans(R_t)
  var_hat <- apply(R_t, 2, var)
  CI_lo   <- apply(R_t, 2, quantile, 0.025)
  CI_hi   <- apply(R_t, 2, quantile, 0.975)
  
  ## 2.   VaR & ES
  losses <- -R_t                               # ≤ 0 is a gain
  VaR <- sapply(alphas, function(a)
    apply(losses, 2, quantile, probs = a))
  ES  <- sapply(seq_along(alphas), function(k) {
    thr <- VaR[k, ]
    vapply(seq_len(ncol(losses)), \(j)
           mean(losses[losses[,j] >= thr[j], j]),
           numeric(1))
  })
  
  ## 3.   flatten to a single named vector  (for rbindlist)
  as.vector(c(mu_hat, var_hat, CI_lo, CI_hi,
              as.vector(VaR), as.vector(ES)))
}

risk_col_names <- function(d, alphas) {
  sn <- paste0("S", seq_len(d))
  c(paste0("mean_", sn),
    paste0("var_",  sn),
    paste0("ci_lo_", sn),
    paste0("ci_hi_", sn),
    unlist(lapply(alphas, \(a) sprintf("VaR%g_%s", a*100, sn))),
    unlist(lapply(alphas, \(a) sprintf("ES%g_%s",  a*100, sn))))
}

# ───────────────────────────────────────────────────────────────────────────
#  1.   run SMC  (keeps *all* history arrays)
# ───────────────────────────────────────────────────────────────────────────
smc_full <- function(data, cfg) {
  
  U      <- data$U
  mu_fc  <- data$mu_fc
  sig_fc <- data$sig_fc
  
  ## (a) structure fixed on a training window ------------------------------
  t_train <- cfg$W_predict
  skeleton <- make_skeleton(U[1:t_train, ])
  
  ## (b) cluster ------------------------------------------------------------
  exports <- c(
    ## -------------- constants & templates -----------------
    "FAM_INFO", "FAM_INDEP", "FAM_GAUSS", "FAM_BB1",
    "T_INDEP",  "T_GAUSS",   "T_BB1",
    ## -------------- helper functions ----------------------
    "active_fams", "sanitize_bb1",
    ## -------------- core SMC kernels ----------------------
    "mh_step", "mh_step_in_tree",
    "update_weights", "ESS", "systematic_resample", "resample_move",
    ## -------------- log-target & proposals ----------------
    "log_prior", "bb1_tail2par", "bb1_par2tail", "bb1_log_jacobian",
    "rtnorm_vec", "log_prior_edge",
    ## -------------- likelihood helpers --------------------
    "bicop_dist", "vinecop_dist", "dvinecop", "fast_vine_from_particle",
    ## -------------- shared data objects ------------------
    "skeleton", "cfg",
    ## -------------- diagnostics & prediction --------------
    "diagnostic_report", "compute_predictive_metrics",
    "compute_log_incr",
    ## -------------- small utilities -----------------------
    "w_mean", "w_var", "mc_se", "w_quantile"
  )
  cl <- make_cluster(cfg$nc, cfg$seed, exports)
  
  ## (c) pre-allocation -----------------------------------------------------
  M <- cfg$M; K <- cfg$K; N <- nrow(U)
  out <- list(
    log_pred    = numeric(N),
    diag_log    = data.table(t = integer(N), ESS = numeric(N),
                             unique = integer(N), euc = numeric(N),
                             sparsity = numeric(N)),
    mh_acc_pct  = rep(NA_real_, N),
    step_sd_hist= rep(NA_real_, N),
    fam_hist  = array(NA_integer_, dim = c(M,N,K)),
    par1_hist = array(NA_real_,    dim = c(M,N,K)),
    par2_hist = array(NA_real_,    dim = c(M,N,K)),
    ancestorIndices = matrix(0L, M, N),
    risk_log  = data.table()                           # will rbind() rows
  )
  
  ## particles
  particles <- replicate(M, new_particle(cfg), simplify = FALSE)
  out$ancestorIndices[,1] <- seq_len(M)
  
  ## column names for risk_log
  d <- ncol(U); risk_cols <- risk_col_names(d, cfg$alphas)
  
  # ── main loop ------------------------------------------------------------
  for (t in seq_len(N)) {
    
    u_t <- U[t,,drop=FALSE]
    
    ## 1. weight-update
    log_inc <- compute_log_incr(particles, u_t, skeleton, cfg)
    particles <- update_weights(particles, log_inc)
    w <- vapply(particles, `[[`, numeric(1), "w")
    
    ## 2. prediction block
    if (t > cfg$W_predict) {
      
      ## 2a. log predictive density
      out$log_pred[t] <- compute_predictive_metrics(
        u_t, particles, skel, w/sum(w), cfg)$log_pred_density
      
      ## 2b. 10 000 predictive draws
      draws <- smc_predictive_sample(particles, skel,
                                     w/sum(w), L = 10000, cl = cl)
      
      ## 2c. risk metrics & store
      risk_row <- as.list(risk_stats(draws,
                                     mu_fc[t,,drop=TRUE],
                                     sig_fc[t,,drop=TRUE],
                                     alphas))
      out$risk_log <- rbind(out$risk_log,
                            data.table(date  = t,
                                       !!!setNames(risk_row, risk_cols)))
    }
    
    ## 3. diagnostics
    dg <- diagnostic_report(t, 0, U, particles, w, cfg)
    out$diag_log[t, `:=`(t        = t,
                         ESS      = dg$ESS,
                         unique   = dg$unique,
                         euc      = dg$euc,
                         sparsity = dg$sparsity)]
    
    ## 4. resample / move
    if (ESS(w) < cfg$ess_thr * M && t < N) {
      newAnc <- stratified_resample(w)
      
      data_up_to_t <- U[max(1, t - cfg$W + 1):t, , drop = FALSE]
      move_out <- resample_move(particles, newAnc,
                                 data_up_to_t,
                                 cl, cfg$type, cfg, skeleton = skeleton)
      
      particles <- move_out$particles
      out$mh_acc_pct[t] <- move_out$acc_pct
      if (cfg$adapt_step_sd) {
        {cfg$step_sd <- compute_adapt_step_sd(cfg, move_out$acc_pct)}
        out$step_sd_hist[t] <- cfg$step_sd
      }
    } else {
      step_prev <- t_idx - 1L
      newAnc    <- if (step_prev < 1L) seq_len(M) else out$ancestorIndices[, step_prev]
    }
    out$ancestorIndices[, t] <- newAnc
    
    ## 5. history arrays
    out$fam_hist [ , t,] <- t(vapply(particles, `[[`, integer(K),"fam"))
    out$par1_hist[, t,]  <- t(vapply(particles, `[[`, numeric(K),"th1"))
    out$par2_hist[, t,]  <- t(vapply(particles, `[[`, numeric(K),"th2"))
  }
  
  out$particles_final    <- particles
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  out
}




# ───────────────────────────────────────────────────────────────────────────
#  2.  driver over cfg variants
# ───────────────────────────────────────────────────────────────────────────
run_empirical <- function() {
  
  dat <- list(
    U      = readRDS("data/PIT.rds"),
    mu_fc  = readRDS("data/returns_mean_forecast.rds")[,-1],  # drop date col
    sig_fc = readRDS("data/returns_vol_forecast.rds")[,-1]
  )
  
  cfg_variants <- list(list(label = "test"))          # add as you like
  for (v in cfg_variants) {
    cfg <- modifyList(build_cfg(ncol(dat$U)), v[ setdiff(names(v),"label") ])
    cfg$label <- v$label %||% "cfg"
    res <- smc_full(dat, cfg)
    res$alg <- "standard"
    save_result(res)
  }
}

# ── fire -------------------------------------------------------------------
run_empirical()



















