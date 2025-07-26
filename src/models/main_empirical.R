
smc_full <- function(data, cfg) {
  
  U      <- data$U
  mu_fc  <- data$mu_fc
  sig_fc <- data$sig_fc
  
  t_train <- cfg$W_predict
  skeleton <- make_skeleton_CVM(U[1:t_train, ])
  
  exports <- c(
    # constants & templates 
    "FAM_INFO", "FAM_INDEP", "FAM_GAUSS", "FAM_BB1",
    "T_INDEP",  "T_GAUSS",   "T_BB1",
    # helper functions 
    "active_fams", "sanitize_bb1", "mh_worker_standard", "mh_worker_block",
    # core SMC kernels 
    "mh_step", "mh_step_in_tree",
    "update_weights", "ESS", "systematic_resample", "resample_move",
    # log-target & proposals 
    "log_prior", "bb1_tail2par", "bb1_par2tail", "bb1_log_jacobian",
    "rtnorm_vec", "log_prior_edge",
    # likelihood helpers 
    "bicop_dist", "vinecop_dist", "dvinecop", "fast_vine_from_particle",
    "rvinecop",
    # shared data objects 
    "skeleton", "cfg",
    # diagnostics & prediction 
    "diagnostic_report", "compute_predictive_metrics",
    "compute_log_incr",
    # small utilities 
    "w_mean", "w_var", "mc_se", "w_quantile"
  )
  cl <- make_cluster(cfg$nc, cfg$seed, exports)
  
  # pre-allocation 
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers    <- colnames(U); A <- length(cfg$alphas)
  risk_cols <- risk_col_names(tickers, cfg$alphas)

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
    risk    =  list(
                      dates = integer(n_oos),                                        # 1‑D
                      mean  = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
                      var   = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
                      ci_lo = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
                      ci_hi = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
                      VaR   = array(NA_real_, dim = c(n_oos, d, A),
                                    dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
                      ES    = array(NA_real_, dim = c(n_oos, d, A),
                                    dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
                    )                     # will rbind() rows
  )
  
  # init particles
  particles <- replicate(M, new_particle(cfg), simplify = FALSE)
  out$ancestorIndices[,1] <- seq_len(M)
  

  for (t in seq_len(N)) {
    u_t <- U[t,,drop=FALSE]
    
    # weight-update
    log_inc <- compute_log_incr(particles, u_t, skeleton, cfg)
    particles <- update_weights(particles, log_inc)
    w <- vapply(particles, `[[`, numeric(1), "w")
    
    # prediction block
    if (t > cfg$W_predict) {
      # log predictive density, predictive draws, risk metrics
      idx <- t - cfg$W_predict
      out$log_pred[t] <- compute_predictive_metrics(u_t, particles, skeleton, w/sum(w), cfg)$log_pred_density
      draws <- smc_predictive_sample(particles, skeleton, w/sum(w), L = 10000, cl = cl)
      rs <- risk_stats_full(            # same helper as before
        draws,
        mu_fc[t, , drop = TRUE],
        sig_fc[t, , drop = TRUE],
        cfg$alphas)

      out$risk$dates[idx]   <- t
      out$risk$mean [idx, ] <- rs$mean
      out$risk$var  [idx, ] <- rs$var
      out$risk$ci_lo[idx, ] <- rs$ci["lo", ]
      out$risk$ci_hi[idx, ] <- rs$ci["hi", ]
      out$risk$VaR [idx, , ] <- rs$VaR          # dim = d × A
      out$risk$ES  [idx, , ] <- rs$ES
    }
    
    # diagnostics
    dg <- diagnostic_report(t, 0, U, particles, w, cfg)
    out$diag_log[t, `:=`(t        = t,
                         ESS      = dg$ESS,
                         unique   = dg$unique,
                         euc      = dg$euc,
                         sparsity = dg$sparsity)]
    
    # resample / move
    if (ESS(w) < cfg$ess_thr * M && t < N) {
      newAnc <- stratified_resample(w)
      data_up_to_t <- U[max(1, t - cfg$W + 1):t, , drop = FALSE]
      move_out <- resample_move(particles, newAnc, data_up_to_t,
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
    
    # history arrays
    out$fam_hist [ , t,] <- t(vapply(particles, `[[`, integer(K),"fam"))
    out$par1_hist[, t,]  <- t(vapply(particles, `[[`, numeric(K),"th1"))
    out$par2_hist[, t,]  <- t(vapply(particles, `[[`, numeric(K),"th2"))
  }
  
  out$particles_final    <- particles
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  out
}


risk_stats <- function(R_pred,   # L × d matrix  (future PIT-draws ⇒ returns)
                       mu_fc,    # length-d mean forecast
                       sig_fc,   # length-d scale forecast
                       alphas = c(.05, .025)) {
  
  mu_fc <- as.numeric(mu_fc)
  sig_fc <- as.numeric(sig_fc)
  
  # 1.   back-transform to returns
  R_t <- sweep(R_pred, 2, sig_fc, `*`) + rep(mu_fc, each = nrow(R_pred))
  
  mu_hat  <- colMeans(R_t)
  var_hat <- apply(R_t, 2, var)
  CI_lo   <- apply(R_t, 2, quantile, 0.025)
  CI_hi   <- apply(R_t, 2, quantile, 0.975)
  
  # 2.   VaR & ES
  losses <- -R_t                               # ≤ 0 is a gain
  VaR <- sapply(alphas, function(a)
    apply(losses, 2, quantile, probs = a))
  ES  <- sapply(seq_along(alphas), function(k) {
    thr <- VaR[k, ]
    vapply(seq_len(ncol(losses)), \(j)
           mean(losses[losses[,j] >= thr[j], j]),
           numeric(1))
  })
  
  # 3.   flatten to a single named vector  (for rbindlist)
  as.vector(c(mu_hat, var_hat, CI_lo, CI_hi,
              as.vector(VaR), as.vector(ES)))
}

risk_stats_full <- function(R_pred,           # L × d matrix  (future PIT‑draws ⇒ returns)
                            mu_fc,            # length‑d numeric
                            sig_fc,           # length‑d numeric
                            alphas = c(.05, .025))
{
  mu_fc  <- as.numeric(mu_fc)
  sig_fc <- as.numeric(sig_fc)
  
  ## 1. back‑transform to returns  ---------------------------------------------
  R_t <- sweep(R_pred, 2, sig_fc, `*`) + rep(mu_fc, each = nrow(R_pred))
  
  mu_hat  <- colMeans(R_t)
  var_hat <- apply(R_t, 2, var)
  CI_lo   <- apply(R_t, 2, quantile, 0.025)
  CI_hi   <- apply(R_t, 2, quantile, 0.975)
  
  ## 2.  VaR & ES  -------------------------------------------------------------
  losses <- R_t
  VaR <- sapply(alphas, function(a)
    apply(losses, 2, quantile, probs = a))  
  ES  <- sapply(seq_along(alphas), function(k) {       # loop over α’s
    vapply(seq_len(ncol(losses)), function(j) {  # loop over assets
      thr <- VaR[j, k]                           # ← correct index
      mean(losses[losses[, j] >= thr, j])        # tail average
    }, numeric(1))
  })
  
  ## 3.  Return as a clean list  ----------------------------------------------
  list(
    mean = mu_hat,                       # length‑d
    var  = var_hat,                      # length‑d
    ci   = rbind(lo = CI_lo, hi = CI_hi),  # 2 × d
    VaR  = VaR,                          # d × length(alphas)
    ES   = ES                            # d × length(alphas)
  )
}














