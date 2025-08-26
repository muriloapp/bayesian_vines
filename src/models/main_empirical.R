

smc_full <- function(data, cfg) {
  
  U      <- data$U
  mu_fc  <- data$mu_fc
  sig_fc <- data$sig_fc
  df_fc     <- data$df_fc 
  shape_fc     <- data$shape_fc 
  
  y_real <- data$y_real
  
  
  t_train <- cfg$W_predict
  skeleton <- make_skeleton_CVM(U[1:t_train, ], trunc_tree = cfg$trunc_tree)
  cfg <- add_first_tree_map(cfg, skeleton)
  
  exports <- c(
    # constants & templates 
    "FAM_INFO", "FAM_INDEP", "FAM_GAUSS", "FAM_BB1", "FAM_BB1R180",
    "T_INDEP",  "T_GAUSS",   "T_BB1", "T_BB1R180", "FAM_BB8R180",
    "FAM_BB7","FAM_BB7R180","T_BB8R180","T_BB7","T_BB7R180",
    # helper functions 
    "active_fams", "sanitize_bb1", "mh_worker_standard", "mh_worker_block",
    "bb1r180_tail2par", "bb1r180_par2tail", "bb1r180_log_jacobian",
    "bb7_tail2par","bb7_par2tail","bb7_log_jacobian","bb7r180_tail2par","bb7r180_par2tail",
    "bb7r180_log_jacobian","bb8r180_tail2par","bb8r180_par2tail","bb8r180_log_jacobian_1d",
    "sanitize_bb7","sanitize_bb8",
    # core SMC kernels 
    "mh_step", "mh_step_in_tree",
    "update_weights", "ESS", "systematic_resample", "resample_move",
    # log-target & proposals 
    "log_prior", "bb1_tail2par", "bb1_par2tail", "bb1_log_jacobian",
    "rtnorm_vec", "log_prior_edge",
    # likelihood helpers 
    "bicop_dist", "vinecop_dist", "dvinecop", "fast_vine_from_particle",
    "rvinecop","bicop",
    # shared data objects 
    "skeleton", "cfg",
    # diagnostics & prediction 
    "diagnostic_report", "compute_predictive_metrics",
    "compute_log_incr",
    # small utilities    
    "w_mean", "w_var", "mc_se", "w_quantile", "fillna_neg", "fam_spec","get_tails","clamp01","init_from_tails",
    "tail_weights", "safe_logdens",
    "logit","ilogit","dlogitnorm",
    "emp_tails_FRAPO","seed_family_from_emp",
    "log_prior_edge_strong",".tip_means_for_edge_t","log_prior_with_tip_time","log_prior_with_tip_cached"
  )
  cl <- make_cluster(cfg$nc, cfg$seed, exports)
  parallel::clusterEvalQ(cl, { library(rvinecopulib); library(FRAPO) })
  
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers    <- colnames(U); A <- length(cfg$alphas)

  out <- list(
    log_pred    = numeric(n_oos),
    diag_log    = data.table(t = integer(N), ESS = numeric(N), unique = integer(N), euc = numeric(N), sparsity = numeric(N)),
    mh_acc_pct  = rep(NA_real_, N),
    step_sd_hist= rep(NA_real_, N),
    fam_hist  = array(NA_integer_, dim = c(M,N,K)),
    par1_hist = array(NA_real_,    dim = c(M,N,K)),
    par2_hist = array(NA_real_,    dim = c(M,N,K)),
    rotation_hist = array(NA_real_, dim = c(M,N,K)),
    ancestorIndices = matrix(0L, M, N),
    risk    =  list(
                      dates = integer(n_oos),                                        
                      VaR   = array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
                      ES    = array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
                    ),      
    
    QL    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    FZL   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    wCRPS = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    
    port = list(
                    dates  = integer(n_oos),
                    VaR    = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
                    ES     = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
                    QL    = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
                    FZL   = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
                    wCRPS = numeric(n_oos)
                  ),
    
    CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
  )

  
  # init particles
  U_init <- U[1:cfg$W, , drop = FALSE]
  particles <- replicate(M, new_particle(cfg, U_init = U_init), simplify = FALSE)
  out$ancestorIndices[,1] <- seq_len(M)
  

  for (t in 127:N) {
    
    #if (t==10){break}
    u_t <- U[t,,drop=FALSE]
    
    
    # weight-update
    log_inc <- compute_log_incr(particles, u_t, skeleton, cfg)
    particles <- update_weights(particles, log_inc)
    w <- vapply(particles, `[[`, numeric(1), "w")
    
    # prediction 
    if (t > cfg$W_predict) {
      # log predictive density, predictive draws, risk metrics
      idx <- t - cfg$W_predict
      
      y_real_t <- y_real[idx,]
      out$log_pred[idx] <- compute_predictive_metrics(u_t, particles, skeleton, w/sum(w), cfg)$log_pred_density
      
      #draws <- smc_predictive_sample(particles, skeleton, w/sum(w), L = 10000, cl = cl)
      
      
      draws <- smc_predictive_sample_fast(particles, skeleton, w/sum(w), L = 10000, cl = cl)
      
      
      Z_pred <- st_inv_fast(draws, shape_fc[idx, ], df_fc[idx, ])  
      R_t  <- sweep(Z_pred, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])          # L × d
      
      rs <- risk_stats_full(R_t, cfg$alphas)

      out$risk$dates[idx]   <- t
      out$risk$VaR [idx, , ] <- rs$VaR          # d × A
      out$risk$ES  [idx, , ] <- rs$ES
      
      # EW-portfolio metrics 
      r_p  <- rowMeans(R_t)                                           
      ps   <- port_stats(r_p, cfg$alphas)    
      
      out$port$dates[idx]   <- t                                       
      out$port$VaR [idx, ]  <- ps$VaR
      out$port$ES  [idx, ]  <- ps$ES 
      
      r_p_real <- mean(as.numeric(y_real_t))
      out$port$QL[idx, ]   <- vapply(seq_along(cfg$alphas), function(k) pinball_loss(r_p_real, ps$VaR[k], cfg$alphas[k]), numeric(1))
      out$port$FZL[idx, ]  <- vapply(seq_along(cfg$alphas), function(k) fzl_pzc_scalar(r_p_real, ps$VaR[k], ps$ES[k], cfg$alphas[k]), numeric(1))
      out$port$wCRPS[idx]  <- wcrps_gr_scalar(r_p, r_p_real)
      
      out$QL[idx, , ]  <- pinball_matrix(as.matrix(y_real_t[1, ]), rs$VaR, cfg$alphas) # Quantile loss per asset & alpha using the VaR you already computed
      out$FZL[idx, , ] <- fzl_pzc_matrix(as.matrix(y_real_t[1, ]), rs$VaR, rs$ES, cfg$alphas) # FZL joint loss for (VaR, ES)
      out$wCRPS[idx, ] <- wcrps_gr_matrix(R_t, as.matrix(y_real_t[1, ])) # Weighted CRPS from predictive draws 'R_t' and realization
      
      
      # CoVaR
      k5  <- which.min(abs(cfg$alphas - 0.05))
      k10 <- which.min(abs(cfg$alphas - 0.10))
      
      VaRj_5  <- rs$VaR[, k5]   # d-vector
      VaRj_10 <- rs$VaR[, k10]
      covar5  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.05, minN = 50)
      covar10 <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.10, minN = 50)
      
      out$CoVaR_tail[idx, , "a0.05"] <- covar5
      out$CoVaR_tail[idx, , "a0.10"] <- covar10
      
      
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
      newAncestors <- stratified_resample(w)
      data_up_to_t <- U[max(1, t - cfg$W + 1):t, , drop = FALSE]
      move_out <- resample_move_old(particles, newAncestors, data_up_to_t, cl,
                                  cfg$type, cfg, skeleton = skeleton)
      
      particles <- move_out$particles
      out$mh_acc_pct[t] <- move_out$acc_pct
      if (cfg$adapt_step_sd) {
        {cfg$step_sd <- compute_adapt_step_sd(cfg, move_out$acc_pct)}
        out$step_sd_hist[t] <- cfg$step_sd
      }
    } else {
      step_prev <- t - 1L
      newAncestors    <- if (step_prev < 1L) seq_len(M) else out$ancestorIndices[, step_prev]
    }
    out$ancestorIndices[, t] <- newAncestors
    
    # history arrays
    out$fam_hist [ , t,] <- t(vapply(particles, `[[`, integer(K),"fam"))
    out$par1_hist[, t,]  <- t(vapply(particles, `[[`, numeric(K),"th1"))
    out$par2_hist[, t,]  <- t(vapply(particles, `[[`, numeric(K),"th2"))
  }
  
  out$particles_final    <- particles
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  out
}






