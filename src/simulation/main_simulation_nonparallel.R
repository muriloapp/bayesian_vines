
ci_quantile_boot <- function(x, a, B = 500, level = 0.95) {
  stopifnot(level > 0 && level < 1, a > 0 && a < 1, B >= 1)
  
  # helper for one vector
  one <- function(v) {
    v <- as.numeric(v)
    L <- length(v)
    vals <- replicate(B, {
      vb <- sample(v, L, replace = TRUE)
      as.numeric(quantile(vb, probs = a, type = 8))
    })
    as.numeric(quantile(vals,
                        probs = c((1 - level)/2, 1 - (1 - level)/2),
                        names = FALSE))
  }
  
  if (is.matrix(x)) {
    # return 2 x ncol matrix: rows = (lower, upper), cols = original columns
    out <- sapply(seq_len(ncol(x)), function(j) one(x[, j]))
    rownames(out) <- c("lower", "upper")
    colnames(out) <- colnames(x) %||% paste0("col", seq_len(ncol(x)))
    return(out)
  } else {
    return(one(x))
  }
}


ci_covar_boot_both <- function(R, alpha = 0.05, beta = 0.05,
                               B = 1000, level = 0.95, type = 8,
                               min_tail = 50) {
  stopifnot(is.matrix(R), ncol(R) == 2)
  L <- nrow(R)
  qs <- c((1 - level)/2, 1 - (1 - level)/2)
  
  covar_2_given_1 <- function(Rb) {
    thr <- as.numeric(quantile(Rb[, 1], probs = beta, type = type))
    sel <- (Rb[, 1] <= thr)
    if (sum(sel) < min_tail) return(NA_real_)
    as.numeric(quantile(Rb[sel, 2], probs = alpha, type = type))
  }
  
  covar_1_given_2 <- function(Rb) {
    thr <- as.numeric(quantile(Rb[, 2], probs = beta, type = type))
    sel <- (Rb[, 2] <= thr)
    if (sum(sel) < min_tail) return(NA_real_)
    as.numeric(quantile(Rb[sel, 1], probs = alpha, type = type))
  }
  
  vals21 <- numeric(B)
  vals12 <- numeric(B)
  
  for (b in seq_len(B)) {
    ii <- sample.int(L, L, replace = TRUE)  # resample ROWS (keeps dependence)
    Rb <- R[ii, , drop = FALSE]
    vals21[b] <- covar_2_given_1(Rb)
    vals12[b] <- covar_1_given_2(Rb)
  }
  
  vals21 <- vals21[is.finite(vals21)]
  vals12 <- vals12[is.finite(vals12)]
  
  ci21 <- if (length(vals21) >= max(30, 0.1 * B)) {
    setNames(as.numeric(quantile(vals21, probs = qs, names = FALSE)), c("lower", "upper"))
  } else c(lower = NA_real_, upper = NA_real_)
  
  ci12 <- if (length(vals12) >= max(30, 0.1 * B)) {
    setNames(as.numeric(quantile(vals12, probs = qs, names = FALSE)), c("lower", "upper"))
  } else c(lower = NA_real_, upper = NA_real_)
  
  list(
    CoVaR_2_given_1 = ci21,
    CoVaR_1_given_2 = ci12
  )
}


ci_covar_boot_both_moutn <- function(R, alpha=0.05, beta=0.05, B=100, level=0.95, type=8,
                                     min_tail=50, m_frac=0.5) {
  L <- nrow(R); m <- max(200L, floor(m_frac * L))
  qs <- c((1-level)/2, 1-(1-level)/2)
  
  one <- function(Rb, cond_col, target_col) {
    thr <- as.numeric(quantile(Rb[,cond_col], beta, type=type))
    sel <- (Rb[,cond_col] <= thr)
    if (sum(sel) < min_tail) return(NA_real_)
    as.numeric(quantile(Rb[sel, target_col], alpha, type=type))
  }
  
  v21 <- v12 <- numeric(B)
  for (b in 1:B) {
    ii <- sample.int(L, m, replace=TRUE)
    Rb <- R[ii, , drop=FALSE]
    v21[b] <- one(Rb, 1, 2)
    v12[b] <- one(Rb, 2, 1)
  }
  v21 <- v21[is.finite(v21)]; v12 <- v12[is.finite(v12)]
  
  list(
    CoVaR_2_given_1 = setNames(as.numeric(quantile(v21, qs, names=FALSE)), c("lower","upper")),
    CoVaR_1_given_2 = setNames(as.numeric(quantile(v12, qs, names=FALSE)), c("lower","upper"))
  )
}




smc_simul_serial <- function(data, cfg, dgp) {
  
  true_bases <- dgp$true_bases
  U      <- data$U
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers    <- colnames(U); A <- length(cfg$alphas); t_train <- cfg$W_predict
  
  mu_fc   <- data$mu_fc[(nrow(data$mu_fc)   - n_oos + 1):nrow(data$mu_fc), , drop = FALSE]
  sig_fc  <- data$sig_fc[(nrow(data$sig_fc) - n_oos + 1):nrow(data$sig_fc), , drop = FALSE]
  df_fc   <- data$df_fc[(nrow(data$df_fc)   - n_oos + 1):nrow(data$df_fc), , drop = FALSE]
  shape_fc<- data$shape_fc[(nrow(data$shape_fc) - n_oos + 1):nrow(data$shape_fc), , drop = FALSE]
  y_real  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  
  
  skeleton <- make_skeleton_CVM(U[1:t_train, ], trunc_tree = cfg$trunc_tree, structure = dgp$vc$structure)
  cfg <- add_first_tree_map(cfg, skeleton)
  # 
  # exports <- c(
  #   # constants & templates 
  #   "FAM_INFO", "FAM_INDEP", "FAM_GAUSS", "FAM_BB1", "FAM_BB1R180",
  #   "T_INDEP",  "T_GAUSS",   "T_BB1", "T_BB1R180", "FAM_BB8R180",
  #   "FAM_BB7","FAM_BB7R180","T_BB8R180","T_BB7","T_BB7R180", "FAM_T","T_T",
  #   # helper functions 
  #   "active_fams", "sanitize_bb1", "mh_worker_standard", "mh_worker_block",
  #   "bb1r180_tail2par", "bb1r180_par2tail", "bb1r180_log_jacobian",
  #   "bb7_tail2par","bb7_par2tail","bb7_log_jacobian","bb7r180_tail2par","bb7r180_par2tail",
  #   "bb7r180_log_jacobian","bb8r180_tail2par","bb8r180_par2tail","bb8r180_log_jacobian_1d",
  #   "sanitize_bb7","sanitize_bb8",
  #   "t_par2tail","t_tail2rho","t_log_jacobian","sanitize_t",
  #   # core SMC kernels 
  #   "mh_step", "mh_step_in_tree",
  #   "update_weights", "ESS", "systematic_resample",
  #   # log-target & proposals 
  #   "log_prior", "bb1_tail2par", "bb1_par2tail", "bb1_log_jacobian",
  #   "rtnorm_vec", "log_prior_edge",
  #   # likelihood helpers 
  #   "bicop_dist", "vinecop_dist", "dvinecop",
  #   "rvinecop","bicop",
  #   # shared data objects 
  #   "skeleton", "cfg",
  #   # diagnostics & prediction 
  #   "diagnostic_report", "compute_predictive_metrics",
  #   "compute_log_incr",
  #   # small utilities    
  #   "w_mean", "w_var", "mc_se", "w_quantile", "fillna_neg", "fam_spec","get_tails","clamp01","init_from_tails",
  #   "tail_weights", "safe_logdens",
  #   "logit","ilogit","dlogitnorm",
  #   "emp_tails_FRAPO","seed_family_from_emp",
  #   "log_prior_edge_strong",".tip_means_for_edge_t","log_prior_with_tip_time","log_prior_with_tip_cached",
  #   ".safe_logdens1","fast_vine_from_row",".build_vine_from_vectors",
  #   ".as_particle_vectors","K_of_skeleton","safe_sample",
  #   "true_bases", "dgp"
  # )


  out <- list(
    log_pred    = numeric(n_oos),
    pit_joint = numeric(n_oos),
    pit_rosen = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
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
    
    CoVaR_tail = array(
      NA_real_, c(n_oos, d, length(SCEN_COVAR)),
      dimnames = list(NULL, tickers, SCEN_COVAR)
    )
    )
  
  
  # init particles
  U_init <- U[1:(cfg$W-1), , drop = FALSE]
  particles <- new_particles_mats(cfg, U_init, true_bases = true_bases) # TRUE BASES
  out$ancestorIndices[,1] <- seq_len(M)
  
  

for (t in (cfg$W+1):N) {
    u_t_1 <- U[t-1,,drop=FALSE]
    u_t <- U[t,,drop=FALSE]

    log_inc <- compute_log_incr(particles, u_t_1, skeleton)
    particles <- update_weights(particles, log_inc)
    w <- particles$w
    
    if (t > cfg$W_predict) {
      idx <- t - cfg$W_predict
      
      y_real_t <- y_real[idx,]
      out$log_pred[idx] <- compute_predictive_metrics(u_t, particles, skeleton, w/sum(w), cfg)$log_pred_density

      #system.time(
      #draws <- smc_predictive_sample_fast2_scoped(particles, skeleton, w, L = 20000, cl = cl[1:8])
      #)
      # system.time(
      #   R_draws <- smc_predictive_sample_fast2_grouped_epoch(
      #     w   = particles$w,
      #     L   = 20000,         # or whatever you need
      #     cl  = cl,
      #     round_digits = 6,    # more aggressive (e.g., 5) groups more; try 5–6
      #     top_k = NULL,        # optional
      #     min_count = 0,       # optional
      #     nc = 1               # keep 1 unless W*nc <= physical cores
      #   )
      # )
      #system.time(
      draws <- smc_predictive_sample_fast2_scoped2_serial(particles, skeleton, w, L = 10000)
      #)
      #cmp  <- sweep(draws, 2, as.numeric(u_t), FUN = "<=")  # L x d logical
      #pitV <- matrixStats::rowAlls(cmp)   
      #out$pit_joint[idx] <- mean(pitV)   
      
      #out$pit_rosen[idx, ] <- empirical_rosenblatt_from_draws(as.numeric(u_t), draws, K = floor(sqrt(nrow(draws))))
      

      Z_pred <- st_inv_fast(draws, shape_fc[idx, ], df_fc[idx, ])  
      R_t  <- sweep(Z_pred, 2, as.numeric(sqrt(sig_fc[idx, ])), `*`) + as.numeric(mu_fc[idx, ])          # L × d
  
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
      
      out$QL[idx, , ]  <- pinball_matrix(t(as.matrix(y_real_t)), rs$VaR, cfg$alphas) # Quantile loss per asset & alpha using the VaR you already computed
      out$FZL[idx, , ] <- fzl_pzc_matrix(t(as.matrix(y_real_t)), rs$VaR, rs$ES, cfg$alphas) # FZL joint loss for (VaR, ES)
      out$wCRPS[idx, ] <- wcrps_gr_matrix(R_t, t(as.matrix(y_real_t))) # Weighted CRPS from predictive draws 'R_t' and realization
      
      
      # CoVaR
      k5  <- which.min(abs(cfg$alphas - 0.05))
      k10 <- which.min(abs(cfg$alphas - 0.10))
      k025 <- which.min(abs(cfg$alphas - 0.025))
      
      VaRj_025 <- rs$VaR[, k025]
      VaRj_5  <- rs$VaR[, k5]   # d-vector
      VaRj_10 <- rs$VaR[, k10]
      
      
      covar5  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.05, minN = 50)
      covar5b10  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.1, minN = 50)
      covar10 <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.10, minN = 50)
      covar10b5 <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.05, minN = 50)
      covar5b0025 <- covar_tail_vec(R_t, r_p, VaRj_5,   port_alpha = 0.025, minN = 50)
      covar025b5  <- covar_tail_vec(R_t, r_p, VaRj_025, port_alpha = 0.05,  minN = 50)
      
      
      out$CoVaR_tail[idx, , "a0.05b0.05"] <- covar5
      out$CoVaR_tail[idx, , "a0.05b0.1"] <- covar5b10
      out$CoVaR_tail[idx, , "a0.1b0.1"] <- covar10
      out$CoVaR_tail[idx, , "a0.1b0.05"] <- covar10b5
      out$CoVaR_tail[idx, , "a0.05b0.025"]  <- covar5b0025
      out$CoVaR_tail[idx, , "a0.025b0.05"]  <- covar025b5
      
      
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
      data_up_to_t <- U[max(1, t - cfg$W + 1):(t-1), , drop = FALSE]
      
      move_out <- resample_move_old_serial(particles, newAncestors, data_up_to_t, NULL,
                                    cfg$type, cfg, skeleton = skeleton, true_bases=true_bases) #true_bases
      
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
    
    # after MH/resample at time t:
    #push_epoch_state(cl, particles, skeleton, epoch_id = t)
    out$ancestorIndices[, t] <- newAncestors
    
    # history arrays
    out$fam_hist [ , t,] <- particles$fam_mat
    out$par1_hist[, t,]  <- particles$th1_mat
    out$par2_hist[, t,]  <- particles$th2_mat

}
  
  out$cfg <- cfg
  out$particles_final    <- particles
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  out
}










