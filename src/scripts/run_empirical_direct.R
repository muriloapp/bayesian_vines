SCEN_COVAR <- c(
  "a0.05b0.05", "a0.05b0.1", "a0.1b0.1", "a0.1b0.05",
  "a0.05b0.025", "a0.025b0.05"
)

# ============================================================
# run_empirical_OPEN_NO_SMCFULL.R
# - NO run_empirical()
# - NO smc_full() call
# - NO eval(body(smc_full))
# - Everything runs "open" in GlobalEnv
# - Runs cfg_variants sequentially and saves
# ============================================================

library(here)

# Load project into GLOBAL
source(here("src/R", "config.R"), local = .GlobalEnv)

# Data (once)
dat <- import_data(drop_first_col = TRUE, n_assets = 7)

# Variants
cfg_variants <- list(
  list(label = "tip_w252_M8000", use_tail_informed_prior = TRUE, W = 252L, M = 8000L),
  list(label = "tip_w126_M8000", use_tail_informed_prior = TRUE, W = 126L, M = 8000L),
  list(label = "tip_w252_M10000", use_tail_informed_prior = TRUE, W = 252L, M = 10000L),
  list(label = "tip_w126_M10000", use_tail_informed_prior = TRUE, W = 126L, M = 10000L)
  #list(label = "tip_w252_M3000_tip", use_tail_informed_prior = TRUE, W = 252L, M = 3000L),
  #list(label = "tip_w504_M3000_tip", use_tail_informed_prior = TRUE, W = 504L, M = 3000L)
  #list(label = "tip_w126_M3000_tip", use_tail_informed_prior = TRUE, W = 126L, M = 3000L)
)

out_dir <- here("empirical_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# OPEN RUNS (copy/paste the CONTENT of main_empirical.R / smc_full
# into the block below). Replace `data$...` with `dat$...` and use `cfg`.
#
# IMPORTANT:
# - This block MUST end by assigning a list named `res`
# - Then we save_result(res, out_dir)
# ============================================================

for (v in cfg_variants) {
  
  cat("\n==============================\n")
  cat("Running:", v$label, "\n")
  cat("==============================\n")
  
  # ---- cfg for this run ----
  cfg <- build_cfg(d = ncol(dat$U), trunc_tree = 4)
  cfg <- modifyList(cfg, v[setdiff(names(v), "label")])
  cfg$label <- v$label
  set.seed(cfg$seed)
  
  # ============================================================
  # OPEN EMPIRICAL RUN (paste your smc_full body here)
  # Use:
  #   U <- dat$U
  #   mu_fc <- dat$mu_fc[...]
  #   ...
  # and keep everything else as-is.
  # At the end, set:  res <- out
  # ============================================================
  
  # --- BEGIN OPEN BLOCK --------------------------------------
  U <- dat$U
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers <- colnames(U); A <- length(cfg$alphas); t_train <- cfg$W_predict
  
  mu_fc    <- dat$mu_fc[(nrow(dat$mu_fc) - n_oos + 1):nrow(dat$mu_fc), , drop = FALSE]
  sig_fc   <- dat$sig_fc[(nrow(dat$sig_fc) - n_oos + 1):nrow(dat$sig_fc), , drop = FALSE]
  df_fc    <- dat$df_fc[(nrow(dat$df_fc) - n_oos + 1):nrow(dat$df_fc), , drop = FALSE]
  shape_fc <- dat$shape_fc[(nrow(dat$shape_fc) - n_oos + 1):nrow(dat$shape_fc), , drop = FALSE]
  y_real   <- dat$y_real[(nrow(dat$y_real) - n_oos + 1):nrow(dat$y_real), , drop = FALSE]
  
  skeleton <- make_skeleton_CVM(U[1:t_train, ], trunc_tree = cfg$trunc_tree)
  cfg <- add_first_tree_map(cfg, skeleton)
  
  exports <- c(
    "FAM_INFO","FAM_INDEP","FAM_GAUSS","FAM_BB1","FAM_BB1R180",
    "T_INDEP","T_GAUSS","T_BB1","T_BB1R180","FAM_BB8R180",
    "FAM_BB7","FAM_BB7R180","T_BB8R180","T_BB7","T_BB7R180","FAM_T","T_T",
    "active_fams","sanitize_bb1",
    "bb1r180_tail2par","bb1r180_par2tail","bb1r180_log_jacobian",
    "bb7_tail2par","bb7_par2tail","bb7_log_jacobian","bb7r180_tail2par","bb7r180_par2tail",
    "bb7r180_log_jacobian","bb8r180_tail2par","bb8r180_par2tail","bb8r180_log_jacobian_1d",
    "sanitize_bb7","sanitize_bb8",
    "t_par2tail","t_tail2rho","t_log_jacobian","sanitize_t",
    "mh_step",
    "update_weights","ESS","systematic_resample",
    "log_prior","bb1_tail2par","bb1_par2tail","bb1_log_jacobian",
    "rtnorm_vec","log_prior_edge",
    "bicop_dist","vinecop_dist","dvinecop",
    "rvinecop","bicop",
    "diagnostic_report","compute_predictive_metrics", "compute_tip_sd_logit",
    "compute_log_incr",
    "w_mean","w_var","mc_se","w_quantile","fillna_neg","fam_spec","get_tails","clamp01","init_from_tails",
    "tail_weights","safe_logdens",
    "logit","ilogit","dlogitnorm",
    "emp_tails_FRAPO","seed_family_from_emp",
    "log_prior_edge_strong",".tip_means_for_edge_t","log_prior_with_tip_time","log_prior_with_tip_cached",
    ".safe_logdens1","fast_vine_from_row",".build_vine_from_vectors",
    ".as_particle_vectors","K_of_skeleton","safe_sample","r_lam"
  )
  
  cl <- make_cluster(cfg$nc, cfg$seed, exports, envir = .GlobalEnv)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
      VECLIB_MAXIMUM_THREADS="1", GOTO_NUM_THREADS="1"
    )
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    NULL
  })
  
  cl20 <- structure(cl[1:20], class = class(cl)) #*ADJUST*
  
  out <- list(
    log_pred = numeric(n_oos),
    pit_joint = numeric(n_oos),
    pit_rosen = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    diag_log = data.table(t = integer(N), ESS = numeric(N), unique = integer(N),
                          euc = numeric(N), sparsity = numeric(N)),
    mh_acc_pct = rep(NA_real_, N),
    step_sd_hist = rep(NA_real_, N),
    #fam_hist = array(NA_integer_, dim = c(M, N, K)),
    #par1_hist = array(NA_real_, dim = c(M, N, K)),
    #par2_hist = array(NA_real_, dim = c(M, N, K)),
    #rotation_hist = array(NA_real_, dim = c(M, N, K)),
    ancestorIndices = matrix(0L, M, N),
    risk = list(
      dates = integer(n_oos),
      VaR = array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
      ES  = array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
    ),
    QL = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    FZL = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    wCRPS = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    port = list(
      dates = integer(n_oos),
      VaR = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
      ES  = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
      QL  = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
      FZL = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
      wCRPS = numeric(n_oos)
    ),
    CoVaR_tail = array(NA_real_, c(n_oos, d, length(SCEN_COVAR)),
                       dimnames = list(NULL, tickers, SCEN_COVAR))
  )
  
  # init particles
  U_init <- U[1:(cfg$W - 1), , drop = FALSE]
  particles <- new_particles_mats(cfg, U_init)
  out$ancestorIndices[, 1] <- seq_len(M)
  
  for (t in (cfg$W + 1):N) {
    
    u_t_1 <- U[t - 1, , drop = FALSE]
    u_t   <- U[t, , drop = FALSE]
    
    log_inc <- compute_log_incr(particles, u_t_1, skeleton)
    particles <- update_weights(particles, log_inc)
    w <- particles$w
    
    if (t > cfg$W_predict) {
      idx <- t - cfg$W_predict
      
      y_real_t <- y_real[idx, ]
      out$log_pred[idx] <- compute_predictive_metrics(u_t, particles, skeleton, w / sum(w), cfg)$log_pred_density
      
      draws <- smc_predictive_sample_fast2_scoped2(particles, skeleton, w, L = 30000, cl = cl20)
      
      cmp  <- sweep(draws, 2, as.numeric(u_t), FUN = "<=")
      pitV <- matrixStats::rowAlls(cmp)
      out$pit_joint[idx] <- mean(pitV)
      
      out$pit_rosen[idx, ] <- empirical_rosenblatt_from_draws(as.numeric(u_t), draws, K = floor(sqrt(nrow(draws))))
      
      Z_pred <- st_inv_fast(draws, shape_fc[idx, ], df_fc[idx, ])
      R_t <- sweep(Z_pred, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
      
      rs <- risk_stats_full(R_t, cfg$alphas)
      out$risk$dates[idx] <- t
      out$risk$VaR[idx, , ] <- rs$VaR
      out$risk$ES[idx, , ]  <- rs$ES
      
      r_p <- rowMeans(R_t)
      ps <- port_stats(r_p, cfg$alphas)
      out$port$dates[idx] <- t
      out$port$VaR[idx, ] <- ps$VaR
      out$port$ES[idx, ]  <- ps$ES
      
      r_p_real <- mean(as.numeric(y_real_t))
      out$port$QL[idx, ]  <- vapply(seq_along(cfg$alphas), function(k) pinball_loss(r_p_real, ps$VaR[k], cfg$alphas[k]), numeric(1))
      out$port$FZL[idx, ] <- vapply(seq_along(cfg$alphas), function(k) fzl_pzc_scalar(r_p_real, ps$VaR[k], ps$ES[k], cfg$alphas[k]), numeric(1))
      out$port$wCRPS[idx] <- wcrps_gr_scalar(r_p, r_p_real)
      
      out$QL[idx, , ]  <- pinball_matrix(t(as.matrix(y_real_t)), rs$VaR, cfg$alphas)
      out$FZL[idx, , ] <- fzl_pzc_matrix(t(as.matrix(y_real_t)), rs$VaR, rs$ES, cfg$alphas)
      out$wCRPS[idx, ] <- wcrps_gr_matrix(R_t, t(as.matrix(y_real_t)))
      
      k025 <- which.min(abs(cfg$alphas - 0.025))
      k5   <- which.min(abs(cfg$alphas - 0.05))
      k10  <- which.min(abs(cfg$alphas - 0.10))
      
      VaRj_025 <- rs$VaR[, k025]
      VaRj_5   <- rs$VaR[, k5]
      VaRj_10  <- rs$VaR[, k10]
      
      out$CoVaR_tail[idx, , "a0.05b0.05"]   <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.05,  minN = 50)
      out$CoVaR_tail[idx, , "a0.05b0.1"]    <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.1,   minN = 50)
      out$CoVaR_tail[idx, , "a0.1b0.1"]     <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.10,  minN = 50)
      out$CoVaR_tail[idx, , "a0.1b0.05"]    <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.05,  minN = 50)
      out$CoVaR_tail[idx, , "a0.05b0.025"]  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.025, minN = 50)
      out$CoVaR_tail[idx, , "a0.025b0.05"]  <- covar_tail_vec(R_t, r_p, VaRj_025,port_alpha = 0.05,  minN = 50)
    }
    
    if (t %% 250 == 0L){
    dg <- diagnostic_report(t, 0, U, particles, w, cfg)
    out$diag_log[t, `:=`(t = t, ESS = dg$ESS, unique = dg$unique, euc = dg$euc, sparsity = dg$sparsity)]
    }
    if (ESS(w) < cfg$ess_thr * M && t < N) {
      newAncestors <- stratified_resample(w)
      data_up_to_t <- U[max(1, t - cfg$W + 1):(t - 1), , drop = FALSE]
      
      move_out <- resample_move_old(particles, newAncestors, data_up_to_t, cl,
                                    cfg$type, cfg, skeleton = skeleton)
      
      particles <- move_out$particles
      out$mh_acc_pct[t] <- move_out$acc_pct
      if (cfg$adapt_step_sd) {
        cfg$step_sd <- compute_adapt_step_sd(cfg, move_out$acc_pct)
        out$step_sd_hist[t] <- cfg$step_sd
      }
    } else {
      step_prev <- t - 1L
      newAncestors <- if (step_prev < 1L) seq_len(M) else out$ancestorIndices[, step_prev]
    }
    
    out$ancestorIndices[, t] <- newAncestors
    #out$fam_hist[, t, ] <- particles$fam_mat
    #out$par1_hist[, t, ] <- particles$th1_mat
    #out$par2_hist[, t, ] <- particles$th2_mat
  }
  
  out$cfg <- cfg
  out$particles_final <- particles
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  
  res <- out
  # --- END OPEN BLOCK ----------------------------------------
  
  # Save
  res$cfg <- cfg
  save_result(res, out_dir)
  
  # Cleanup between variants
  rm(res, out, particles, skeleton, cl, cl20)
  gc()
}

cat("\nDONE.\n")
