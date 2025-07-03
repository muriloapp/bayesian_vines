
# Block vine Block_vine
# refactor trees and adjust mhstep

# ================================================================
#  SMC para C-vine com spike-and-slab 
# ================================================================


assignInNamespace("assert_that", function(...) invisible(TRUE), ns = "assertthat")
assignInNamespace("see_if", function(...) invisible(TRUE), ns = "assertthat")

source(here('src','core_functions.R'))
source(here('src','simulation.R'))


#Rcpp::sourceCpp(here::here("src", "calc_loglik_tree_tr.cpp"))


run_block_smc <- function(data,
                          cfg,
                          type       = "block",
                          n_cores    = max(parallel::detectCores() - 1, 1)) {
  
  U <- data$U
  
  N  <- nrow(U);  d <- ncol(U)
  M  <- cfg$M;    K <- cfg$K;     G <- cfg$G
  S  <- N * G                         #  total sub-steps

  full_skeleton <- vinecop(U, family_set = "gaussian", structure = data$RVM$Matrix[nrow(data$RVM$Matrix):1, ])        
  skeletons_by_tr <- lapply(seq_len(d - 1L), function(tr) {
    vinecop(U,
            family_set = "gaussian",
            structure  = full_skeleton$structure,
            trunc_lvl  = tr)
  })
  
  particles <- replicate(M, new_particle(cfg), simplify = FALSE)
  
  n_pairs <- N * cfg$G
  out <- list(
    log_pred   = numeric(N),
    theta_mean = matrix(NA_real_, N, K),
    theta_se   = matrix(NA_real_, N, K),
    gamma_mean = matrix(NA_real_, N, K),
    gamma_se   = matrix(NA_real_, N, K),
    diag_log   = data.table::data.table(
      t      = integer(n_pairs),
      tr     = integer(n_pairs),
      ESS    = numeric(n_pairs),
      unique = integer(n_pairs),
      euc    = numeric(n_pairs),
      pi_mean = numeric(n_pairs),      # NEW
      pi_sd   = numeric(n_pairs)       # NEW
    ),
    mh_acc_pct      = rep(NA_real_, n_pairs),
    step_sd_hist      = rep(NA_real_, N),
    theta_hist      = array(NA_real_,    dim = c(M, S, K)),
    #gamma_hist      = array(NA_integer_, dim = c(M, S, K)),
    ancestorIndices = matrix(0L, M, S),
    incl_hist = matrix(NA_real_, S, K),
    theta_q025   = matrix(NA_real_, N, K),
    theta_q975   = matrix(NA_real_, N, K)
  )
  out$ancestorIndices[, 1L] <- seq_len(M)
  
  # ── 3 • parallel backend ─────────────────────────────────────────────────
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, cfg$seed)
  
  parallel::clusterExport(
    cl, c("mh_step_in_tree", "vine_from_particle", "log_prior", "slab_sd_from_tau", "spike_sd_from_tau", "update_tau2", "rinvgamma", "dinvgamma", "update_pi",
          "bicop_dist", "vinecop_dist", "dvinecop", "full_skeleton", "fast_vine_from_particle",
          "skeletons_by_tr", "cfg", "ESS", "propagate_particles",
          "update_weights", "diagnostic_report", "systematic_resample",
          "resample_move", "calculate_log_lik_tree_tr",
          "compute_predictive_metrics"),
    envir = environment()
  )
  
  step_id <- 0L    # flat counter for sub-steps
  for (t_idx in seq_len(N)) {
    
    u_row     <- U[t_idx, , drop = FALSE]
    #particles <- propagate_particles(particles, cfg)   # state advance
    
    # ---- Predictive metrics once per row ----------------------------------
    if (t_idx > cfg$W_predict) {
      w_prev <- vapply(particles, `[[`, numeric(1), "w")
      w_prev <- w_prev / sum(w_prev)
      
      pm <- compute_predictive_metrics(u_row, particles,
                                       full_skeleton, w_prev, cfg)
      
      out$log_pred[t_idx]    <- pm$log_pred_density
      out$theta_mean[t_idx,] <- pm$theta_mean
      out$theta_se[t_idx,]   <- pm$theta_se
      out$gamma_mean[t_idx,] <- pm$gamma_mean
      out$gamma_se[t_idx,]   <- pm$gamma_se
    }
    
    # ---- Inner loop over G “tree moves” -----------------------------------
    for (tr_idx in seq_len(G)) { #CHECK THIS
      step_id <- step_id + 1L
      
      tmp_skel <- if (tr_idx == G) {  # final move uses deepest skeleton
        skeletons_by_tr[[length(skeletons_by_tr)]]
      } else {
        skeletons_by_tr[[tr_idx]]
      } 
      
      tr_prev <- tr_idx - 1L
      log_incr <- vapply(
        particles,
        function(p) calculate_log_lik_tree_tr(
          p, tmp_skel, u_row, t_idx, tr_idx, tr_prev,
          skeletons_by_tr, cfg),
        numeric(1)
      )
      particles <- update_weights(particles, log_incr)
      w_new     <- vapply(particles, `[[`, numeric(1), "w")
      
      # Diagnostics ---------------------------------------------------------
      dg <- diagnostic_report(t_idx, tr_idx, U, particles, w_new, cfg)
      out$diag_log[step_id, `:=`(
        t      = t_idx,
        tr     = tr_idx,
        ESS    = dg$ESS,
        unique = dg$unique,
        euc    = dg$euc,
        tau_mean = dg$tau_mean,        
        tau_sd   = dg$tau_sd,           
        pi_mean = dg$pi_mean,
        pi_sd = dg$pi_sd 
      )]
      
      # Resample + block-MH move -------------------------------------------
      if (ESS(w_new) < cfg$ess_thr * M && t_idx < N) {
        
        data_up_to_t <- U[max(1, t_idx - cfg$W + 1) : t_idx, , drop = FALSE]
        newAnc <- stratified_resample(w_new)
        move_out    <- resample_move(particles, newAnc, data_up_to_t,
                                      cl, type, cfg, tr=tr_idx, temp_skel = tmp_skel)
        particles <- move_out$particles
        out$mh_acc_pct[step_id] <- move_out$acc_pct
        
        if (cfg$adapt_step_sd) {
          {cfg$step_sd <- compute_adapt_step_sd(cfg, move_out$acc_pct)}
          out$step_sd_hist[t_idx] <- cfg$step_sd
        }
      } else {
        prev_step <- step_id - 1L
        newAnc    <- if (prev_step < 1L) seq_len(M) else out$ancestorIndices[, prev_step]
      }
      
      # Genealogy + history -------------------------------------------------
      out$ancestorIndices[, step_id] <- newAnc
      out$theta_hist[ , step_id, ] <- t(vapply(particles, `[[`, numeric(K), "theta"))
      #out$gamma_hist[ , step_id, ] <- t(vapply(particles, `[[`, integer(K), "gamma"))
      
      theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
      tau_vec   <- sqrt(vapply(particles, function(p) p$tau2, numeric(1)))
      
      pi_vec  <- vapply(particles, function(p) p$pi, numeric(1))
      slab_w  <- responsibility(theta_mat, tau_vec, pi_vec, cfg)  
      
      out$incl_hist[step_id, ] <- colSums(slab_w * w_new) 
      out$theta_q025[step_id, ] <- dg$edges$q025[[1]]
      out$theta_q975[step_id, ] <- dg$edges$q975[[1]]
    }  # end tr
  } # end t 
  
  out$log_model_evidence <- sum(out$log_pred, na.rm = TRUE)
  out$particles_final    <- particles
  return(out)
}



