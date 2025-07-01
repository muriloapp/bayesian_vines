library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(Rcpp)
library(matrixStats)


spike_sd_from_tau <- function(tau){tau}

slab_sd_from_tau  <- function(tau, cfg){cfg$c_slab * tau}

log_prior <- function(p,cfg){
  tau <- sqrt(p$tau2)
  pi  <- p$pi                             
  
  ll_mix <- sum(
    log( pi      * dnorm(p$theta,0, spike_sd_from_tau(tau)) +
           (1-pi)  * dnorm(p$theta,0, slab_sd_from_tau(tau,cfg)) + 1e-300))
  
  ll_tau <- if(cfg$tau_prior=="fixed") 0
  else dinvgamma(p$tau2, cfg$a0, cfg$b0, log=TRUE)
  
  ll_pi  <- if(cfg$pi_prior=="fixed") 0
  else dbeta(pi, cfg$a_pi, cfg$b_pi, log=TRUE)
  
  ll_mix + ll_tau + ll_pi
}

new_particle <- function(cfg){
  # --- τ² block unchanged -----------------------------------
  tau2 <- if(cfg$tau_prior=="fixed") {cfg$tau0^2}
  else {rinvgamma(1, cfg$a0, cfg$b0)}
  tau  <- sqrt(tau2)
  
  # --- NEW π draw ------------------------------------------
  pi   <- if(cfg$pi_prior=="fixed") cfg$pi0
  else rbeta(1, cfg$a_pi, cfg$b_pi)
  
  # --- θ draws from *mixture* -------------------------------
  is_slab <- rbinom(cfg$K, 1, 1 - pi)
  
  aux_slab <- slab_sd_from_tau(tau,cfg)
  aux_spike <- spike_sd_from_tau(tau)
  theta   <- rnorm(cfg$K, 0,
                   ifelse(is_slab,
                          aux_slab,
                          aux_spike))
  
  list(theta = theta,
       tau2  = tau2,
       pi    = pi,        # << store π
       w     = 1/cfg$M,
       last_accept = FALSE,
       aux_slab = aux_slab,
       aux_spike = aux_spike,
       is_slab = is_slab
       )
}


vine_from_particle <- function(p, skel, cfg) {
  with(cfg, {
    pcs <- lapply(skel$pair_copulas, function(lvl) lapply(lvl, identity))
    idx <- 1L
    for (tr in seq_along(pcs)) {
      for (ed in seq_along(pcs[[tr]])) {
        pcs[[tr]][[ed]] <- bicop_dist("gaussian", parameters = tanh(p$theta[idx]))
        idx <- idx + 1L
      }
    }
    vinecop_dist(pcs, structure = skel$structure)
  })
}


update_pi <- function(p,cfg){
  if(cfg$pi_prior=="fixed") return(p)
  
  tau <- sqrt(p$tau2); c2 <- cfg$c_slab^2
  # responsibilities for 'slab'
  den_sp <- p$pi      * dnorm(p$theta,0, tau)
  den_sl <- (1-p$pi)  * dnorm(p$theta,0, sqrt(c2)*tau)
  w_sl   <- den_sl / (den_sp + den_sl + 1e-300)
  
  m1 <- sum(w_sl)           # expected slab count
  m0 <- cfg$K - m1
  
  p$pi <- rbeta(1, cfg$a_pi + m0, cfg$b_pi + m1)
  p
}

rinvgamma <- function(n, shape, scale)  1/rgamma(n, shape, rate = scale)

dinvgamma <- function(x, shape, scale, log = FALSE){
  lg <- shape*log(scale) - lgamma(shape) -
    (shape+1)*log(x) - scale/x
  if(log) lg else exp(lg)
}

update_tau2 <- function(p, cfg){
  if (cfg$tau_prior == "fixed") return(p)     # no-op for fixed model
  
  tau  <- sqrt(p$tau2)
  c2   <- cfg$c_slab^2
  v_sp <- tau^2
  v_sl <- c2 * v_sp
  
  ## 1. posterior slab responsibility for each θ
  den_sp <- p$pi      * dnorm(p$theta,0, sqrt(v_sp))
  den_sl <- (1-p$pi)  * dnorm(p$theta,0, sqrt(v_sl))
  w_sl   <- den_sl/(den_sp+den_sl+1e-300)
  
  ## 2. sum of squares with correct scaling
  s0 <- sum((1 - w_sl) * p$theta^2)        # spike part
  s1 <- sum(w_sl      * p$theta^2) / c2    # slab part
  
  ## 3. IG update
  a_post <- cfg$a0 + cfg$K/2
  b_post <- cfg$b0 + 0.5 * (s0 + s1)
  
  p$tau2 <- rinvgamma(1, a_post, b_post)
  p
}


ESS <- function(w) 1 / sum(w^2)

mh_step_in_tree <- function(p, tr, data_up_to_t, temp_skel, cfg) {
  
  ## ──────────────────────────────
  ## 1.  Identify the edges in tree tr
  ## ──────────────────────────────
  if (tr == 1) {
    edges_in_previous_trees <- 0L
  } else {
    edges_in_previous_trees <- cfg$d * (tr - 1) - (tr - 1) * tr / 2
  }
  first_edge <- edges_in_previous_trees + 1L
  last_edge  <- edges_in_previous_trees + (cfg$d - tr)
  if (tr == cfg$G) {
    last_edge <- cfg$K
  }
  idx_tree   <- first_edge:last_edge        # vector of edge indices
  
  ## ──────────────────────────────
  ## 2.  Propose γ / θ on those edges
  ## ──────────────────────────────
  prop <- p
  
  ## 2b.  Random-walk on active θ

  prop$theta <- p$theta + rnorm(cfg$K, 0, cfg$step_sd)
  
  ## ──────────────────────────────
  ## 3.  Log-likelihood + prior
  ## ──────────────────────────────
  vine_prop <- vine_from_particle(prop, temp_skel, cfg)
  prop_ll   <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
                            .Machine$double.eps)))
  
  vine_curr <- vine_from_particle(p, temp_skel, cfg)
  curr_ll   <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
                            .Machine$double.eps)))
  
  ## ──────────────────────────────
  ## 5.  MH acceptance
  ## ──────────────────────────────
  log_acc <- (prop_ll + log_prior(prop, cfg)) -
    (curr_ll + log_prior(p, cfg))  
  
  if (log(runif(1)) < log_acc) {
    p$theta       <- prop$theta
    p$last_accept <- TRUE
  } else {
    p$last_accept <- FALSE
  }
  ## NEW ---------------------------------------------------------------
  p <- update_tau2(p, cfg)
  p <- update_pi(p,cfg)        # << NEW
  return(p)
}


mh_step <- function(p, data_up_to_t, skeleton, cfg) {
  prop <- p
  prop$theta <- p$theta + rnorm(cfg$K, 0, cfg$step_sd)
  
  vine_prop <- vine_from_particle(prop, skeleton, cfg)
  vine_curr <- vine_from_particle(p,    skeleton, cfg)
  
  ll_prop <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
                          .Machine$double.eps)))
  ll_curr <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
                          .Machine$double.eps)))
  
  # log_acc <- (ll_prop + log_prior(prop$theta, cfg)) -
  #   (ll_curr + log_prior(p$theta, cfg))
  log_acc <- (ll_prop + log_prior(prop, cfg)) -      # ← pass WHOLE prop
    (ll_curr + log_prior(p,    cfg))       # ← pass WHOLE p
  
  if (log(runif(1)) < log_acc) {
    p$theta <- prop$theta
    p$last_accept <- TRUE
  } else p$last_accept <- FALSE
  
  ## NEW ---------------------------------------------------------------
  p <- update_tau2(p, cfg)
  p <- update_pi(p,cfg)        # << NEW
  p
}


# # Adaptive Proposal SD
# adaptive_proc_sd <- function(theta_mat, w_new, scale = 0.3, floor_val = 1e-3) {
#   # Posterior SD per component of theta
#   post_sd <- apply(theta_mat, 2, function(col)
#     sqrt(w_var(col, w_new)))
#   
#   # Adaptively scaled proposal SD
#   proc_sd <- scale * pmax(post_sd, floor_val)
#   
#   return(proc_sd)
# }

calculate_log_lik_tree_tr <- function(particle, skel_tr, u_row, t, tr, tr_prev, skeletons_by_tr, cfg) {
  # If it's the first tree, calculate its log-density
  if (tr == 1) {
    vine_j <- vine_from_particle(particle, skel_tr, cfg)
    dens_val <- dvinecop(u_row, vine_j)
    return(log(dens_val + 1e-100))
  }

  # For subsequent trees, we need the density of the current and previous model
  # The pre-computation of skeletons_by_tr makes this much faster.
  skel_tr_minus_1 <- skeletons_by_tr[[tr_prev]]

  vine_j_tr <- vine_from_particle(particle, skel_tr, cfg)
  dens_val_tr <- dvinecop(u_row, vine_j_tr)

  vine_j_prev <- vine_from_particle(particle, skel_tr_minus_1, cfg)
  dens_val_prev <- dvinecop(u_row, vine_j_prev)

  # Return the log-difference
  return(log(dens_val_tr + 1e-100) - log(dens_val_prev + 1e-100))
}

propagate_particles <- function(particles, cfg) {
  lapply(particles, function(p) {
    # Random‑walk on active thetas
    active <- which(p$gamma == 1L)
    if (length(active))
      p$theta[active] <- p$theta[active] + rnorm(length(active), 0, cfg$proc_sd)
    
    # Potential gamma flips
    flip <- which(runif(cfg$K) < cfg$p_dyn_flip)
    if (length(flip)) {
      for (e in flip) {
        if (p$gamma[e]) { p$gamma[e] <- 0L; p$theta[e] <- 0 }
        #else { p$gamma[e] <- 1L; p$theta[e] <- rnorm(1, 0, cfg$slab_sd) }
        else {
          tau          <- sqrt(p$tau2)
          p$gamma[e]   <- 1L
          p$theta[e]   <- rnorm(1, 0, slab_sd_from_tau(tau, cfg))
        }
      }
    }
    p
  })
}

compute_log_incr <- function(particles, u_row, skeleton, cfg) {
  log_incr <- numeric(cfg$M)
  for (j in seq_along(particles)) {
      p <- particles[[j]]
      vine_j     <- vine_from_particle(p, skeleton, cfg)
      dens_val     <- dvinecop(u_row, vine_j)
      li           <- ifelse(dens_val <= 0, -1e100, log(dens_val))
      log_incr[j]  <- li
    }
 
  log_incr
}


update_weights <- function(particles, log_lik) {
  
  w_prev <- vapply(particles, `[[`, numeric(1), 'w')
  log_w  <- log(w_prev) + log_lik
  log_w  <- log_w - max(log_w)
  w_new  <- exp(log_w); w_new <- w_new / sum(w_new)
  
  for (j in seq_along(particles))
    particles[[j]]$w <- w_new[j]
  
  particles
}

need_resample <- function(w, cfg, t, T) {
  ess(w) < cfg$ess_thr * cfg$M && t < T
}

resample_move <- function(particles, newAncestors, data_up_to_t, cl, type, cfg, skeleton=NULL, tr=NULL, temp_skel=NULL) {
  
  #idx       <- sample.int(M, M, TRUE, prob = w_new)
  particles <- particles[newAncestors]                  # cópia simples
  for (p in particles) {                       # reinicia pesos
    p$w <- 1 / cfg$M
  }
  
  mh_n_prop <- cfg$M * cfg$n_mh 
  mh_n_acc  <- 0
  clusterSetRNGStream(cl, 42) 
  
  tic("mh time")
  if (type=='standard'){
    mh_results <- parLapply(cl, seq_along(particles), function(i, particles_local, data_up_to_t, skeleton, cfg) {
      p <- particles_local[[i]]
      local_acc <- 0L
      for (k in seq_len(cfg$n_mh)) {
        p <- mh_step(p, data_up_to_t, skeleton, cfg) 
        if (isTRUE(p$last_accept)) {
          local_acc <- local_acc + 1L
        }
      }
      list(p = p, acc = local_acc)
    },
    particles, data_up_to_t, skeleton, cfg)
  } else if (type == "block") {
    mh_results <- parLapply(cl, seq_along(particles), function(i, particles_local, data_up_to_t, temp_skel, tr, cfg) {
      p <- particles_local[[i]]
      local_acc <- 0L
      for (k in seq_len(cfg$n_mh)) {
        p <- mh_step_in_tree(p, tr, data_up_to_t, temp_skel, cfg) # Use the new function signature
        if (isTRUE(p$last_accept)) {
          local_acc <- local_acc + 1L
        }
      }
      list(p = p, acc = local_acc)
    },
    particles, data_up_to_t, temp_skel, tr, cfg)
  } else {
    stop(sprintf("Unknown type: '%s'. Use 'standard' or 'block'.", type))
  }
  
  toc()
  #stopCluster(cl)
  
  mh_n_acc <- sum(vapply(mh_results, `[[`, integer(1), "acc"))
  particles <- lapply(mh_results, `[[`, "p")
  
  acc_pct <- 100 * mh_n_acc / mh_n_prop
  cat(sprintf(
    "MH acceptance = %4d / %4d  =  %.2f%%\n\n",
    mh_n_acc, mh_n_prop, acc_pct
  ))
  
  return(list(particles = particles, acc_pct = acc_pct))
}

## ───────────── Monte-Carlo utilities ────────────────────────────
w_mean <- function(x, w)  sum(w * x)
w_var  <- function(x, w, mu = w_mean(x, w))  sum(w * (x - mu)^2)
mc_se  <- function(x, w, ess = 1 / sum(w^2))
  sqrt(w_var(x, w) / ess)

w_quantile <- function(x, w, probs = c(0.025, 0.975)) {
  o <- order(x)
  x <- x[o]; w <- w[o] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) x[which.max(cw >= p)], numeric(1))
}

responsibility <- function(theta_mat, tau_vec, pi_vec, cfg){
  M <- nrow(theta_mat); K <- ncol(theta_mat)
  
  spike_sd <- matrix(tau_vec, nrow=M, ncol=K)
  slab_sd  <- spike_sd * cfg$c_slab
  pi_mat   <- matrix(pi_vec, nrow=M, ncol=K)
  
  spike <- pi_mat      * dnorm(theta_mat,0,spike_sd)
  slab  <- (1-pi_mat)  * dnorm(theta_mat,0,slab_sd)
  
  slab / (spike + slab + 1e-300)
}

## ───────────── diagnostic & plotting ────────────────────────────
diagnostic_report <- function(t, tr, U, particles, w_new,
                              cfg, q_probs  = c(0.025, 0.975)) {
  
  k_step <- cfg$k_step
  M <- cfg$M
  
  if (t %% k_step != 0L && t != nrow(U)) return(invisible())
  ## 1. unpack particle state -------------------------------------------------
  #gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
  theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
  rho_mat   <- tanh(theta_mat)
  
  #slab_w <- responsibility(theta_mat, cfg)   # M × K
  tau_vec <- sqrt(vapply(particles, function(p) p$tau2, numeric(1)))
  pi_vec  <- vapply(particles, function(p) p$pi, numeric(1))
  slab_w  <- responsibility(theta_mat, tau_vec, pi_vec, cfg)
  
  ## 2. normalise weights (essential) -----------------------------------------
  w_new <- w_new / sum(w_new)
  ess_t <- 1 / sum(w_new^2)
  
  tau_vec <- sqrt(vapply(particles, `[[`, numeric(1), "tau2"))
  mu_tau  <- w_mean(tau_vec, w_new)
  sd_tau  <- sqrt( w_var(tau_vec, w_new, mu_tau) )
  
  # print
  pi_mean <- w_mean(pi_vec, w_new)
  pi_sd   <- sqrt(w_var(pi_vec, w_new, pi_mean))
  cat(sprintf(" | π = %.3f ± %.3f", pi_mean, pi_sd))
  
  
  cat(sprintf(" | τ = %.4f ± %.4f\n", mu_tau, sd_tau))
  
  ## 3. global diagnostics ----------------------------------------------------
  edge_ct    <- rowSums(slab_w)
  mean_edges <- w_mean(edge_ct, w_new)
  se_edges   <- mc_se(edge_ct, w_new, ess_t)
  
  key_vec   <- apply(cbind(slab_w, theta_mat), 1L,
                     \(row) paste(row, collapse = ","))
  n_unique  <- length(unique(key_vec))
  
  dists     <- as.matrix(dist(theta_mat))
  avg_dist  <- mean(dists[lower.tri(dists)])
  
  cat(sprintf(
    "t = %4d | tr = %4d | ESS/M = %.3f | mean #edges = %.3f ± %.3f | unique = %d\n | Euclidean dist = %.4f\n",
    t, tr, ess_t / M, mean_edges, se_edges, n_unique, avg_dist))
  
  ## 4. per-edge summaries ----------------------------------------------------
  inc_prob <- colSums(slab_w  * w_new)
  
  edge_summ <- lapply(seq_along(inc_prob), function(e) {
    
    ## MC-SE of inclusion probability
    se_inc <- sqrt(inc_prob[e] * (1 - inc_prob[e]) / ess_t)
    
    ## conditional posterior of rho if edge present at all
    gamma_e <- slab_w[, e]
    rho_e   <- rho_mat[,  e]
    w_cond  <- w_new * gamma_e
    if (sum(w_cond) < 1e-12) {
      cat(sprintf("  Edge %2d : P(dep)=%.3f ± %.3f | never present\n",
                  e, inc_prob[e], se_inc))
      return(list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
                  mu_rho = NA, sd_rho = NA,
                  q025 = NA, q975 = NA, se_rho = NA))
    }
    
    w_cond <- w_cond / sum(w_cond)
    mu_rho <- sum(w_cond * rho_e)
    sd_rho <- sqrt( sum(w_cond * (rho_e - mu_rho)^2) )
    qs     <- w_quantile(rho_e, w_cond, q_probs)
    
    ## MC-SE of the mean (zero-padded variance / ESS)
    var_rho_zpad <- w_var(rho_e * slab_w[, e], w_new,
                          mu_rho * inc_prob[e])          # PATCH
    se_rho       <- sqrt(var_rho_zpad / ess_t)
    
    cat(sprintf(
      "  Edge %2d : P(dep)=%.3f ± %.3f | ρ = %.3f (SD %.3f, MC-SE %.3f) | 95%% CI = [%.3f, %.3f]\n",
      e, inc_prob[e], se_inc, mu_rho, sd_rho, se_rho, qs[1], qs[2]))
    
    list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
         mu_rho = mu_rho, sd_rho = sd_rho,
         q025 = qs[1], q975 = qs[2], se_rho = se_rho)
  })
  
  edge_df <- do.call(rbind.data.frame, edge_summ)
  
  
  
  invisible(list(ESS  = ess_t,
                 unique = n_unique,
                 euc = avg_dist,
                 tau_mean = mu_tau,      # pass back
                 tau_sd   = sd_tau,
                 edges = edge_df,
                 pi_mean = pi_mean,
                 pi_sd = pi_sd))
}



## Average L1 (absolute) distance between parameter vectors
avg_L1 <- function(theta_mat, w) {
  K  <- ncol(theta_mat)
  d  <- as.matrix(dist(theta_mat, method = "manhattan")) / K   # normalise
  ## Rao–Blackwellised average:  E_w  E_{w’≠w}  |θ−θ’|
  sum(w * colSums(d) / (1 - w))
}

systematic_resample <- function(w) {
  M <- length(w)
  u0 <- runif(1)/M
  cumw <- cumsum(w)
  (u0 + 0:(M-1)/M) |> 
    (\(u) findInterval(u, cumw)+1L)()
}


# Stratified resampling (O(M))
stratified_resample <- function(w) {
  M <- length(w)
  cumw <- cumsum(w)                     # cumulative weights
  u    <- ((seq_len(M) - 1) + runif(M)) / M   # one uniform draw per stratum
  findInterval(u, cumw) + 1L            # ancestor indices (1-based)
}



# ================================================================
#  EFFICIENT Predictive Function 
# ================================================================

compute_predictive_metrics <- function(u_obs, particles, skel, w_prev_for_prediction, cfg) {
  M <- length(particles)
  d <- ncol(u_obs)
  
  # --- 1. Compute Log Predictive Density (using original weighted particles) ---
  likelihoods <- sapply(particles, function(p) {
    vine_p <- vine_from_particle(p, skel, cfg)
    dvinecop(u_obs, vine_p)
  })
  weighted_likelihood <- sum(w_prev_for_prediction * likelihoods)
  log_pred_density <- log(weighted_likelihood + 1e-30)
  
  
  ## -- 2. posterior mean + s.e. of θ -----------------------------------
  theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
  rho_mat <- tanh(theta_mat)
  tau_vec   <- sqrt(vapply(particles, function(p) p$tau2, numeric(1)))
  #gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
  pi_vec <- vapply(particles, function(p) p$pi, numeric(1))
  slab_w <- responsibility(theta_mat, tau_vec, pi_vec, cfg)
  
  gamma_mean_now <- colSums(slab_w * w_prev_for_prediction)
  gamma_sd_now   <- sqrt(gamma_mean_now * (1 - gamma_mean_now))
  ## mean
  theta_mean_now <- colSums(rho_mat * w_prev_for_prediction)
  
  ## weighted var  = Σ w (θ-μ)²
  theta_var <- colSums(w_prev_for_prediction * (rho_mat - rep(theta_mean_now, each = nrow(rho_mat)))^2)
  #gamma_mean_now <- colSums(gamma_mat * w_prev_for_prediction)
  
  ## ESS for each step (scalar)  -> s.e. = √(var / ESS)
  ess <- 1 / sum(w_prev_for_prediction^2)
  theta_sd_now <- sqrt(theta_var)
  #gamma_sd_now  <- sqrt(gamma_mean_now * (1 - gamma_mean_now))
  # --- 4. Return all metrics in a single list ---
  return(list(
    log_pred_density = log_pred_density,
    theta_mean       = theta_mean_now,
    theta_se         = theta_sd_now,
    gamma_mean       = gamma_mean_now,
    gamma_se      = gamma_sd_now
  ))
}



