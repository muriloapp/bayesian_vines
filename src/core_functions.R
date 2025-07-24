library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(matrixStats)


# spike_sd_from_tau <- function(tau){tau}
# 
# slab_sd_from_tau  <- function(tau, cfg){cfg$c_slab * tau}
# 
# log_prior <- function(p,cfg){
#   tau <- sqrt(p$tau2)
#   pi  <- p$pi                             
#   
#   ll_mix <- sum(
#     log( pi      * dnorm(p$theta,0, spike_sd_from_tau(tau)) +
#            (1-pi)  * dnorm(p$theta,0, slab_sd_from_tau(tau,cfg)) + 1e-300))
#   
#   ll_tau <- if(cfg$tau_prior=="fixed") 0
#   else dinvgamma(p$tau2, cfg$a0, cfg$b0, log=TRUE)
#   
#   ll_pi  <- if(cfg$pi_prior=="fixed") 0
#   else dbeta(pi, cfg$a_pi, cfg$b_pi, log=TRUE)
#   
#   ll_mix + ll_tau + ll_pi
# }

## ---------------------------------------------------------------
##  Exponential prior   p(γ,θ) ∝ exp(-λ Σγ)  ×  ∏_{γ=1}  N(0,σ²)
## ---------------------------------------------------------------
## ---------------------------------------------------------------
##  Truncated normal sampler  N(μ,σ²) ∩ (a,b)
## ---------------------------------------------------------------
rtnorm1 <- function(mu, sd, a = -0.99, b = 0.99) {
  repeat {
    x <- rnorm(1, mu, sd)
    if (x > a && x < b) return(x)
  }
}



log_prior <- function(p, cfg) {
  ##  Uniform prior / proposal density for ρ  ∈ (-0.99,0.99)
  log_unif_rho <- -log(1.98)
  
  -cfg$lambda * sum(p$gamma) +     # complexity penalty
    sum(p$gamma) * log_unif_rho       # uniform density for active ρ
}


# new_particle <- function(cfg){
#   # --- τ² block unchanged -----------------------------------
#   tau2 <- if(cfg$tau_prior=="fixed") {cfg$tau0^2}
#   else {rinvgamma(1, cfg$a0, cfg$b0)}
#   tau  <- sqrt(tau2)
#   
#   # --- NEW π draw ------------------------------------------
#   pi   <- if(cfg$pi_prior=="fixed") cfg$pi0
#   else rbeta(1, cfg$a_pi, cfg$b_pi)
#   
#   # --- θ draws from *mixture* -------------------------------
#   is_slab <- rbinom(cfg$K, 1, 1 - pi)
#   
#   aux_slab <- slab_sd_from_tau(tau,cfg)
#   aux_spike <- spike_sd_from_tau(tau)
#   theta   <- rnorm(cfg$K, 0,
#                    ifelse(is_slab,
#                           aux_slab,
#                           aux_spike))
#   
#   list(theta = theta,
#        tau2  = tau2,
#        pi    = pi,        # << store π
#        w     = 1/cfg$M,
#        last_accept = FALSE,
#        aux_slab = aux_slab,
#        aux_spike = aux_spike,
#        is_slab = is_slab
#        )
# }

## core_functions.R  ──────────────────────────────────────────────────
##  Family codes: 0 = independence, 1 = Gaussian  (extend later)

## ---------------------------------------------------------------
##  Prior over models:   P(γ)  ∝  exp(-λ Σ γ)
##  For equal “no-preference” switching probs we still need γ
## ---------------------------------------------------------------
new_particle <- function(cfg) {
  ##  Draw γ_i independently with P(γ_i = 1) =  exp(-λ)/(1+exp(-λ))
  p_gauss <- exp(-cfg$lambda) / (1 + exp(-cfg$lambda))
  gamma <- rbinom(cfg$K, 1, p_gauss)
  
  theta <- numeric(cfg$K)
  if (any(gamma == 1L))
    theta[gamma == 1L] <- runif(sum(gamma), -0.99, 0.99)
  
  list(theta = theta,    # these ARE the correlations ρ
       gamma = gamma,
       w     = 1 / cfg$M,
       last_accept = FALSE)
}


vine_from_particle <- function(p, skel) {
    pcs <- lapply(skel$pair_copulas, function(lvl) lapply(lvl, identity))
    idx <- 1L
    for (tr in seq_along(pcs)) {
      for (ed in seq_along(pcs[[tr]])) {
        pcs[[tr]][[ed]] <- bicop_dist("gaussian", parameters = tanh(p$theta[idx]))
        idx <- idx + 1L
      }
    }
    vinecop_dist(pcs, structure = skel$structure)
}


# fast_vine_from_particle <- function(p, skel) {
#   
#   rho <- tanh(p$theta)
#   pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)  # or FALSE in parallel
#   
#   idx <- 1L
#   for (tr in seq_along(pcs)){
#     for (ed in seq_along(pcs[[tr]])) {
#       pcs[[tr]][[ed]]$parameters <- matrix(rho[idx], 1, 1)    # ← key change
#       pcs[[tr]][[ed]]$npars      <- 1L                        # be explicit
#       idx <- idx + 1L
#     }}
#   
#   vinecop_dist(pcs, skel$structure)
# }

fast_vine_from_particle <- function(p, skel) {
  
  pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)  # << fast
  idx <- 1L
  for (tr in seq_along(pcs)) {
    for (ed in seq_along(pcs[[tr]])) {
      
      if (p$gamma[idx] == 0L) {
        pcs[[tr]][[ed]] <- bicop_dist("indep")
      } else {
        pcs[[tr]][[ed]]$family     <- "gaussian"
        pcs[[tr]][[ed]]$parameters <- matrix(p$theta[idx], 1, 1)
        pcs[[tr]][[ed]]$npars      <- 1L
      }
      idx <- idx + 1L
    }
  }
  vinecop_dist(pcs, skel$structure)
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
  vine_prop <- fast_vine_from_particle(prop, temp_skel)
  prop_ll   <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
                            .Machine$double.eps)))
  
  vine_curr <- fast_vine_from_particle(p, temp_skel)
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

  return(p)
}


# mh_step <- function(p, data_up_to_t, skeleton, cfg) {
#   prop <- p
#   prop$theta <- p$theta + rnorm(cfg$K, 0, cfg$step_sd)
#   
#   vine_prop <- fast_vine_from_particle(prop, skeleton)
#   vine_curr <- fast_vine_from_particle(p,    skeleton)
#   
#   ll_prop <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
#                           .Machine$double.eps)))
#   ll_curr <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
#                           .Machine$double.eps)))
#   
#   # log_acc <- (ll_prop + log_prior(prop$theta, cfg)) -
#   #   (ll_curr + log_prior(p$theta, cfg))
#   log_acc <- (ll_prop + log_prior(prop, cfg)) -      # ← pass WHOLE prop
#     (ll_curr + log_prior(p,    cfg))       # ← pass WHOLE p
#   
#   if (log(runif(1)) < log_acc) {
#     p$theta <- prop$theta
#     p$last_accept <- TRUE
#   } else p$last_accept <- FALSE
#   
#   ## NEW ---------------------------------------------------------------
#   p <- update_tau2(p, cfg)
#   p <- update_pi(p,cfg)        # << NEW
#   p
# }
rtnorm_vec <- function(mu_vec, sd, a = -0.99, b = 0.99) {
  x <- mu_vec + sd * rnorm(length(mu_vec))
  pmin(pmax(x, a), b)                     # hard clip (reflection is slower)
}

mh_step <- function(p, data_up_to_t, skeleton, cfg) {
  
  ## ---------- 1. PROPOSAL (vectorised) ------------------------------
  prop <- p
  
  ## (a)  γ-flips:  Bernoulli(q_flip)  for every edge
  flip_mask          <- runif(cfg$K) < cfg$q_flip
  prop$gamma[flip_mask] <- 1L - prop$gamma[flip_mask]
  
  ## (b)  handle BIRTHS & DEATHS
  births <- (prop$gamma == 1L & p$gamma == 0L)
  deaths <- (prop$gamma == 0L & p$gamma == 1L)
  if (any(births))
    prop$theta[births] <- runif(sum(births), -0.99, 0.99)
  if (any(deaths))
    prop$theta[deaths] <- 0
  
  ## (c)  Random-walk on all active θ  (after birth/death)
  active <- which(prop$gamma == 1L)
  if (length(active))
    prop$theta[active] <- rtnorm_vec(prop$theta[active], cfg$step_sd)
  
  ## ---------- 2. MH ACCEPTANCE --------------------------------------
  vine_prop <- fast_vine_from_particle(prop, skeleton)
  vine_curr <- fast_vine_from_particle(p,    skeleton)
  
  ll_prop <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
                          .Machine$double.eps)))
  ll_curr <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
                          .Machine$double.eps)))
  
  ## Symmetric proposal:  q(m→m′) = q(m′→m)
  ##   – γ-flips independent Bernoulli(q_flip) in both directions
  ##   – RW step is symmetric around current θ
  ##   – births draw ρ from Uniform identical to its prior
  log_acc <- (ll_prop + log_prior(prop, cfg)) -
    (ll_curr + log_prior(p,    cfg))
  
  if (log(runif(1)) < log_acc) {
    prop$last_accept <- TRUE
    return(prop)
  } else {
    p$last_accept <- FALSE
    return(p)
  }
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

compute_adapt_step_sd <- function(cfg, acc_pct, lambda = 0.25, target_acc = 0.30, sd_min = 0.01, sd_max=0.05){
  log_sd_new <- log(cfg$step_sd) + lambda * (acc_pct/100 - target_acc)
  step_sd <- pmin(pmax(exp(log_sd_new), sd_min), sd_max)
  cat(sprintf(
    "step_sd = %4f \n\n",
    step_sd
  ))
  return(step_sd)
}


calculate_log_lik_tree_tr <- function(particle, skel_tr, u_row, t, tr, tr_prev, skeletons_by_tr, cfg) {
  # If it's the first tree, calculate its log-density
  if (tr == 1) {
    vine_j <- fast_vine_from_particle(particle, skel_tr)
    dens_val <- dvinecop(u_row, vine_j)
    return(log(dens_val + 1e-100))
  }

  # For subsequent trees, we need the density of the current and previous model
  # The pre-computation of skeletons_by_tr makes this much faster.
  skel_tr_minus_1 <- skeletons_by_tr[[tr_prev]]

  vine_j_tr <- fast_vine_from_particle(particle, skel_tr)
  dens_val_tr <- dvinecop(u_row, vine_j_tr)

  vine_j_prev <- fast_vine_from_particle(particle, skel_tr_minus_1)
  dens_val_prev <- dvinecop(u_row, vine_j_prev)

  # Return the log-difference
  return(log(dens_val_tr + 1e-100) - log(dens_val_prev + 1e-100))
}



compute_log_incr <- function(particles, u_row, skeleton, cfg) {
  log_incr <- numeric(cfg$M)
  
  for (j in seq_along(particles)) {
      p <- particles[[j]]
      #vine_from_particle(p, skeleton)
      vine_j     <- fast_vine_from_particle(p, skeleton)
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


## ───────────── diagnostic & plotting ────────────────────────────
# diagnostic_report <- function(t, tr, U, particles, w_new,
#                               cfg, q_probs  = c(0.025, 0.975)) {
#   
#   k_step <- cfg$k_step
#   M <- cfg$M
#   
#   if (t %% k_step != 0L && t != nrow(U)) return(invisible())
#   ## 1. unpack particle state -------------------------------------------------
#   #gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
#   theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
#   rho_mat   <- tanh(theta_mat)
#   
#   #slab_w <- responsibility(theta_mat, cfg)   # M × K
#   tau_vec <- sqrt(vapply(particles, function(p) p$tau2, numeric(1)))
#   pi_vec  <- vapply(particles, function(p) p$pi, numeric(1))
#   slab_w  <- responsibility(theta_mat, tau_vec, pi_vec, cfg)
#   
#   ## 2. normalise weights (essential) -----------------------------------------
#   w_new <- w_new / sum(w_new)
#   ess_t <- 1 / sum(w_new^2)
#   
#   tau_vec <- sqrt(vapply(particles, `[[`, numeric(1), "tau2"))
#   mu_tau  <- w_mean(tau_vec, w_new)
#   sd_tau  <- sqrt( w_var(tau_vec, w_new, mu_tau) )
#   
#   # print
#   pi_mean <- w_mean(pi_vec, w_new)
#   pi_sd   <- sqrt(w_var(pi_vec, w_new, pi_mean))
#   cat(sprintf(" | π = %.3f ± %.3f", pi_mean, pi_sd))
#   
#   
#   cat(sprintf(" | τ = %.4f ± %.4f\n", mu_tau, sd_tau))
#   
#   ## 3. global diagnostics ----------------------------------------------------
#   edge_ct    <- rowSums(slab_w)
#   mean_edges <- w_mean(edge_ct, w_new)
#   se_edges   <- mc_se(edge_ct, w_new, ess_t)
#   
#   key_vec   <- apply(cbind(theta_mat), 1L,
#                      \(row) paste(row, collapse = ","))
#   n_unique  <- length(unique(key_vec))
#   
#   dists     <- as.matrix(dist(theta_mat))
#   avg_dist  <- mean(dists[lower.tri(dists)])
#   
#   cat(sprintf(
#     "t = %4d | tr = %4d | ESS/M = %.3f | mean #edges = %.3f ± %.3f | unique = %d\n | Euclidean dist = %.4f\n",
#     t, tr, ess_t / M, mean_edges, se_edges, n_unique, avg_dist))
#   
#   ## 4. per-edge summaries ----------------------------------------------------
#   inc_prob <- colSums(slab_w  * w_new)
#   
#   edge_summ <- lapply(seq_along(inc_prob), function(e) {
#     
#     ## MC-SE of inclusion probability
#     se_inc <- sqrt(inc_prob[e] * (1 - inc_prob[e]) / ess_t)
#     
#     ## conditional posterior of rho if edge present at all
#     gamma_e <- slab_w[, e]
#     rho_e   <- rho_mat[,  e]
#     w_cond  <- w_new * gamma_e
#     if (sum(w_cond) < 1e-4) {
#       cat(sprintf("  Edge %2d : P(dep)=%.3f ± %.3f | never present\n",
#                   e, inc_prob[e], se_inc))
#       return(list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
#                   mu_rho = NA, sd_rho = NA,
#                   q025 = NA, q975 = NA, se_rho = NA))
#     }
#     
#     w_cond <- w_cond / sum(w_cond)
#     mu_rho <- sum(w_cond * rho_e)
#     sd_rho <- sqrt( sum(w_cond * (rho_e - mu_rho)^2) )
#     qs     <- w_quantile(rho_e, w_cond, q_probs)
#     
#     ## MC-SE of the mean (zero-padded variance / ESS)
#     var_rho_zpad <- w_var(rho_e * slab_w[, e], w_new,
#                           mu_rho * inc_prob[e])          # PATCH
#     se_rho       <- sqrt(var_rho_zpad / ess_t)
#     
#     cat(sprintf(
#       "  Edge %2d : P(dep)=%.3f ± %.3f | ρ = %.3f (SD %.3f, MC-SE %.3f) | 95%% CI = [%.3f, %.3f]\n",
#       e, inc_prob[e], se_inc, mu_rho, sd_rho, se_rho, qs[1], qs[2]))
#     
#     list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
#          mu_rho = mu_rho, sd_rho = sd_rho,
#          q025 = qs[1], q975 = qs[2], se_rho = se_rho)
#   })
#   
#   edge_df <- do.call(rbind.data.frame, edge_summ)
#   
#   
#   
#   invisible(list(ESS  = ess_t,
#                  unique = n_unique,
#                  euc = avg_dist,
#                  tau_mean = mu_tau,      # pass back
#                  tau_sd   = sd_tau,
#                  edges = edge_df,
#                  pi_mean = pi_mean,
#                  pi_sd = pi_sd))
# }


# diagnostic_report <- function(t, tr, U, particles, w_new,
#                               cfg, q_probs = c(0.025, 0.975)) {
#   k_step     <- cfg$k_step
#   M          <- cfg$M
#   print_flag <- (t %% k_step == 0L) || (t == nrow(U))
#   
#   ## 1. unpack particle state -------------------------------------------------
#   theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
#   rho_mat   <- (theta_mat)
#   
#   ## 2. normalise weights (essential) ----------------------------------------
#   w_new <- w_new / sum(w_new)
#   ess_t <- 1 / sum(w_new^2)
#   
#   
#   ## 3. global diagnostics ----------------------------------------------------
#   edge_ct    <- rowSums(slab_w)
#   mean_edges <- w_mean(edge_ct, w_new)
#   se_edges   <- mc_se(edge_ct, w_new, ess_t)
#   
#   key_vec   <- apply(theta_mat, 1L, \(row) paste(row, collapse = ","))
#   n_unique  <- length(unique(key_vec))
#   
#   dists     <- as.matrix(stats::dist(theta_mat))
#   avg_dist  <- mean(dists[lower.tri(dists)])
#   
#   if (print_flag) {
#     cat(sprintf(
#       "t = %4d | tr = %4d | ESS/M = %.3f | mean #edges = %.3f ± %.3f | unique = %d\n | Euclidean dist = %.4f\n",
#       t, tr, ess_t / M, mean_edges, se_edges, n_unique, avg_dist))
#   }
#   
#   ## 4. per-edge summaries ----------------------------------------------------
#   inc_prob <- colSums(slab_w * w_new)
#   
#   edge_summ <- lapply(seq_along(inc_prob), function(e) {
#     se_inc  <- sqrt(inc_prob[e] * (1 - inc_prob[e]) / ess_t)
#     
#     gamma_e <- slab_w[, e]
#     rho_e   <- rho_mat[,  e]
#     w_cond  <- w_new * gamma_e
#     
#     if (sum(w_cond) < 1e-4) {
#       if (print_flag)
#         cat(sprintf("  Edge %2d : P(dep)=%.3f ± %.3f | never present\n",
#                     e, inc_prob[e], se_inc))
#       return(list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
#                   mu_rho = NA, sd_rho = NA,
#                   q025 = NA, q975 = NA, se_rho = NA))
#     }
#     
#     w_cond <- w_cond / sum(w_cond)
#     mu_rho <- sum(w_cond * rho_e)
#     sd_rho <- sqrt(sum(w_cond * (rho_e - mu_rho)^2))
#     qs     <- w_quantile(rho_e, w_cond, q_probs)
#     
#     var_rho_zpad <- w_var(rho_e * slab_w[, e], w_new,
#                           mu_rho * inc_prob[e])
#     se_rho <- sqrt(var_rho_zpad / ess_t)
#     
#     if (print_flag) {
#       cat(sprintf(
#         "  Edge %2d : P(dep)=%.3f ± %.3f | ρ = %.3f (SD %.3f, MC-SE %.3f) | 95%% CI = [%.3f, %.3f]\n",
#         e, inc_prob[e], se_inc, mu_rho, sd_rho, se_rho, qs[1], qs[2]))
#     }
#     
#     list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
#          mu_rho = mu_rho, sd_rho = sd_rho,
#          q025 = qs[1], q975 = qs[2], se_rho = se_rho)
#   })
#   
#   edge_df <- do.call(rbind.data.frame, edge_summ)
#   
#   ## 5. invisibly return everything ------------------------------------------
#   invisible(list(ESS       = ess_t,
#                  unique    = n_unique,
#                  euc       = avg_dist,
#                  edges     = edge_df))
# }

diagnostic_report <- function(t, tr, U, particles, w_new,
                              cfg, q_probs = c(0.025, 0.975)) {
  
  k_step     <- cfg$k_step
  M          <- cfg$M
  print_flag <- (t %% k_step == 0L) || (t == nrow(U))
  
  ## ---------- 1. unpack particle state -------------------------------
  theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
  rho_mat   <- theta_mat                       # correlations
  gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
  
  ## ---------- 2. normalise weights -----------------------------------
  w_new <- w_new / sum(w_new)
  ess_t <- 1 / sum(w_new^2)
  
  ## ---------- 3. global diagnostics ----------------------------------
  edge_ct    <- rowSums(gamma_mat)                 # #Gaussian edges / particle
  mean_edges <- w_mean(edge_ct, w_new)
  se_edges   <- mc_se(edge_ct, w_new, ess_t)
  sparsity   <- 1 - mean_edges / cfg$K            # 1 = fully sparse
  
  key_vec   <- apply(theta_mat, 1, \(row) paste(row, collapse = ","))
  n_unique  <- length(unique(key_vec))
  
  dists     <- as.matrix(stats::dist(theta_mat))
  avg_dist  <- mean(dists[lower.tri(dists)])
  
  if (print_flag) {
    cat(sprintf(
      "t=%4d | tr=%2d | ESS/M=%.3f | dep.edges=%.2f ± %.2f | sparsity=%.3f | unique=%d | L2=%.4f\n",
      t, tr, ess_t / M, mean_edges, se_edges, sparsity, n_unique, avg_dist))
  }
  
  ## ---------- 4. per-edge summaries ----------------------------------
  inc_prob <- colSums(gamma_mat * w_new)          # P(γ_e = 1)
  
  edge_summ <- lapply(seq_along(inc_prob), function(e) {
    
    se_inc <- sqrt(inc_prob[e] * (1 - inc_prob[e]) / ess_t)
    
    active  <- gamma_mat[, e] == 1L
    if (!any(active)) {
      if (print_flag)
        cat(sprintf("  Edge %2d : P(dep)=%.3f ± %.3f | never active\n",
                    e, inc_prob[e], se_inc))
      return(list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
                  mu_rho = NA, sd_rho = NA,
                  q025 = NA, q975 = NA))
    }
    
    rho_e <- rho_mat[active, e]
    w_e   <- w_new[active] / sum(w_new[active])
    
    mu_rho <- sum(w_e * rho_e)
    sd_rho <- sqrt(sum(w_e * (rho_e - mu_rho)^2))
    qs     <- w_quantile(rho_e, w_e, q_probs)
    
    if (print_flag)
      cat(sprintf(
        "  Edge %2d : P(dep)=%.3f ± %.3f | ρ = %.3f (SD %.3f) | 95%% CI = [%.3f, %.3f]\n",
        e, inc_prob[e], se_inc, mu_rho, sd_rho, qs[1], qs[2]))
    
    list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
         mu_rho = mu_rho, sd_rho = sd_rho,
         q025 = qs[1], q975 = qs[2])
  })
  
  edge_df <- do.call(rbind.data.frame, edge_summ)
  
  invisible(list(
    ESS      = ess_t,
    unique   = n_unique,
    euc      = avg_dist,
    sparsity = sparsity,
    edges    = edge_df
  ))
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

# compute_predictive_metrics <- function(u_obs, particles, skel, w_prev_for_prediction, cfg) {
#   M <- length(particles)
#   d <- ncol(u_obs)
#   
#   # --- 1. Compute Log Predictive Density (using original weighted particles) ---
#   likelihoods <- sapply(particles, function(p) {
#     vine_p <- fast_vine_from_particle(p, skel)
#     dvinecop(u_obs, vine_p)
#   })
#   weighted_likelihood <- sum(w_prev_for_prediction * likelihoods)
#   log_pred_density <- log(weighted_likelihood + 1e-30)
#   
#   
#   ## -- 2. posterior mean + s.e. of θ -----------------------------------
#   theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
#   rho_mat <- (theta_mat)
#   #gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
# 
#   gamma_mean_now <- colSums(slab_w * w_prev_for_prediction)
#   gamma_sd_now   <- sqrt(gamma_mean_now * (1 - gamma_mean_now))
#   ## mean
#   theta_mean_now <- colSums(rho_mat * w_prev_for_prediction)
#   
#   ## weighted var  = Σ w (θ-μ)²
#   theta_var <- colSums(w_prev_for_prediction * (rho_mat - rep(theta_mean_now, each = nrow(rho_mat)))^2)
#   #gamma_mean_now <- colSums(gamma_mat * w_prev_for_prediction)
#   
#   ## ESS for each step (scalar)  -> s.e. = √(var / ESS)
#   ess <- 1 / sum(w_prev_for_prediction^2)
#   theta_sd_now <- sqrt(theta_var)
#   #gamma_sd_now  <- sqrt(gamma_mean_now * (1 - gamma_mean_now))
#   # --- 4. Return all metrics in a single list ---
#   return(list(
#     log_pred_density = log_pred_density,
#     theta_mean       = theta_mean_now,
#     theta_se         = theta_sd_now,
#     gamma_mean       = gamma_mean_now,
#     gamma_se      = gamma_sd_now
#   ))
# }


compute_predictive_metrics <- function(u_obs, particles, skel,
                                       w_prev_for_prediction, cfg,
                                       q_probs = c(0.025, 0.975)) {
  
  ## ---------- 1. predictive log-density -------------------------------
  lik <- vapply(
    particles,
    function(p) {
      vine_p <- fast_vine_from_particle(p, skel)
      dvinecop(u_obs, vine_p)
    },
    numeric(1)
  )
  log_pred_density <- log(sum(w_prev_for_prediction * lik) + 1e-30)
  
  ## ---------- 2. posterior summaries ---------------------------------
  theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
  rho_mat   <- theta_mat                          # already ρ in (-0.99,0.99)
  
  gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
  
  ## inclusion probs
  gamma_mean_now <- colSums(gamma_mat * w_prev_for_prediction)
  gamma_sd_now   <- sqrt(gamma_mean_now * (1 - gamma_mean_now))
  
  ## θ means & MC-SE
  theta_mean_now <- colSums(rho_mat * w_prev_for_prediction)
  theta_var_now  <- colSums(
    w_prev_for_prediction *
      (rho_mat - rep(theta_mean_now, each = nrow(rho_mat)))^2
  )
  theta_sd_now <- sqrt(theta_var_now)
  
  list(
    log_pred_density = log_pred_density,
    theta_mean = theta_mean_now,
    theta_se   = theta_sd_now,
    gamma_mean = gamma_mean_now,
    gamma_se   = gamma_sd_now
  )
}




