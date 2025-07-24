library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(matrixStats)


FAM_INFO <- data.frame(
  name = c("indep", "gaussian", "bb1"),
  code = c(0L, 1L, 2L),
  npar = c(0L, 1L, 2L)
)

# 0 = indep, 1 = gaussian, 2 = bb1
FAM_INDEP   <- 0L
FAM_GAUSS   <- 1L
FAM_BB1     <- 2L

## ------------------------------------------------------------------
##  prepare ONE template per family 
## ------------------------------------------------------------------
T_INDEP  <- bicop_dist("indep")
T_GAUSS  <- bicop_dist("gaussian", parameters = 0)
T_BB1    <- bicop_dist("bb1",      parameters = c(1, 2))  # dummy θ,δ




active_fams <- function(cfg) FAM_INFO[FAM_INFO$name %in% cfg$families, ]


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



log_prior_edge <- function(fam, th1, th2, cfg) {
  
  if (fam == FAM_INDEP) return(0)           # constant
  
  if (fam == FAM_GAUSS) {
    lp <- log(1/1.98)                       # Uniform ρ
    return(lp - cfg$lambda)                 # + penalty
    
  } else {  # BB1
    lam <- bb1_par2tail(th1, th2)
    lp_tail <- dbeta(lam[1], 2, 2, log=TRUE) +
      dbeta(lam[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb1_log_jacobian(lam[1], lam[2]) - 2*cfg$lambda
    return(lp)
  }
}

log_prior <- function(p, cfg) {
  sum(mapply(log_prior_edge,
             p$fam, p$th1, p$th2,
             MoreArgs = list(cfg = cfg)))
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
# new_particle <- function(cfg) {
#   ##  Draw γ_i independently with P(γ_i = 1) =  exp(-λ)/(1+exp(-λ))
#   p_gauss <- exp(-cfg$lambda) / (1 + exp(-cfg$lambda))
#   gamma <- rbinom(cfg$K, 1, p_gauss)
#   
#   theta <- numeric(cfg$K)
#   if (any(gamma == 1L))
#     theta[gamma == 1L] <- runif(sum(gamma), -0.99, 0.99)
#   
#   list(theta = theta,    # these ARE the correlations ρ
#        gamma = gamma,
#        w     = 1 / cfg$M,
#        last_accept = FALSE)
# }

new_particle <- function(cfg) {
  
  fam_tbl <- active_fams(cfg)
  codes   <- fam_tbl$code
  npar    <- fam_tbl$npar
  
  ## complexity weight  exp(-λ d_e)
  weights <- exp(-cfg$lambda * npar)
  prob    <- weights / sum(weights)
  
  fam <- sample(codes, cfg$K, replace = TRUE, prob = prob)
  
  th1 <- numeric(cfg$K); th2 <- numeric(cfg$K)
  
  for (k in seq_len(cfg$K)) {
    if (fam[k] == 1L) {                     # Gaussian
      th1[k] <- runif(1, -0.99, 0.99)
    } else if (fam[k] == 2L) {              # BB1
      lamU <- rbeta(1, 2, 2); lamL <- rbeta(1, 2, 2)
      tp   <- bb1_tail2par(lamL, lamU)
      th1[k] <- tp[1]; th2[k] <- tp[2]
    }
  }
  
  list(fam = fam, th1 = th1, th2 = th2, w = 1 / cfg$M, last_accept = FALSE)
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

## ----- BB1  (λ  ↔  θ,δ) ----------------------------------------------
bb1_tail2par <- function(lambdaL, lambdaU) {
  theta <- log(2-lambdaU)/(-log(lambdaL))          # θ > 0
  delta <- log(2) / log(2-lambdaU)             # δ > 1
  c(theta, delta)
}

bb1_par2tail <- function(theta, delta) {
  lambdaL <- 2^(-1 / (delta*theta))
  lambdaU <- 2 - 2^(1 / delta)
  c(lambdaL, lambdaU)
}

# bb1_log_jacobian <- function(lambdaL, lambdaU) {
#   -log( (2 - lambdaU) * lambdaL * (log(2))^2 )
# }

bb1_log_jacobian <- function(lambdaL, lambdaU) {
  log2 <- log(2)
  denom <- (log(1 / lambdaL))^2 * lambdaL * log(2 - lambdaU) * (2 - lambdaU)
  log(log2) - log(denom)
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

# fast_vine_from_particle <- function(p, skel) {
#   
#   pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)  # << fast
#   idx <- 1L
#   for (tr in seq_along(pcs)) {
#     for (ed in seq_along(pcs[[tr]])) {
#       
#       if (p$gamma[idx] == 0L) {
#         pcs[[tr]][[ed]] <- bicop_dist("indep")
#       } else {
#         pcs[[tr]][[ed]]$family     <- "gaussian"
#         pcs[[tr]][[ed]]$parameters <- matrix(p$theta[idx], 1, 1)
#         pcs[[tr]][[ed]]$npars      <- 1L
#       }
#       idx <- idx + 1L
#     }
#   }
#   vinecop_dist(pcs, skel$structure)
# }
# indep <- bicop_dist("indep")
# 
# fast_vine_from_particle <- function(p, skel) {
#   
#   pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)
#   idx <- 1L
#   for (tr in seq_along(pcs)) {
#     for (ed in seq_along(pcs[[tr]])) {
#       
#       f <- p$fam[idx]
#       if (f == FAM_INDEP) {
#         pcs[[tr]][[ed]] <- indep
#         
#       } else if (f == FAM_GAUSS) {
#         pcs[[tr]][[ed]]$family     <- "gaussian"
#         pcs[[tr]][[ed]]$parameters <- matrix(p$th1[idx], 1, 1)
#         pcs[[tr]][[ed]]$npars      <- 1L
#         
#       } else {  # BB1
#         pcs[[tr]][[ed]]$family     <- "bb1"
#         pcs[[tr]][[ed]]$parameters <- matrix(c(p$th1[idx], p$th2[idx]), 1, 2)
#         pcs[[tr]][[ed]]$npars      <- 2L
#         #pcs[[tr]][[ed]] <-  bicop_dist("bb1", 0, c(p$th1[idx], p$th2[idx]))
#       }
#       idx <- idx + 1L
#     }
#   }
#   
#   vinecop_dist(pcs, skel$structure)
# }

sanitize_bb1 <- function(theta, delta,
                         eps   = 1e-6,   # lower-bound cushion
                         upper = 7 - 1e-6) {
  
  theta <- pmin(pmax(theta, eps),   upper)   # (0 , 7]
  delta <- pmin(pmax(delta, 1+eps), upper)   # (1 , 7]
  c(theta, delta)
}

fast_vine_from_particle <- function(p, skel) {
  
  pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)
  idx <- 1L
  
  for (tr in seq_along(pcs)) {
    for (ed in seq_along(pcs[[tr]])) {
      
      f <- p$fam[idx]
      
      if (f == FAM_INDEP) {
        pcs[[tr]][[ed]] <- T_INDEP           # shallow copy by assignment
        
      } else if (f == FAM_GAUSS) {
        bc <- rlang::duplicate(T_GAUSS, shallow = TRUE)
        bc$parameters[1] <- p$th1[idx]       # update ρ
        pcs[[tr]][[ed]] <- bc
        
      } else {                               # BB1
        pars <- sanitize_bb1(p$th1[idx], p$th2[idx])   # ⟵ NEW
        p$th1[idx] <- pars[1]                          # keep particle in sync
        p$th2[idx] <- pars[2]
        
        bc <- rlang::duplicate(T_BB1, shallow = TRUE)
        bc$parameters[1:2] <- pars
        pcs[[tr]][[ed]] <- bc
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

# mh_step <- function(p, data_up_to_t, skeleton, cfg) {
#   
#   ## ---------- 1. PROPOSAL (vectorised) ------------------------------
#   prop <- p
#   
#   ## (a)  γ-flips:  Bernoulli(q_flip)  for every edge
#   flip_mask          <- runif(cfg$K) < cfg$q_flip
#   prop$gamma[flip_mask] <- 1L - prop$gamma[flip_mask]
#   
#   ## (b)  handle BIRTHS & DEATHS
#   births <- (prop$gamma == 1L & p$gamma == 0L)
#   deaths <- (prop$gamma == 0L & p$gamma == 1L)
#   if (any(births))
#     prop$theta[births] <- runif(sum(births), -0.99, 0.99)
#   if (any(deaths))
#     prop$theta[deaths] <- 0
#   
#   ## (c)  Random-walk on all active θ  (after birth/death)
#   active <- which(prop$gamma == 1L)
#   if (length(active))
#     prop$theta[active] <- rtnorm_vec(prop$theta[active], cfg$step_sd)
#   
#   ## ---------- 2. MH ACCEPTANCE --------------------------------------
#   vine_prop <- fast_vine_from_particle(prop, skeleton)
#   vine_curr <- fast_vine_from_particle(p,    skeleton)
#   
#   ll_prop <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
#                           .Machine$double.eps)))
#   ll_curr <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
#                           .Machine$double.eps)))
#   
#   ## Symmetric proposal:  q(m→m′) = q(m′→m)
#   ##   – γ-flips independent Bernoulli(q_flip) in both directions
#   ##   – RW step is symmetric around current θ
#   ##   – births draw ρ from Uniform identical to its prior
#   log_acc <- (ll_prop + log_prior(prop, cfg)) -
#     (ll_curr + log_prior(p,    cfg))
#   
#   if (log(runif(1)) < log_acc) {
#     prop$last_accept <- TRUE
#     return(prop)
#   } else {
#     p$last_accept <- FALSE
#     return(p)
#   }
# }

mh_step <- function(p, data_up_to_t, skeleton, cfg) {
  
  K <- cfg$K
  fam_tbl <- active_fams(cfg)            # data.frame(name, code, npar)
  codes   <- fam_tbl$code
  
  ## ---------- 1. PROPOSAL (vectorised) ------------------------------
  prop <- p
  
  ## (a) family flips with probability q_flip per edge
  flip_mask <- runif(K) < cfg$q_flip
  if (any(flip_mask)) {
    idxs <- which(flip_mask)
    for (idx in idxs) {
      old_code <- prop$fam[idx]
      
      ## choose a *different* family uniformly among the active set
      new_code <- sample(codes[codes != old_code], 1L)
      prop$fam[idx] <- new_code
      
      ## draw parameters from that family's PRIOR
      if (new_code == FAM_INDEP) {
        prop$th1[idx] <- 0; prop$th2[idx] <- 0
        
      } else if (new_code == FAM_GAUSS) {
        prop$th1[idx] <- runif(1, -0.99, 0.99)
        prop$th2[idx] <- 0
        
      } else if (new_code == FAM_BB1) {
        lambdaU <- rbeta(1, 2, 2)
        lambdaL <- rbeta(1, 2, 2)
        tp      <- bb1_tail2par(lambdaL, lambdaU)
        prop$th1[idx] <- tp[1]      # θ
        prop$th2[idx] <- tp[2]      # δ
      }
    }
  }
  
  ## (b) Random-walk on parameters of *current* families
  ## Gaussian  — truncated normal on ρ
  ga_idx <- which(prop$fam == FAM_GAUSS)
  if (length(ga_idx))
    prop$th1[ga_idx] <- rtnorm_vec(prop$th1[ga_idx], cfg$step_sd,
                                   -0.99, 0.99)
  
  ## BB1 — RW on (λL, λU) then map back
  bb_idx <- which(prop$fam == FAM_BB1)
  if (length(bb_idx)) {
    lam <- t(mapply(bb1_par2tail,
                    prop$th1[bb_idx], prop$th2[bb_idx]))
    lam[,1] <- rtnorm_vec(lam[,1], cfg$step_sd, 1e-3, 0.99)   # λL
    lam[,2] <- rtnorm_vec(lam[,2], cfg$step_sd, 1e-3, 0.99)   # λU
    par <- t(apply(lam, 1, \(x) bb1_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
    prop$th1[bb_idx] <- par[,1]    # θ
    prop$th2[bb_idx] <- par[,2]    # δ
  }
  
  ## ---------- 2. MH ACCEPTANCE --------------------------------------
  vine_prop <- fast_vine_from_particle(prop, skeleton)
  vine_curr <- fast_vine_from_particle(p,    skeleton)
  
  ll_prop <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
                          .Machine$double.eps)))
  ll_curr <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
                          .Machine$double.eps)))
  
  ## proposal symmetric by construction
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

# diagnostic_report <- function(t, tr, U, particles, w_new,
#                               cfg, q_probs = c(0.025, 0.975)) {
#   
#   k_step     <- cfg$k_step
#   M          <- cfg$M
#   print_flag <- (t %% k_step == 0L) || (t == nrow(U))
#   
#   ## ---------- 1. unpack particle state -------------------------------
#   theta_mat <- do.call(rbind, lapply(particles, `[[`, "theta"))
#   rho_mat   <- theta_mat                       # correlations
#   gamma_mat <- do.call(rbind, lapply(particles, `[[`, "gamma"))
#   
#   ## ---------- 2. normalise weights -----------------------------------
#   w_new <- w_new / sum(w_new)
#   ess_t <- 1 / sum(w_new^2)
#   
#   ## ---------- 3. global diagnostics ----------------------------------
#   edge_ct    <- rowSums(gamma_mat)                 # #Gaussian edges / particle
#   mean_edges <- w_mean(edge_ct, w_new)
#   se_edges   <- mc_se(edge_ct, w_new, ess_t)
#   sparsity   <- 1 - mean_edges / cfg$K            # 1 = fully sparse
#   
#   key_vec   <- apply(theta_mat, 1, \(row) paste(row, collapse = ","))
#   n_unique  <- length(unique(key_vec))
#   
#   dists     <- as.matrix(stats::dist(theta_mat))
#   avg_dist  <- mean(dists[lower.tri(dists)])
#   
#   if (print_flag) {
#     cat(sprintf(
#       "t=%4d | tr=%2d | ESS/M=%.3f | dep.edges=%.2f ± %.2f | sparsity=%.3f | unique=%d | L2=%.4f\n",
#       t, tr, ess_t / M, mean_edges, se_edges, sparsity, n_unique, avg_dist))
#   }
#   
#   ## ---------- 4. per-edge summaries ----------------------------------
#   inc_prob <- colSums(gamma_mat * w_new)          # P(γ_e = 1)
#   
#   edge_summ <- lapply(seq_along(inc_prob), function(e) {
#     
#     se_inc <- sqrt(inc_prob[e] * (1 - inc_prob[e]) / ess_t)
#     
#     active  <- gamma_mat[, e] == 1L
#     if (!any(active)) {
#       if (print_flag)
#         cat(sprintf("  Edge %2d : P(dep)=%.3f ± %.3f | never active\n",
#                     e, inc_prob[e], se_inc))
#       return(list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
#                   mu_rho = NA, sd_rho = NA,
#                   q025 = NA, q975 = NA))
#     }
#     
#     rho_e <- rho_mat[active, e]
#     w_e   <- w_new[active] / sum(w_new[active])
#     
#     mu_rho <- sum(w_e * rho_e)
#     sd_rho <- sqrt(sum(w_e * (rho_e - mu_rho)^2))
#     qs     <- w_quantile(rho_e, w_e, q_probs)
#     
#     if (print_flag)
#       cat(sprintf(
#         "  Edge %2d : P(dep)=%.3f ± %.3f | ρ = %.3f (SD %.3f) | 95%% CI = [%.3f, %.3f]\n",
#         e, inc_prob[e], se_inc, mu_rho, sd_rho, qs[1], qs[2]))
#     
#     list(edge = e, p_dep = inc_prob[e], se_inc = se_inc,
#          mu_rho = mu_rho, sd_rho = sd_rho,
#          q025 = qs[1], q975 = qs[2])
#   })
#   
#   edge_df <- do.call(rbind.data.frame, edge_summ)
#   
#   invisible(list(
#     ESS      = ess_t,
#     unique   = n_unique,
#     euc      = avg_dist,
#     sparsity = sparsity,
#     edges    = edge_df
#   ))
# }


diagnostic_report <- function(t, tr, U, particles, w_new,
                              cfg, q_probs = c(0.025, 0.975)) {
  
  k_step     <- cfg$k_step
  M          <- cfg$M
  print_flag <- (t %% k_step == 0L) || (t == nrow(U))
  
  ## ---------- 1. unpack particle state ----------------------------
  fam_mat <- do.call(rbind, lapply(particles, `[[`, "fam"))
  th1_mat <- do.call(rbind, lapply(particles, `[[`, "th1"))
  th2_mat <- do.call(rbind, lapply(particles, `[[`, "th2"))
  
  IND <- 0L; GAU <- 1L; BB1 <- 2L
  
  ## ---------- 2. weights & ESS ------------------------------------
  w_new <- w_new / sum(w_new)
  ess_t <- 1 / sum(w_new^2)
  
  ## ---------- 3. sparsity / uniqueness ----------------------------
  dep_mask   <- fam_mat != IND
  edge_ct    <- rowSums(dep_mask)
  mean_edges <- w_mean(edge_ct, w_new)
  se_edges   <- mc_se(edge_ct, w_new, ess_t)
  sparsity   <- 1 - mean_edges / cfg$K
  
  key_vec  <- apply(th1_mat, 1, \(row) paste(row, collapse = ","))
  n_unique <- length(unique(key_vec))
  
  dists    <- as.matrix(stats::dist(th1_mat))
  avg_dist <- mean(dists[lower.tri(dists)])
  
  ## ---------- 4. console header -----------------------------------
  if (print_flag) {
    cat(sprintf(
      "t=%4d | tr=%2d | ESS/M=%.3f | dep.edges=%.2f ± %.2f | sparsity=%.3f | unique=%d | L2=%.4f\n",
      t, tr, ess_t / M, mean_edges, se_edges, sparsity, n_unique, avg_dist))
  }
  
  ## ---------- 5. family proportions -------------------------------
  fam_codes <- sapply(cfg$families,
                      switch, indep = IND, gaussian = GAU, bb1 = BB1)
  prop_line <- vapply(seq_along(fam_codes), function(j) {
    code <- fam_codes[j]
    p_j  <- w_mean(rowSums(fam_mat == code), w_new) / cfg$K
    sprintf("%s %.2f", cfg$families[j], p_j)
  }, character(1))
  if (print_flag)
    cat("   Family proportions | ", paste(prop_line, collapse = " | "), "\n")
  
  ## ---------- 6. per-edge summaries -------------------------------
  edge_summ <- vector("list", cfg$K)
  
  for (e in seq_len(cfg$K)) {
    
    p_indep <- sum(w_new[fam_mat[, e] == IND])
    p_gauss <- sum(w_new[fam_mat[, e] == GAU])
    p_bb1   <- sum(w_new[fam_mat[, e] == BB1])
    
    res <- list(edge = e,
                p_indep = p_indep,
                p_gauss = p_gauss,
                p_bb1   = p_bb1)
    
    if (print_flag)
      cat(sprintf("     Edge %2d  Indep   : P=%.3f\n", e, p_indep))
    
    ## ---- Gaussian stats ------------------------------------------
    if (p_gauss > 0) {
      mask  <- fam_mat[, e] == GAU
      rho_e <- th1_mat[mask, e]
      w_g   <- w_new[mask] / p_gauss
      mu_rho <- sum(w_g * rho_e)
      sd_rho <- sqrt(sum(w_g * (rho_e - mu_rho)^2))
      qs     <- w_quantile(rho_e, w_g, q_probs)
      
      res <- c(res,
               list(mu_rho = mu_rho, sd_rho = sd_rho,
                    rho_q025 = qs[1], rho_q975 = qs[2]))
      if (print_flag)
        cat(sprintf(
          "                Gaussian: P=%.3f | ρ=%.3f (SD %.3f) CI=[%.3f,%.3f]\n",
          p_gauss, mu_rho, sd_rho, qs[1], qs[2]))
    }
    
    ## ---- BB1 stats -----------------------------------------------
    if (p_bb1 > 0) {
      mask  <- fam_mat[, e] == BB1
      th1_e <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      lam   <- t(mapply(bb1_par2tail, th1_e, th2_e))
      w_b   <- w_new[mask] / p_bb1
      
      stats <- function(x) {
        mu <- sum(w_b * x)
        sd <- sqrt(sum(w_b * (x - mu)^2))
        ci <- w_quantile(x, w_b, q_probs)
        c(mu, sd, ci)
      }
      s_lamL <- stats(lam[,1]); s_lamU <- stats(lam[,2])
      s_th1  <- stats(th1_e);   s_th2  <- stats(th2_e)
      
      res <- c(res,
               list(mu_lambdaL = s_lamL[1], sd_lambdaL = s_lamL[2],
                    lamL_q025  = s_lamL[3], lamL_q975  = s_lamL[4],
                    mu_lambdaU = s_lamU[1], sd_lambdaU = s_lamU[2],
                    lamU_q025  = s_lamU[3], lamU_q975  = s_lamU[4],
                    mu_theta   = s_th1[1],  sd_theta   = s_th1[2],
                    theta_q025 = s_th1[3],  theta_q975 = s_th1[4],
                    mu_delta   = s_th2[1],  sd_delta   = s_th2[2],
                    delta_q025 = s_th2[3],  delta_q975 = s_th2[4]))
      
      if (print_flag)
        cat(sprintf(
          "                BB1     : P=%.3f | λU=%.3f±%.3f CI=[%.3f,%.3f] | λL=%.3f±%.3f CI=[%.3f,%.3f] | θ=%.2f±%.2f | δ=%.2f±%.2f\n",
          p_bb1,
          s_lamU[1], s_lamU[2], s_lamU[3], s_lamU[4],
          s_lamL[1], s_lamL[2], s_lamL[3], s_lamL[4],
          s_th1[1],  s_th1[2],
          s_th2[1],  s_th2[2]))
    }
    
    edge_summ[[e]] <- res
  }
  
  edge_df <- data.table::rbindlist(edge_summ, fill = TRUE)
  
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




