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

FAM_INDEP   <- 0L
FAM_GAUSS   <- 1L
FAM_BB1     <- 2L

T_INDEP  <- bicop_dist("indep")
T_GAUSS  <- bicop_dist("gaussian", parameters = 0)
T_BB1    <- bicop_dist("bb1",      parameters = c(1, 2))  # dummy θ,δ

active_fams <- function(cfg) FAM_INFO[FAM_INFO$name %in% cfg$families, ]

edge_tree_map <- function(d) {
  K   <- d * (d - 1) / 2
  map <- integer(K)
  idx <- 1L
  for (tr in 1:(d - 1)) {          # tree index
    n_edges <- d - tr              # edges in that tree
    map[idx:(idx + n_edges - 1)] <- tr
    idx <- idx + n_edges
  }
  map
}

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


new_particle <- function(cfg) {
  
  K        <- cfg$K
  fam      <- integer(K)
  th1      <- numeric(K)
  th2      <- numeric(K)
  
  for (k in seq_len(K)) {
    
    tr_k <- cfg$edge_tree[k]                 # 1, 2, …
    
    allowed_names <- if (tr_k == 1)
      cfg$families_first else cfg$families_deep
    
    fam_tbl <- FAM_INFO[FAM_INFO$name %in% allowed_names, ]
    weights <- exp(-cfg$lambda * fam_tbl$npar)
    code_k  <- sample(fam_tbl$code, 1, prob = weights)
    fam[k]  <- code_k
    
    if (code_k == FAM_GAUSS) {
      th1[k] <- runif(1, -0.99, 0.99)
      
    } else if (code_k == FAM_BB1) {
      lamU <- rbeta(1, 2, 2); lamL <- rbeta(1, 2, 2)
      tp   <- bb1_tail2par(lamL, lamU)
      th1[k] <- tp[1]; th2[k] <- tp[2]
    }
  }
  
  list(fam = fam, th1 = th1, th2 = th2,
       w = 1 / cfg$M, last_accept = FALSE)
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


rtnorm_vec <- function(mu_vec, sd, a = -0.99, b = 0.99) {
  x <- mu_vec + sd * rnorm(length(mu_vec))
  pmin(pmax(x, a), b)                     # hard clip (reflection is slower)
}


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
      
      tr_idx <- cfg$edge_tree[idx]
      
      allowed_names <- if (tr_idx == 1)
        cfg$families_first else cfg$families_deep
      
      allowed_codes <- FAM_INFO$code[FAM_INFO$name %in% allowed_names]
      
      ## choose a different family *within* the allowed set
      new_code <- sample(allowed_codes[allowed_codes != old_code], 1L)
      
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




diagnostic_report <- function(t, tr, U, particles, w_new,
                              cfg, q_probs = c(0.025, 0.975)) {
  
  k_step     <- cfg$k_step
  M          <- cfg$M
  print_flag <- (t %% k_step == 0L) || (t == nrow(U))
  #print_flag <- FALSE
  
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
  
  list(
    log_pred_density = log_pred_density
  )
}



# -------------------------------------------------------------------------
#  smc_predictive_sample()
#  Draw L iid observations from the posterior predictive C-vine mixture
# -------------------------------------------------------------------------
smc_predictive_sample <- function(particles,          # list of M particles
                                  skel,               # fixed structure
                                  w,                  # numeric(M) normalised
                                  L       = 10000,    # # draws wanted
                                  cache   = FALSE,    # reuse vine objects?
                                  cl      = NULL) {   # existing cluster
  
  stopifnot(abs(sum(w) - 1) < 1e-8)
  M  <- length(particles)
  d  <- ncol(skel$structure)
  
  ## 0️⃣  optional caching -------------------------------------------------
  vines <- if (isTRUE(cache)) {
    lapply(particles, fast_vine_from_particle, skel = skel)
  }
  
  ## 1️⃣  which particle generates each row? ------------------------------
  id <- sample.int(M, L, replace = TRUE, prob = w)
  
  ## 2️⃣  simulation helper (one row) --------------------------------------
  sim_row <- function(j) {               # j = row index (1…L)
    p_idx <- id[j]                       # which particle?
    vine  <- #if (isTRUE(cache)) vines[[p_idx]]
    #else          
      fast_vine_from_particle(particles[[p_idx]], skel)
    rvinecop(1, vine)                    # 1 × d numeric
  }
  
  ## 3️⃣  run in parallel or serial ----------------------------------------
  if (!is.null(cl)) {
    parallel::clusterExport(
      cl,
      c("id", "particles", "vines", "skel", "sim_row"),
      envir = environment())
    
    ## split the indices 1:L as evenly as possible across the workers
    idx_split <- parallel::splitIndices(L, length(cl))
    
    sims <- parallel::parLapply(cl, idx_split, function(chunk) {
      do.call(rbind, lapply(chunk, sim_row))
    })
    out <- do.call(rbind, sims)          # L × d  (rows already in order)
    
  } else {                               # serial fallback
    out <- do.call(rbind, lapply(seq_len(L), sim_row))
  }
  
  dimnames(out) <- NULL
  out                              # matrix(L , d)
}





# ------------------------------------------------------------------
#  update_risk_log()   –– append one row of risk metrics
# ------------------------------------------------------------------
#
#  Arguments
#  ----------
#  risk_log   data.table | initially empty (0 rows) and kept by caller
#  U_pred     L × d matrix of uniform draws from smc_predictive_sample()
#  mu_fc      numeric-vector length d   –– conditional means  (AR part)
#  sig_fc     numeric-vector length d   –– conditional stdevs (GARCH part)
#  t_idx      scalar, time index (stored in the log)
#  alphas     tail probabilities for VaR / ES   (default 95 % & 97.5 %)
#
#  Returns
#  -------
#  Same data.table, one extra row per call
# ------------------------------------------------------------------
update_risk_log <- function(risk_log,
                            U_pred,
                            mu_f,
                            sig_f,
                            t_idx,
                            alphas = c(0.95, 0.975))
{
  stopifnot(is.matrix(U_pred),
            length(mu_f)  == ncol(U_pred),
            length(sig_f) == ncol(U_pred))
  
  
  sig_f = as.numeric(sig_f)
  mu_f = as.numeric(mu_f)
  
  ## 1. quantile transform & rescale -------------------------------
  Z      <- qnorm(U_pred)                           # L × d
  R_pred <- sweep(Z,  2, sig_f, `*`)
  R_pred <- sweep(R_pred, 2, mu_f,  `+`)
  
  L <- nrow(R_pred);  d <- ncol(R_pred)
  series_names <- paste0("S", seq_len(d))           # generic column labels
  dimnames(R_pred) <- list(NULL, series_names)
  
  ## 2. point moments ----------------------------------------------
  mu_hat  <- colMeans(R_pred)
  var_hat <- colMeans((R_pred - rep(mu_hat, each = L))^2)
  
  ## 3. tail risk ---------------------------------------------------
  losses <- -R_pred                          # L × d matrix
  
  VaR <- vapply(alphas,                          # (K × d)
                function(a)
                  apply(losses, 2, quantile, probs = a),
                numeric(d))
  
  ES  <- vapply(seq_along(alphas),   # k = 1 … K
                function(k) {
                  thr <- VaR[ , k]   # one threshold per series
                  vapply(seq_len(d), function(j) {
                    lj <- losses[ , j]
                    exc <- lj[lj >= thr[j]]
                    if (length(exc)) mean(exc) else NA_real_
                  }, numeric(1))
                },
                numeric(d))
  dimnames(VaR) <- dimnames(ES) <- list(colnames(losses), paste0("a", alphas))
  
  ## 4. 95 % predictive interval -----------------------------------
  CI_lower <- apply(R_pred, 2, quantile, probs = 0.025)
  CI_upper <- apply(R_pred, 2, quantile, probs = 0.975)
  
  ## 5. gather as one long row (wide format) ------------------------
  row <- data.table::data.table(t_idx = t_idx)
  
  ## names of the d series, e.g.  c("S1","S2", …)
  series_names <- colnames(R_pred)
  
  ## row is the single-row data.table that accumulates the statistics
  ## ----------------------------------------------------------------
  add_vec <- function(x, prefix) {
    cols <- paste0(prefix, "_", series_names)   # e.g.  mean_S1, mean_S2 …
    stopifnot(length(x) == length(cols))        # sanity check
    row[, (cols) := as.list(unname(x))]         # <- list = one entry / col
  }
  
  ## point estimates ------------------------------------------------
  add_vec(mu_hat,  "mean")
  add_vec(var_hat, "var")
  add_vec(CI_lower,"ci_lo")
  add_vec(CI_upper,"ci_hi")
  
  ## tail–risk measures ---------------------------------------------
  for (k in seq_along(alphas)) {
    add_vec(VaR[ , k], sprintf("VaR%g", alphas[k] * 100))
    add_vec(ES [ , k], sprintf("ES%g",  alphas[k] * 100))
  }
  
  ## 6. append & return --------------------------------------------
  data.table::rbindlist(list(risk_log, row), use.names = TRUE, fill = TRUE)
}



























