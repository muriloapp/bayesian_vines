library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(matrixStats)


FAM_INFO <- data.frame(
  # name: how you refer to it in cfg$families / families_first / families_deep
  # rv_name: rvinecopulib family string
  # rotation: 0 or 180
  # code: internal integer id you sample/compare on
  # npar: number of free parameters (used by sparsity penalty)
  name      = c("indep", "gaussian", "bb1", "bb1r180", "bb8r180"),
  rv_name   = c("indep","gaussian","bb1","bb1","bb8"),
  rotation  = c(0L, 0L, 0L, 180L, 180L),
  code      = c(0L, 1L, 3L, 13L, 16L),
  npar      = c(0L, 1L, 2L, 2L, 2L),
  stringsAsFactors = FALSE
)

## Shortcuts (optional)
FAM_INDEP    <- FAM_INFO$code[FAM_INFO$name == "indep"]
FAM_GAUSS    <- FAM_INFO$code[FAM_INFO$name == "gaussian"]
FAM_BB1      <- FAM_INFO$code[FAM_INFO$name == "bb1"]
FAM_BB1R180  <- FAM_INFO$code[FAM_INFO$name == "bb1r180"]
FAM_BB8R180  <- FAM_INFO$code[FAM_INFO$name == "bb8r180"]

## Minimal templates (we’ll build bicops from rv_name/rotation on the fly too)
T_INDEP     <- bicop_dist("indep")
T_GAUSS     <- bicop_dist("gaussian", parameters = 0)

T_BB1       <- bicop_dist("bb1", rotation = 0,   parameters = c(1, 2))
T_BB1R180   <- bicop_dist("bb1", rotation = 180, parameters = c(1, 2))

# T_BB8R180   <- bicop_dist("bb8", rotation = 180, parameters = c(1, 2))  # dummy two-par



## Priors

# log_prior_edge <- function(fam, th1, th2, cfg) {
#   
#   if (fam == FAM_INDEP) return(0)           # constant
#   
#   if (fam == FAM_GAUSS) {
#     lp <- log(1/1.98)                       # Uniform ρ
#     return(lp - cfg$lambda)                 # + penalty
#     
#   } else {  # BB1
#     lam <- bb1_par2tail(th1, th2)
#     lp_tail <- dbeta(lam[1], 2, 2, log=TRUE) +
#       dbeta(lam[2], 2, 2, log=TRUE)
#     lp <- lp_tail + bb1_log_jacobian(lam[1], lam[2]) - 2*cfg$lambda
#     return(fillna_neg(lp))
#   }
# }

log_prior_edge <- function(fam, th1, th2, cfg) {
  if (fam == FAM_INDEP) return(0)
  
  if (fam == FAM_GAUSS) {
    # Uniform rho in (-0.99,0.99) approx -> constant; keep your complexity penalty
    return(log(1/1.98) - cfg$lambda)
  }
  
  if (fam == FAM_BB1) {
    lam <- bb1_par2tail(th1, th2)                                  # (λL, λU)
    lp_tail <- dbeta(lam[1], 2, 2, log=TRUE) + dbeta(lam[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb1_log_jacobian(lam[1], lam[2]) - 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  if (fam == FAM_BB1R180) {
    lam_rot <- bb1r180_par2tail(th1, th2)                           # (λL_rot, λU_rot)
    lp_tail <- dbeta(lam_rot[1], 2, 2, log=TRUE) + dbeta(lam_rot[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb1r180_log_jacobian(lam_rot[1], lam_rot[2]) - 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  if (fam == FAM_BB8R180) {
    # Simple complexity prior only (like Gaussian); adjust if you later design a tail prior
    return(- 2 * cfg$lambda)
  }
  
  stop("Unknown family code in log_prior_edge: ", fam)
}


log_prior <- function(p, cfg) {
  sum(mapply(log_prior_edge,
             p$fam, p$th1, p$th2,
             MoreArgs = list(cfg = cfg)))
}

## Intialize particles

# new_particle <- function(cfg) {
#   
#   K        <- cfg$K
#   fam      <- integer(K)
#   th1      <- numeric(K)
#   th2      <- numeric(K)
#   
#   for (k in seq_len(K)) {
#     
#     tr_k <- cfg$edge_tree[k]                 # 1, 2, …
#     
#     allowed_names <- if (tr_k == 1)
#       cfg$families_first else cfg$families_deep
#     
#     fam_tbl <- FAM_INFO[FAM_INFO$name %in% allowed_names, ]
#     weights <- exp(-cfg$lambda * fam_tbl$npar)
#     code_k  <- sample(fam_tbl$code, 1, prob = weights)
#     fam[k]  <- code_k
#     
#     if (code_k == FAM_GAUSS) {
#       th1[k] <- runif(1, -0.99, 0.99)
#       
#     } else if (code_k == FAM_BB1) {
#       lamU <- rbeta(1, 2, 2); lamL <- rbeta(1, 2, 2)
#       tp   <- bb1_tail2par(lamL, lamU)
#       pars <- sanitize_bb1(tp[1], tp[2])
#       th1[k] <- pars[1]; th2[k] <- pars[2]
#     }
#   }
#   
#   list(fam = fam, th1 = th1, th2 = th2,
#        w = 1 / cfg$M, last_accept = FALSE)
# }

new_particle <- function(cfg) {
  K   <- cfg$K
  fam <- integer(K); th1 <- numeric(K); th2 <- numeric(K)
  
  for (k in seq_len(K)) {
    tr_k <- cfg$edge_tree[k]
    allowed_names <- if (tr_k == 1) cfg$families_first else cfg$families_deep
    fam_tbl <- FAM_INFO[FAM_INFO$name %in% allowed_names, ]
    #weights <- exp(-cfg$lambda * fam_tbl$npar)
    weights <- rep(1/dim(fam_tbl)[1], dim(fam_tbl)[1])
    code_k  <- sample(fam_tbl$code, 1, prob = weights)
    fam[k]  <- code_k
    
    if (code_k == FAM_INDEP) {
      th1[k] <- 0; th2[k] <- 0
      
    } else if (code_k == FAM_GAUSS) {
      th1[k] <- runif(1, -0.99, 0.99); th2[k] <- 0
      
    } else if (code_k == FAM_BB1) {
      lamU <- rbeta(1, 2, 2); lamL <- rbeta(1, 2, 2)
      pars <- sanitize_bb1(bb1_tail2par(lamL, lamU)[1],
                           bb1_tail2par(lamL, lamU)[2])
      th1[k] <- pars[1]; th2[k] <- pars[2]
      
    } else if (code_k == FAM_BB1R180) {
      # sample in rotated tail space, then map with rotated converter
      lamL_rot <- rbeta(1, 2, 2); lamU_rot <- rbeta(1, 2, 2)
      pars <- sanitize_bb1(bb1r180_tail2par(lamL_rot, lamU_rot)[1],
                           bb1r180_tail2par(lamL_rot, lamU_rot)[2])
      th1[k] <- pars[1]; th2[k] <- pars[2]
      
    } else if (code_k == FAM_BB8R180) {
      # lightweight init: random but sane
      th1[k] <- runif(1, 0.1, 3.0)  # θ rough box (adjust later if you like)
      th2[k] <- runif(1, 1.1, 5.0)  # δ rough box
    }
  }
  
  list(fam = fam, th1 = th1, th2 = th2, w = 1 / cfg$M, last_accept = FALSE)
}


## Vine for each particle
# 
# fast_vine_from_particle <- function(p, skel) {
#   
#   pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)
#   idx <- 1L
#   
#   for (tr in seq_along(pcs)) {
#     for (ed in seq_along(pcs[[tr]])) {
#       
#       f <- p$fam[idx]
#       
#       if (f == FAM_INDEP) {
#         pcs[[tr]][[ed]] <- T_INDEP           # shallow copy by assignment
#         
#       } else if (f == FAM_GAUSS) {
#         bc <- rlang::duplicate(T_GAUSS, shallow = TRUE)
#         bc$parameters[1] <- p$th1[idx]       # update ρ
#         pcs[[tr]][[ed]] <- bc
#         
#       } else {                               # BB1
#         pars <- sanitize_bb1(p$th1[idx], p$th2[idx])   # ⟵ NEW
#         p$th1[idx] <- pars[1]                          # keep particle in sync
#         p$th2[idx] <- pars[2]
#         
#         bc <- rlang::duplicate(T_BB1, shallow = TRUE)
#         bc$parameters[1:2] <- pars
#         pcs[[tr]][[ed]] <- bc
#       }
#       idx <- idx + 1L
#     }
#   }
#   
#   vinecop_dist(pcs, skel$structure)
# }

fast_vine_from_particle <- function(p, skel) {
  pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)
  idx <- 1L
  
  for (tr in seq_along(pcs)) {
    for (ed in seq_along(pcs[[tr]])) {
      code <- p$fam[idx]
      spec <- fam_spec(code)
      
      if (code == FAM_INDEP) {
        pcs[[tr]][[ed]] <- T_INDEP
        
      } else if (code == FAM_GAUSS) {
        bc <- rlang::duplicate(T_GAUSS, shallow = TRUE)
        bc$parameters[1] <- p$th1[idx]
        pcs[[tr]][[ed]] <- bc
        
      } else if (code %in% c(FAM_BB1, FAM_BB1R180)) {
        pars <- sanitize_bb1(p$th1[idx], p$th2[idx])
        p$th1[idx] <- pars[1]; p$th2[idx] <- pars[2]
        
        tmpl <- if (code == FAM_BB1) T_BB1 else T_BB1R180
        bc <- rlang::duplicate(tmpl, shallow = TRUE)
        bc$parameters[1:2] <- pars
        pcs[[tr]][[ed]] <- bc
        
      } else if (code == FAM_BB8R180) {
        bc <- rlang::duplicate(T_BB8R180, shallow = TRUE)
        bc$parameters[1:2] <- c(p$th1[idx], p$th2[idx])
        pcs[[tr]][[ed]] <- bc
        
      } else {
        stop("fast_vine_from_particle: unsupported family code: ", code)
      }
      
      idx <- idx + 1L
    }
  }
  
  vinecop_dist(pcs, skel$structure)
}


compute_log_incr <- function(particles, u_row, skeleton, cfg) {
  log_incr <- numeric(cfg$M)
  
  for (j in seq_along(particles)) {
    p <- particles[[j]]
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



mh_step <- function(p, data_up_to_t, skeleton, cfg) {
  K    <- cfg$K
  prop <- p
  
  ## (a) Family flips — prior-only init (no MLE, no data-length gating)
  end_first_tr <- sum(cfg$edge_tree == 1L)  # #edges in Tree 1 (by your ordering)
  if (end_first_tr > 0L) {
    flip_mask <- runif(end_first_tr) < cfg$q_flip
    if (any(flip_mask)) {
      idxs <- which(flip_mask)
      for (idx in idxs) {
        tr_idx <- cfg$edge_tree[idx]
        if (tr_idx != 1L) next  # keep flips on Tree 1 only
        
        old_code <- prop$fam[idx]
        allowed_names <- if (tr_idx == 1L) cfg$families_first else cfg$families_deep
        allowed_codes <- FAM_INFO$code[FAM_INFO$name %in% allowed_names]
        allowed_codes <- setdiff(allowed_codes, old_code)  # avoid no-op flip
        if (!length(allowed_codes)) next
        
        new_code <- sample(allowed_codes, 1L)
        prop$fam[idx] <- new_code
        
        ## Initialize from priors (no MLE)
        if (new_code == FAM_INDEP) {
          prop$th1[idx] <- 0;                       prop$th2[idx] <- 0
          
        } else if (new_code == FAM_GAUSS) {
          prop$th1[idx] <- runif(1, -0.99, 0.99);   prop$th2[idx] <- 0
          
        } else if (new_code == FAM_BB1) {
          lamL <- rbeta(1, 2, 2); lamU <- rbeta(1, 2, 2)
          tp   <- bb1_tail2par(lamL, lamU)
          tp   <- sanitize_bb1(tp[1], tp[2])
          prop$th1[idx] <- tp[1];                  prop$th2[idx] <- tp[2]
          
        } else if (new_code == FAM_BB1R180) {
          lamLr <- rbeta(1, 2, 2); lamUr <- rbeta(1, 2, 2)
          tp    <- bb1r180_tail2par(lamLr, lamUr)
          tp    <- sanitize_bb1(tp[1], tp[2])
          prop$th1[idx] <- tp[1];                  prop$th2[idx] <- tp[2]
          
        } else if (new_code == FAM_BB8R180) {
          # diffuse box prior (adjust if you add a tailored prior later)
          prop$th1[idx] <- runif(1, 0.1, 3.0)
          prop$th2[idx] <- runif(1, 1.1, 5.0)
        }
      }
    }
  }
  
  ## (b) Random-walk on parameters of current families
  # Gaussian — truncated normal on ρ
  ga_idx <- which(prop$fam == FAM_GAUSS)
  if (length(ga_idx))
    prop$th1[ga_idx] <- rtnorm_vec(prop$th1[ga_idx], cfg$step_sd, -0.99, 0.99)
  
  # BB1 (0°) — RW in (λL, λU), then map back and sanitize
  bb1_idx <- which(prop$fam == FAM_BB1)
  if (length(bb1_idx)) {
    lam <- t(mapply(bb1_par2tail, prop$th1[bb1_idx], prop$th2[bb1_idx]))
    lam[, 1] <- rtnorm_vec(lam[, 1], cfg$step_sd, 1e-3, 0.99)
    lam[, 2] <- rtnorm_vec(lam[, 2], cfg$step_sd, 1e-3, 0.99)
    par <- t(apply(lam, 1, \(x) bb1_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
    prop$th1[bb1_idx] <- par[, 1]; prop$th2[bb1_idx] <- par[, 2]
  }
  
  # BB1^180 — RW in (λL_rot, λU_rot), map back with rotated converter
  bb1r_idx <- which(prop$fam == FAM_BB1R180)
  if (length(bb1r_idx)) {
    lamr <- t(mapply(bb1r180_par2tail, prop$th1[bb1r_idx], prop$th2[bb1r_idx]))
    lamr[, 1] <- rtnorm_vec(lamr[, 1], cfg$step_sd, 1e-3, 0.99)
    lamr[, 2] <- rtnorm_vec(lamr[, 2], cfg$step_sd, 1e-3, 0.99)
    par <- t(apply(lamr, 1, \(x) bb1r180_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
    prop$th1[bb1r_idx] <- par[, 1]; prop$th2[bb1r_idx] <- par[, 2]
  }
  
  # BB8^180 — simple box RW in parameter space
  bb8r_idx <- which(prop$fam == FAM_BB8R180)
  if (length(bb8r_idx)) {
    prop$th1[bb8r_idx] <- pmin(pmax(prop$th1[bb8r_idx] + rnorm(length(bb8r_idx), 0, cfg$step_sd), 0.05), 6.0)
    prop$th2[bb8r_idx] <- pmin(pmax(prop$th2[bb8r_idx] + rnorm(length(bb8r_idx), 0, cfg$step_sd), 1.05), 7.0)
  }
  
  ## (c) MH ratio
  vine_prop <- fast_vine_from_particle(prop, skeleton)
  vine_curr <- fast_vine_from_particle(p,    skeleton)
  
  ll_prop <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_prop), .Machine$double.eps))))
  ll_curr <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_curr), .Machine$double.eps))))
  
  log_acc <- (ll_prop + log_prior(prop, cfg)) - (ll_curr + log_prior(p, cfg))
  
  if (log(runif(1)) < log_acc) {
    prop$last_accept <- TRUE
    return(prop)
  } else {
    p$last_accept <- FALSE
    return(p)
  }
}



## Adaptive RW step

compute_adapt_step_sd <- function(cfg, acc_pct, lambda = 0.25, target_acc = 0.30, sd_min = 0.01, sd_max=0.05){
  log_sd_new <- log(cfg$step_sd) + lambda * (acc_pct/100 - target_acc)
  step_sd <- pmin(pmax(exp(log_sd_new), sd_min), sd_max)
  cat(sprintf(
    "step_sd = %4f \n\n",
    step_sd
  ))
  return(step_sd)
}

## Resample move

mh_worker_standard <- function(idx, n_mh, slice, data_up_to_t, skeleton, cfg)
{
  p         <- slice[[idx]]
  acc_local <- 0L
  
  for (k in seq_len(n_mh)) {
    p <- mh_step(p, data_up_to_t, skeleton, cfg)
    if (isTRUE(p$last_accept)) acc_local <- acc_local + 1L
  }
  
  list(p = p, acc = acc_local)
}

mh_worker_block <- function(idx, n_mh, slice, data_up_to_t, temp_skel, tr, cfg)
{
  p         <- slice[[idx]]
  acc_local <- 0L
  
  for (k in seq_len(n_mh)) {
    p <- mh_step_in_tree(p, tr, data_up_to_t, temp_skel, cfg)
    if (isTRUE(p$last_accept)) acc_local <- acc_local + 1L
  }
  
  list(p = p, acc = acc_local)
}


resample_move <- function(particles, newAncestors, data_up_to_t,
                          cl, type = "standard",
                          cfg, skeleton = NULL, tr = NULL, temp_skel = NULL)
{
  ## --- 1a.  Resample ------------------------------------------------
  particles <- particles[newAncestors]        # simple copy
  for (p in particles) p$w <- 1 / cfg$M       # reset weights
  
  ## --- 1b.  Ship particles ONCE to each worker ----------------------
  clusterExport(cl, "particles")
  
  ## --- 1c.  Parallel MH moves --------------------------------------
  tic("MH move")
  mh_results <- switch(type,
                       standard = parLapply(
                         cl, seq_along(particles),
                         mh_worker_standard,
                         n_mh         = cfg$n_mh,
                         slice        = particles,
                         data_up_to_t = data_up_to_t,
                         skeleton     = skeleton,
                         cfg          = cfg
                       ),
                       block = parLapply(
                         cl, seq_along(particles),
                         mh_worker_block,
                         n_mh         = cfg$n_mh,
                         slice        = particles,
                         data_up_to_t = data_up_to_t,
                         temp_skel    = temp_skel,
                         tr           = tr,
                         cfg          = cfg
                       )
  )
  toc()
  
  ## --- 1d.  Collect outputs ----------------------------------------
  particles   <- lapply(mh_results, `[[`, "p")
  acc_vec     <- vapply(mh_results, `[[`, integer(1), "acc")
  mh_n_prop   <- cfg$M * cfg$n_mh
  mh_n_acc    <- sum(acc_vec)
  acc_pct     <- 100 * mh_n_acc / mh_n_prop
  
  cat(sprintf(
    "MH acceptance = %4d / %4d = %.2f%%\n\n",
    mh_n_acc, mh_n_prop, acc_pct
  ))
  
  list(particles = particles, acc_pct = acc_pct)
}


resample_move_old <- function(particles, newAncestors, data_up_to_t, cl, type, cfg, skeleton=NULL, tr=NULL, temp_skel=NULL) {

  #idx       <- sample.int(M, M, TRUE, prob = w_new)
  particles <- particles[newAncestors]                  # cópia simples
  for (p in particles) {                       # reinicia pesos
    p$w <- 1 / cfg$M
  }

  mh_n_prop <- cfg$M * cfg$n_mh
  mh_n_acc  <- 0
  #clusterSetRNGStream(cl, 43)

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

  mh_n_acc <- sum(vapply(mh_results, `[[`, integer(1), "acc"))
  particles <- lapply(mh_results, `[[`, "p")

  acc_pct <- 100 * mh_n_acc / mh_n_prop
  cat(sprintf(
    "MH acceptance = %4d / %4d  =  %.2f%%\n\n",
    mh_n_acc, mh_n_prop, acc_pct
  ))

  return(list(particles = particles, acc_pct = acc_pct))
}

resample_move_serialized <- function(particles,
                              newAncestors,
                              data_up_to_t,
                              type,           
                              cfg,
                              skeleton = NULL,
                              tr = NULL,
                              temp_skel = NULL) {
  
  particles <- particles[newAncestors]
  for (p in particles) p$w <- 1 / cfg$M
  
  mh_n_prop <- cfg$M * cfg$n_mh
  mh_results <- vector("list", length(particles))   
  
  for (i in seq_along(particles)) {
    if (i==220){break}
    p    <- particles[[i]]
    acc  <- 0L
    
    for (k in seq_len(cfg$n_mh)) {
      if (type == "standard") {
        p <- mh_step(p, data_up_to_t, skeleton, cfg)
      } else if (type == "block") {
        p <- mh_step_in_tree(p, tr, data_up_to_t, temp_skel, cfg)
      } else {
        stop(sprintf("Unknown type '%s' (use 'standard' or 'block')", type))
      }
      if (isTRUE(p$last_accept)) acc <- acc + 1L
    }
        mh_results[[i]] <- list(p = p, acc = acc)
  }
  
  mh_n_acc <- sum(vapply(mh_results, `[[`, integer(1), "acc"))
  particles <- lapply(mh_results, `[[`, "p")
  
  acc_pct <- 100 * mh_n_acc / mh_n_prop
  cat(sprintf("MH acceptance = %4d / %4d  =  %.2f%%\n\n",
              mh_n_acc, mh_n_prop, acc_pct))
  
  list(particles = particles, acc_pct = acc_pct)
}


# ## Diagnostics
# 
# diagnostic_report <- function(t, tr, U, particles, w_new,
#                               cfg, q_probs = c(0.025, 0.975)) {
#   
#   k_step     <- cfg$k_step
#   M          <- cfg$M
#   print_flag <- (t %% k_step == 0L) || (t == nrow(U))
#   #print_flag <- FALSE
#   
#   fam_mat <- do.call(rbind, lapply(particles, `[[`, "fam"))
#   th1_mat <- do.call(rbind, lapply(particles, `[[`, "th1"))
#   th2_mat <- do.call(rbind, lapply(particles, `[[`, "th2"))
#   
#   IND  <- FAM_INDEP
#   GAU  <- FAM_GAUSS
#   BB1c <- FAM_BB1
#   BB1s <- FAM_BB1R180
#   BB8s <- FAM_BB8R180
#   
#   w_new <- w_new / sum(w_new)
#   ess_t <- 1 / sum(w_new^2)
#   
#   dep_mask   <- fam_mat != IND
#   edge_ct    <- rowSums(dep_mask)
#   mean_edges <- w_mean(edge_ct, w_new)
#   se_edges   <- mc_se(edge_ct, w_new, ess_t)
#   sparsity   <- 1 - mean_edges / cfg$K
#   
#   key_vec  <- apply(th1_mat, 1, \(row) paste(row, collapse = ","))
#   n_unique <- length(unique(key_vec))
#   
#   dists    <- as.matrix(stats::dist(th1_mat))
#   avg_dist <- mean(dists[lower.tri(dists)])
#   
#   if (print_flag) {
#     cat(sprintf(
#       "t=%4d | tr=%2d | ESS/M=%.3f | dep.edges=%.2f ± %.2f | sparsity=%.3f | unique=%d | L2=%.4f\n",
#       t, tr, ess_t / M, mean_edges, se_edges, sparsity, n_unique, avg_dist))
#   }
#   
#   # family proportions 
#   fam_codes <- sapply(cfg$families,
#                       switch, indep = IND, gaussian = GAU, bb1 = BB1)
#   prop_line <- vapply(seq_along(fam_codes), function(j) {
#     code <- fam_codes[j]
#     p_j  <- w_mean(rowSums(fam_mat == code), w_new) / cfg$K
#     sprintf("%s %.2f", cfg$families[j], p_j)
#   }, character(1))
#   if (print_flag)
#     cat("   Family proportions | ", paste(prop_line, collapse = " | "), "\n")
#   
#   # per-edge summaries 
#   edge_summ <- vector("list", cfg$K)
#   
#   for (e in seq_len(cfg$K)) {
#     
#     p_indep <- sum(w_new[fam_mat[, e] == IND])
#     p_gauss <- sum(w_new[fam_mat[, e] == GAU])
#     p_bb1   <- sum(w_new[fam_mat[, e] == BB1])
#     
#     res <- list(edge = e,
#                 p_indep = p_indep,
#                 p_gauss = p_gauss,
#                 p_bb1   = p_bb1)
#     
#     if (print_flag)
#       cat(sprintf("     Edge %2d  Indep   : P=%.3f\n", e, p_indep))
#     
#     # Gaussian stats 
#     if (p_gauss > 0) {
#       mask  <- fam_mat[, e] == GAU
#       rho_e <- th1_mat[mask, e]
#       w_g   <- w_new[mask] / p_gauss
#       mu_rho <- sum(w_g * rho_e)
#       sd_rho <- sqrt(sum(w_g * (rho_e - mu_rho)^2))
#       qs     <- w_quantile(rho_e, w_g, q_probs)
#       
#       res <- c(res,
#                list(mu_rho = mu_rho, sd_rho = sd_rho,
#                     rho_q025 = qs[1], rho_q975 = qs[2]))
#       if (print_flag)
#         cat(sprintf(
#           "                Gaussian: P=%.3f | ρ=%.3f (SD %.3f) CI=[%.3f,%.3f]\n",
#           p_gauss, mu_rho, sd_rho, qs[1], qs[2]))
#     }
#     
#     # BB1 stats 
#     if (p_bb1 > 0) {
#       mask  <- fam_mat[, e] == BB1
#       th1_e <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
#       lam   <- t(mapply(bb1_par2tail, th1_e, th2_e))
#       w_b   <- w_new[mask] / p_bb1
#       
#       stats <- function(x) {
#         mu <- sum(w_b * x)
#         sd <- sqrt(sum(w_b * (x - mu)^2))
#         ci <- w_quantile(x, w_b, q_probs)
#         c(mu, sd, ci)
#       }
#       s_lamL <- stats(lam[,1]); s_lamU <- stats(lam[,2])
#       s_th1  <- stats(th1_e);   s_th2  <- stats(th2_e)
#       
#       res <- c(res,
#                list(mu_lambdaL = s_lamL[1], sd_lambdaL = s_lamL[2],
#                     lamL_q025  = s_lamL[3], lamL_q975  = s_lamL[4],
#                     mu_lambdaU = s_lamU[1], sd_lambdaU = s_lamU[2],
#                     lamU_q025  = s_lamU[3], lamU_q975  = s_lamU[4],
#                     mu_theta   = s_th1[1],  sd_theta   = s_th1[2],
#                     theta_q025 = s_th1[3],  theta_q975 = s_th1[4],
#                     mu_delta   = s_th2[1],  sd_delta   = s_th2[2],
#                     delta_q025 = s_th2[3],  delta_q975 = s_th2[4]))
#       
#       if (print_flag)
#         cat(sprintf(
#           "                BB1     : P=%.3f | λU=%.3f±%.3f CI=[%.3f,%.3f] | λL=%.3f±%.3f CI=[%.3f,%.3f] | θ=%.2f±%.2f | δ=%.2f±%.2f\n",
#           p_bb1,
#           s_lamU[1], s_lamU[2], s_lamU[3], s_lamU[4],
#           s_lamL[1], s_lamL[2], s_lamL[3], s_lamL[4],
#           s_th1[1],  s_th1[2],
#           s_th2[1],  s_th2[2]))
#     }
#     
#     edge_summ[[e]] <- res
#   }
#   
#   edge_df <- data.table::rbindlist(edge_summ, fill = TRUE)
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
  # print_flag <- FALSE
  
  fam_mat <- do.call(rbind, lapply(particles, `[[`, "fam"))
  th1_mat <- do.call(rbind, lapply(particles, `[[`, "th1"))
  th2_mat <- do.call(rbind, lapply(particles, `[[`, "th2"))
  
  ## family codes
  IND  <- FAM_INDEP
  GAU  <- FAM_GAUSS
  BB1c <- FAM_BB1
  BB1s <- FAM_BB1R180
  BB8s <- FAM_BB8R180
  
  ## weights & ESS
  w_new <- w_new / sum(w_new)
  ess_t <- 1 / sum(w_new^2)
  
  ## sparsity etc.
  dep_mask   <- fam_mat != IND
  edge_ct    <- rowSums(dep_mask)
  mean_edges <- w_mean(edge_ct, w_new)
  se_edges   <- mc_se(edge_ct, w_new, ess_t)
  sparsity   <- 1 - mean_edges / cfg$K
  
  key_vec  <- apply(th1_mat, 1, \(row) paste(row, collapse = ","))
  n_unique <- length(unique(key_vec))
  
  dists    <- as.matrix(stats::dist(th1_mat))
  avg_dist <- mean(dists[lower.tri(dists)])
  
  if (print_flag) {
    cat(sprintf(
      "t=%4d | tr=%2d | ESS/M=%.3f | dep.edges=%.2f ± %.2f | sparsity=%.3f | unique=%d | L2=%.4f\n",
      t, tr, ess_t / M, mean_edges, se_edges, sparsity, n_unique, avg_dist))
  }
  
  ## family proportions (driven by cfg$families and FAM_INFO)
  fam_codes <- FAM_INFO$code[match(cfg$families, FAM_INFO$name)]
  prop_line <- vapply(seq_along(fam_codes), function(j) {
    code <- fam_codes[j]
    if (is.na(code)) {
      sprintf("%s n/a", cfg$families[j])
    } else {
      p_j <- w_mean(rowSums(fam_mat == code), w_new) / cfg$K
      sprintf("%s %.2f", cfg$families[j], p_j)
    }
  }, character(1))
  if (print_flag)
    cat("   Family proportions | ", paste(prop_line, collapse = " | "), "\n")
  
  ## per-edge summaries
  edge_summ <- vector("list", cfg$K)
  
  for (e in seq_len(cfg$K)) {
    
    p_indep <- sum(w_new[fam_mat[, e] == IND ])
    p_gauss <- sum(w_new[fam_mat[, e] == GAU ])
    p_bb1   <- sum(w_new[fam_mat[, e] == BB1c])
    p_bb1r  <- sum(w_new[fam_mat[, e] == BB1s])
    p_bb8r  <- sum(w_new[fam_mat[, e] == BB8s])
    
    res <- list(edge   = e,
                p_indep = p_indep,
                p_gauss = p_gauss,
                p_bb1   = p_bb1,
                p_bb1r  = p_bb1r,
                p_bb8r  = p_bb8r)
    
    if (print_flag)
      cat(sprintf("     Edge %2d  Indep   : P=%.3f\n", e, p_indep))
    
    ## Gaussian stats
    if (p_gauss > 0) {
      mask   <- fam_mat[, e] == GAU
      rho_e  <- th1_mat[mask, e]
      w_g    <- w_new[mask] / p_gauss
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
    
    ## BB1 (0°) stats (tails + params)
    if (p_bb1 > 0) {
      mask   <- fam_mat[, e] == BB1c
      th1_e  <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      lam    <- t(mapply(bb1_par2tail, th1_e, th2_e))   # (λL, λU)
      w_b    <- w_new[mask] / p_bb1
      
      stats <- function(x) {
        mu <- sum(w_b * x)
        sd <- sqrt(sum(w_b * (x - mu)^2))
        ci <- w_quantile(x, w_b, q_probs)
        c(mu, sd, ci)
      }
      s_lamL <- stats(lam[, 1]); s_lamU <- stats(lam[, 2])
      s_th1  <- stats(th1_e);    s_th2  <- stats(th2_e)
      
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
    
    ## BB1^180 (survival) stats – rotated tail map
    if (p_bb1r > 0) {
      mask    <- fam_mat[, e] == BB1s
      th1_e   <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      lam_rot <- t(mapply(bb1r180_par2tail, th1_e, th2_e))  # (λL_rot, λU_rot)
      w_br    <- w_new[mask] / p_bb1r
      
      stats <- function(x) {
        mu <- sum(w_br * x)
        sd <- sqrt(sum(w_br * (x - mu)^2))
        ci <- w_quantile(x, w_br, q_probs)
        c(mu, sd, ci)
      }
      s_lamLr <- stats(lam_rot[, 1]); s_lamUr <- stats(lam_rot[, 2])
      s_th1r  <- stats(th1_e);        s_th2r  <- stats(th2_e)
      
      res <- c(res,
               list(mu_lambdaL_rot = s_lamLr[1], sd_lambdaL_rot = s_lamLr[2],
                    lamL_rot_q025  = s_lamLr[3], lamL_rot_q975  = s_lamLr[4],
                    mu_lambdaU_rot = s_lamUr[1], sd_lambdaU_rot = s_lamUr[2],
                    lamU_rot_q025  = s_lamUr[3], lamU_rot_q975  = s_lamUr[4],
                    mu_theta_rot   = s_th1r[1],  sd_theta_rot   = s_th1r[2],
                    theta_rot_q025 = s_th1r[3],  theta_rot_q975 = s_th1r[4],
                    mu_delta_rot   = s_th2r[1],  sd_delta_rot   = s_th2r[2],
                    delta_rot_q025 = s_th2r[3],  delta_rot_q975 = s_th2r[4]))
      
      if (print_flag)
        cat(sprintf(
          "                BB1^180 : P=%.3f | λU(rot)=%.3f±%.3f CI=[%.3f,%.3f] | λL(rot)=%.3f±%.3f CI=[%.3f,%.3f] | θ=%.2f±%.2f | δ=%.2f±%.2f\n",
          p_bb1r,
          s_lamUr[1], s_lamUr[2], s_lamUr[3], s_lamUr[4],
          s_lamLr[1], s_lamLr[2], s_lamLr[3], s_lamLr[4],
          s_th1r[1],  s_th1r[2],
          s_th2r[1],  s_th2r[2]))
    }
    
    ## BB8^180 parameter stats (no tail map here)
    if (p_bb8r > 0) {
      mask   <- fam_mat[, e] == BB8s
      th1_e  <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      w_b8   <- w_new[mask] / p_bb8r
      
      stats <- function(x) {
        mu <- sum(w_b8 * x)
        sd <- sqrt(sum(w_b8 * (x - mu)^2))
        ci <- w_quantile(x, w_b8, q_probs)
        c(mu, sd, ci)
      }
      s_th1 <- stats(th1_e)
      s_th2 <- stats(th2_e)
      
      res <- c(res,
               list(mu_bb8r_theta = s_th1[1], sd_bb8r_theta = s_th1[2],
                    bb8r_theta_q025 = s_th1[3], bb8r_theta_q975 = s_th1[4],
                    mu_bb8r_delta = s_th2[1], sd_bb8r_delta = s_th2[2],
                    bb8r_delta_q025 = s_th2[3], bb8r_delta_q975 = s_th2[4]))
      
      if (print_flag)
        cat(sprintf(
          "                BB8^180 : P=%.3f | θ=%.2f±%.2f CI=[%.2f,%.2f] | δ=%.2f±%.2f CI=[%.2f,%.2f]\n",
          p_bb8r,
          s_th1[1], s_th1[2], s_th1[3], s_th1[4],
          s_th2[1], s_th2[2], s_th2[3], s_th2[4]))
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



## Predictive

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


smc_predictive_sample <- function(particles,          # list of M particles
                                  skel,               # fixed structure
                                  w,                  # numeric(M) normalised
                                  L       = 10000,    # # draws wanted
                                  cache   = FALSE,    # reuse vine objects?
                                  cl      = NULL) {   # existing cluster
  
  stopifnot(abs(sum(w) - 1) < 1e-8)
  M  <- length(particles)
  d  <- ncol(skel$structure)
  
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
      c("id", "particles", "skel", "sim_row"),
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


smc_predictive_sample_fast <- function(particles,  # list of M particles
                                       skel,       # fixed structure
                                       w,          # numeric(M), not necessarily normalized
                                       L = 10000,  # total draws
                                       cl = NULL)  # optional cluster
{
  ## --- normalize weights robustly ---
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0
  sw <- sum(w)
  if (sw <= 0) stop("All weights are zero/invalid.")
  w <- w / sw
  
  M <- length(particles)
  if (length(w) != M) stop("length(w) must equal number of particles.")
  if (L < 1L) stop("L must be >= 1.")
  
  ## --- counts per particle (sum == L) ---
  cnt  <- as.vector(rmultinom(1, size = L, prob = w))  # length M
  used <- which(cnt > 0L)
  if (length(used) == 0L) {
    d0 <- ncol(skel$structure)
    return(matrix(numeric(0), nrow = 0, ncol = d0))
  }
  
  ## --- build each used vine ON MASTER (once) ---
  vines_used <- lapply(used, function(i) fast_vine_from_particle(particles[[i]], skel))
  
  ## helper: reshape rvinecop() output to (ni x d_sim), with d_sim inferred from length(sim)/ni
  coerce_block <- function(sim, ni) {
    vec <- as.numeric(sim)                       # flatten deterministically
    len <- length(vec)
    if (len %% ni != 0L) {
      stop(sprintf("rvinecop length mismatch: got %d, not divisible by ni=%d", len, ni))
    }
    d_sim <- as.integer(len / ni)
    dim(vec) <- c(ni, d_sim)
    vec
  }
  
  ## --- SERIAL PATH ---
  if (is.null(cl)) {
    sims  <- vector("list", length(used))
    d_ref <- NULL
    for (k in seq_along(used)) {
      ni  <- cnt[used[k]]
      sim <- rvinecop(ni, vines_used[[k]])
      blk <- coerce_block(sim, ni)
      d_k <- ncol(blk)
      if (is.null(d_ref)) d_ref <- d_k else if (d_k != d_ref)
        stop(sprintf("Inconsistent simulated dimension across particles: %d vs %d", d_k, d_ref))
      sims[[k]] <- blk
    }
    out <- do.call(rbind, sims)
    dimnames(out) <- NULL
    return(out)
  }
  
  ## --- PARALLEL PATH (load-balanced) ---
  # Pack only vine + ni (keep payload minimal)
  work <- lapply(seq_along(used), function(k) {
    list(vine = vines_used[[k]], ni = cnt[used[k]])
  })
  
  # Make sure workers can call rvinecop()
  parallel::clusterEvalQ(cl, library(rvinecopulib))
  
  sims_list <- parallel::parLapplyLB(cl, work, function(task) {
    sim <- rvinecop(task$ni, task$vine)
    vec <- as.numeric(sim)
    len <- length(vec)
    if (len %% task$ni != 0L)
      stop(sprintf("rvinecop length mismatch on worker: len=%d, ni=%d", len, task$ni))
    d_sim <- as.integer(len / task$ni)
    dim(vec) <- c(task$ni, d_sim)
    vec
  })
  
  ## check all blocks have same d, then rbind
  d_vec <- vapply(sims_list, ncol, integer(1))
  if (length(unique(d_vec)) != 1L)
    stop(sprintf("Inconsistent simulated dimension across workers: %s",
                 paste(unique(d_vec), collapse = ", ")))
  out <- do.call(rbind, sims_list)
  dimnames(out) <- NULL
  out
}









# Risk metrics

# update_risk_log <- function(risk_log,
#                             U_pred,
#                             mu_f,
#                             sig_f,
#                             t_idx,
#                             alphas = c(0.95, 0.975))
# {
#   # ------------------------------------------------------------------
#   #
#   #  Arguments
#   #  ----------
#   #  risk_log   data.table | initially empty (0 rows) and kept by caller
#   #  U_pred     L × d matrix of uniform draws from smc_predictive_sample()
#   #  mu_fc      numeric-vector length d   –– conditional means  (AR part)
#   #  sig_fc     numeric-vector length d   –– conditional stdevs (GARCH part)
#   #  t_idx      scalar, time index (stored in the log)
#   #  alphas     tail probabilities for VaR / ES   (default 95 % & 97.5 %)
#   #
#   #  Returns
#   #  -------
#   #  Same data.table, one extra row per call
#   # ------------------------------------------------------------------
#   stopifnot(is.matrix(U_pred),
#             length(mu_f)  == ncol(U_pred),
#             length(sig_f) == ncol(U_pred))
#   
#   
#   sig_f = as.numeric(sig_f)
#   mu_f = as.numeric(mu_f)
#   
#   # quantile transform & rescale 
#   Z      <- qnorm(U_pred)                           # L × d
#   R_pred <- sweep(Z,  2, sig_f, `*`)
#   R_pred <- sweep(R_pred, 2, mu_f,  `+`)
#   
#   L <- nrow(R_pred);  d <- ncol(R_pred)
#   series_names <- paste0("S", seq_len(d))           # generic column labels
#   dimnames(R_pred) <- list(NULL, series_names)
#   
#   # point moments 
#   mu_hat  <- colMeans(R_pred)
#   var_hat <- colMeans((R_pred - rep(mu_hat, each = L))^2)
#   
#   # tail risk 
#   losses <- -R_pred                          # L × d matrix
#   
#   VaR <- vapply(alphas,                          # (K × d)
#                 function(a)
#                   apply(losses, 2, quantile, probs = a),
#                 numeric(d))
#   
#   ES  <- vapply(seq_along(alphas),   # k = 1 … K
#                 function(k) {
#                   thr <- VaR[ , k]   # one threshold per series
#                   vapply(seq_len(d), function(j) {
#                     lj <- losses[ , j]
#                     exc <- lj[lj >= thr[j]]
#                     if (length(exc)) mean(exc) else NA_real_
#                   }, numeric(1))
#                 },
#                 numeric(d))
#   dimnames(VaR) <- dimnames(ES) <- list(colnames(losses), paste0("a", alphas))
#   
#   # 95 % predictive interval 
#   CI_lower <- apply(R_pred, 2, quantile, probs = 0.025)
#   CI_upper <- apply(R_pred, 2, quantile, probs = 0.975)
#   
#   # gather as one long row (wide format) 
#   row <- data.table::data.table(t_idx = t_idx)
#   
#   # names of the d series, e.g.  c("S1","S2", …)
#   series_names <- colnames(R_pred)
#   
#   # row is the single-row data.table that accumulates the statistics
#   add_vec <- function(x, prefix) {
#     cols <- paste0(prefix, "_", series_names)   # e.g.  mean_S1, mean_S2 …
#     stopifnot(length(x) == length(cols))        # sanity check
#     row[, (cols) := as.list(unname(x))]         # <- list = one entry / col
#   }
#   
#   # point estimates 
#   add_vec(mu_hat,  "mean")
#   add_vec(var_hat, "var")
#   add_vec(CI_lower,"ci_lo")
#   add_vec(CI_upper,"ci_hi")
#   
#   # tail–risk measures 
#   for (k in seq_along(alphas)) {
#     add_vec(VaR[ , k], sprintf("VaR%g", alphas[k] * 100))
#     add_vec(ES [ , k], sprintf("ES%g",  alphas[k] * 100))
#   }
#   
#   # append & return 
#   data.table::rbindlist(list(risk_log, row), use.names = TRUE, fill = TRUE)
# }

## Tree by tree functions

mh_step_in_tree <- function(p, tr, data_up_to_t, temp_skel, cfg) {
  
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
  
  prop <- p
  prop$theta <- p$theta + rnorm(cfg$K, 0, cfg$step_sd)
  
  vine_prop <- fast_vine_from_particle(prop, temp_skel)
  prop_ll   <- sum(log(pmax(dvinecop(data_up_to_t, vine_prop),
                            .Machine$double.eps)))
  
  vine_curr <- fast_vine_from_particle(p, temp_skel)
  curr_ll   <- sum(log(pmax(dvinecop(data_up_to_t, vine_curr),
                            .Machine$double.eps)))
  
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

calculate_log_lik_tree_tr <- function(particle, skel_tr, u_row, t, tr, tr_prev, skeletons_by_tr, cfg) {
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
  
  return(log(dens_val_tr + 1e-100) - log(dens_val_prev + 1e-100))
}
















