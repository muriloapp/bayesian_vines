library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(matrixStats)


FAM_INFO <- data.frame(
  name      = c("indep","gaussian","bb1","bb1r180","bb8r180","bb7","bb7r180"),
  rv_name   = c("indep","gaussian","bb1","bb1","bb8","bb7","bb7"),
  rotation  = c(0L, 0L, 0L, 180L, 180L, 0L, 180L),
  code      = c(0L, 1L, 3L, 13L, 16L, 7L, 17L),  # any unique ints are fine internally
  npar      = c(0L, 1L, 2L, 2L, 2L, 2L, 2L),
  stringsAsFactors = FALSE
)

## Shortcuts (optional)
FAM_INDEP    <- FAM_INFO$code[FAM_INFO$name == "indep"]
FAM_GAUSS    <- FAM_INFO$code[FAM_INFO$name == "gaussian"]
FAM_BB1      <- FAM_INFO$code[FAM_INFO$name == "bb1"]
FAM_BB1R180  <- FAM_INFO$code[FAM_INFO$name == "bb1r180"]
FAM_BB8R180  <- FAM_INFO$code[FAM_INFO$name == "bb8r180"]
FAM_BB7      <- FAM_INFO$code[FAM_INFO$name == "bb7"]          # NEW
FAM_BB7R180  <- FAM_INFO$code[FAM_INFO$name == "bb7r180"]      # NEW

## Minimal templates (we’ll build bicops from rv_name/rotation on the fly too)
T_INDEP     <- bicop_dist("indep")
T_GAUSS     <- bicop_dist("gaussian", parameters = 0)

T_BB1       <- bicop_dist("bb1", rotation = 0,   parameters = c(1, 2))
T_BB1R180   <- bicop_dist("bb1", rotation = 180, parameters = c(1, 2))

T_BB8R180   <- bicop_dist("bb8", rotation = 180, parameters = c(1, 1))  # dummy two-par

T_BB7       <- bicop_dist("bb7", rotation = 0,   parameters = c(1, 2))   # NEW
T_BB7R180   <- bicop_dist("bb7", rotation = 180, parameters = c(1, 2))   # NEW


log_prior_edge <- function(fam, th1, th2, cfg) {
  if (fam == FAM_INDEP) return(0)
  
  if (fam == FAM_GAUSS) {
    return(log(1/1.98)) #- cfg$lambda)
  }
  
  if (fam == FAM_BB1) {
    lam <- bb1_par2tail(th1, th2)   # (λL, λU)
    lp_tail <- dbeta(lam[1], 2, 2, log=TRUE) + dbeta(lam[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb1_log_jacobian(lam[1], lam[2]) #- 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  if (fam == FAM_BB1R180) {
    lamr <- bb1r180_par2tail(th1, th2)   # (λL^rot, λU^rot)
    lp_tail <- dbeta(lamr[1], 2, 2, log=TRUE) + dbeta(lamr[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb1r180_log_jacobian(lamr[1], lamr[2]) #- 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  if (fam == FAM_BB7) {
    lam <- bb7_par2tail(th1, th2)   # (λL, λU)
    lp_tail <- dbeta(lam[1], 2, 2, log=TRUE) + dbeta(lam[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb7_log_jacobian(lam[1], lam[2]) #- 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  if (fam == FAM_BB7R180) {
    lamr <- bb7r180_par2tail(th1, th2) # (λL^rot, λU^rot) = (λU, λL)
    lp_tail <- dbeta(lamr[1], 2, 2, log=TRUE) + dbeta(lamr[2], 2, 2, log=TRUE)
    lp <- lp_tail + bb7r180_log_jacobian(lamr[1], lamr[2]) #- 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  if (fam == FAM_BB8R180) {
    lamLr <- bb8r180_par2tail(th1, th2)[1]  # only λ_L^rot is nonzero
    lp_tail <- dbeta(lamLr, 2, 2, log=TRUE)
    ## add 1D Jacobian for θ(λ_L^rot) + complexity penalty for 2 params
    lp <- lp_tail + bb8r180_log_jacobian_1d(lamLr) #- 2*cfg$lambda
    return(fillna_neg(lp))
  }
  
  stop("Unknown family code in log_prior_edge: ", fam)
}



log_prior <- function(p, cfg) {
  sum(mapply(log_prior_edge,
             p$fam, p$th1, p$th2,
             MoreArgs = list(cfg = cfg)))
}


# new_particle <- function(cfg) {
#   K   <- cfg$K
#   fam <- integer(K); th1 <- numeric(K); th2 <- numeric(K)
#   
#   for (k in seq_len(K)) {
#     tr_k <- cfg$edge_tree[k]
#     allowed_names <- if (tr_k == 1) cfg$families_first else cfg$families_deep
#     fam_tbl <- FAM_INFO[FAM_INFO$name %in% allowed_names, ]
#     #weights <- exp(-cfg$lambda * fam_tbl$npar)
#     weights <- rep(1/dim(fam_tbl)[1], dim(fam_tbl)[1])
#     code_k  <- sample(fam_tbl$code, 1, prob = weights)
#     fam[k]  <- code_k
#     
#     if (code_k == FAM_INDEP) {
#       th1[k] <- 0; th2[k] <- 0
#       
#     } else if (code_k == FAM_GAUSS) {
#       th1[k] <- runif(1, -0.99, 0.99); th2[k] <- 0
#       
#     } else if (code_k == FAM_BB1) {
#       lamU <- rbeta(1, 2, 2); lamL <- rbeta(1, 2, 2)
#       pars <- sanitize_bb1(bb1_tail2par(lamL, lamU)[1],
#                            bb1_tail2par(lamL, lamU)[2])
#       th1[k] <- pars[1]; th2[k] <- pars[2]
#       
#     } else if (code_k == FAM_BB1R180) {
#       # sample in rotated tail space, then map with rotated converter
#       lamL_rot <- rbeta(1, 2, 2); lamU_rot <- rbeta(1, 2, 2)
#       pars <- sanitize_bb1(bb1r180_tail2par(lamL_rot, lamU_rot)[1],
#                            bb1r180_tail2par(lamL_rot, lamU_rot)[2])
#       th1[k] <- pars[1]; th2[k] <- pars[2]
#       
#     } else if (code_k == FAM_BB8R180) {
#       ## λ_L^rot ~ Beta, map to θ; δ (Frank) is tail-free → draw diffuse
#       lamLr <- rbeta(1, 2, 2)
#       pars  <- bb8r180_tail2par(lamLr, delta_free = runif(1, 1.1, 5.0))
#       pars  <- sanitize_bb8(pars[1], pars[2])
#       th1[k] <- pars[1]; th2[k] <- pars[2]
#       
#     } else if (code_k == FAM_BB7) {
#       ## (λL, λU) ~ Beta×Beta, map to (θ, δ)
#       lamL <- rbeta(1, 2, 2); lamU <- rbeta(1, 2, 2)
#       pars <- bb7_tail2par(lamL, lamU)
#       pars <- sanitize_bb7(pars[1], pars[2])
#       th1[k] <- pars[1]; th2[k] <- pars[2]
#       
#     } else if (code_k == FAM_BB7R180) {
#       ## (λL^rot, λU^rot) ~ Beta×Beta, map to (θ, δ) via rotated transform
#       lamLr <- rbeta(1, 2, 2); lamUr <- rbeta(1, 2, 2)
#       pars  <- bb7r180_tail2par(lamLr, lamUr)
#       pars  <- sanitize_bb7(pars[1], pars[2])
#       th1[k] <- pars[1]; th2[k] <- pars[2]
#     }
#   }
#   
#   list(fam = fam, th1 = th1, th2 = th2, w = 1 / cfg$M, last_accept = FALSE)
# }

new_particle <- function(cfg, U_init = NULL) {
  K   <- cfg$K
  fam <- integer(K); th1 <- numeric(K); th2 <- numeric(K)
  
  tip_on  <- isTRUE(cfg$use_tail_informed_prior) && !is.null(U_init)
  if (tip_on) {
    if (is.null(cfg$edge_pair))
      stop("TIP init requested but cfg$edge_pair is NULL. Provide Tree-1 pairs aligned with edge order.")
  }
  
  for (k in seq_len(K)) {
    tr_k <- cfg$edge_tree[k]
    allowed_names <- if (tr_k == 1L) cfg$families_first else cfg$families_deep
    fam_tbl <- FAM_INFO[FAM_INFO$name %in% allowed_names, , drop = FALSE]
    
    ## --------- FAMILY SELECTION (TIP-aware on Tree 1) -----------------
    if (tip_on && tr_k == 1L) {
      # local index within Tree 1 and its (i,j) pair
      local_idx <- sum(cfg$edge_tree[seq_len(k)] == 1L)
      pair      <- cfg$edge_pair[local_idx, ]
      uv        <- U_init[, pair, drop = FALSE]
      
      # empirical tails from FRAPO (uses only U_init)
      emps <- emp_tails_FRAPO(uv, method = cfg$tip_method, k = cfg$tip_k)
      mL   <- as.numeric(emps["L"]);  mU <- as.numeric(emps["U"])
      tail_strength <- max(mL, mU)             # ∈ (0,1)
      
      # weight tail families up with tail_strength; keep others possible
      tail_codes <- c(FAM_BB1, FAM_BB1R180, FAM_BB7, FAM_BB7R180, FAM_BB8R180)
      base_w     <- rep(1.0, nrow(fam_tbl))
      is_tailfam <- fam_tbl$code %in% tail_codes
      # e.g., tails strong -> upweight tail fams by up to 4x; downweight Gauss/Indep
      w          <- base_w * ifelse(is_tailfam, 1 + 3*tail_strength, pmax(1 - 0.9*tail_strength, 0.1))
      w          <- w / sum(w)
      
      code_k  <- sample(fam_tbl$code, 1L, prob = w)
      fam[k]  <- code_k
      
      ## --------- PARAMETER SEEDING from empirical tails ----------------
      if (code_k == FAM_INDEP) {
        th1[k] <- 0; th2[k] <- 0
        
      } else if (code_k == FAM_GAUSS) {
        th1[k] <- runif(1, -0.99, 0.99); th2[k] <- 0
        
      } else {
        # strong seeding for tail families (logit-normal around mL/mU)
        th_new <- seed_family_from_emp(
          new_code     = code_k,
          mL           = mL,
          mU           = mU,
          cur_delta    = runif(1, 1.1, 5.0),   # starting delta if needed (BB8r180)
          step_sd      = cfg$step_sd,
          tip_sd_logit = cfg$tip_sd_logit      # small sd -> strong prior around empirical tails
        )
        th1[k] <- th_new[1]; th2[k] <- th_new[2]
      }
      
    } else {
      ## --------- ORIGINAL (non-TIP or deeper trees) -------------------
      # Equal weights across allowed families (your current behavior)
      weights <- rep(1 / nrow(fam_tbl), nrow(fam_tbl))
      code_k  <- sample(fam_tbl$code, 1L, prob = weights)
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
        lamL_rot <- rbeta(1, 2, 2); lamU_rot <- rbeta(1, 2, 2)
        pars <- sanitize_bb1(bb1r180_tail2par(lamL_rot, lamU_rot)[1],
                             bb1r180_tail2par(lamL_rot, lamU_rot)[2])
        th1[k] <- pars[1]; th2[k] <- pars[2]
        
      } else if (code_k == FAM_BB8R180) {
        lamLr <- rbeta(1, 2, 2)
        pars  <- bb8r180_tail2par(lamLr, delta_free = runif(1, 1.1, 5.0))
        pars  <- sanitize_bb8(pars[1], pars[2])
        th1[k] <- pars[1]; th2[k] <- pars[2]
        
      } else if (code_k == FAM_BB7) {
        lamL <- rbeta(1, 2, 2); lamU <- rbeta(1, 2, 2)
        pars <- bb7_tail2par(lamL, lamU)
        pars <- sanitize_bb7(pars[1], pars[2])
        th1[k] <- pars[1]; th2[k] <- pars[2]
        
      } else if (code_k == FAM_BB7R180) {
        lamLr <- rbeta(1, 2, 2); lamUr <- rbeta(1, 2, 2)
        pars  <- bb7r180_tail2par(lamLr, lamUr)
        pars  <- sanitize_bb7(pars[1], pars[2])
        th1[k] <- pars[1]; th2[k] <- pars[2]
      }
    }
  }
  
  list(fam = fam, th1 = th1, th2 = th2, w = 1 / cfg$M, last_accept = FALSE)
}




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
        
      } else if (code == FAM_BB7) {                       # NEW
        bc <- rlang::duplicate(T_BB7, shallow = TRUE)
        bc$parameters[1:2] <- c(p$th1[idx], p$th2[idx])
        pcs[[tr]][[ed]] <- bc
        
      } else if (code == FAM_BB7R180) {                   # NEW
        bc <- rlang::duplicate(T_BB7R180, shallow = TRUE)
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
    li           <- ifelse(is.na(dens_val), -1e100,
                                      ifelse(dens_val <= 0, -1e100, log(dens_val)))
    
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








## --- helpers used by the flip init ------------------------------------------------
clamp01 <- function(x, lo = 1e-3, hi = 0.99) pmin(pmax(x, lo), hi)

get_tails <- function(code, th1, th2) {
  if (code == FAM_BB1)      { v <- bb1_par2tail(th1, th2);      return(list(L = v[1], U = v[2])) }
  if (code == FAM_BB1R180)  { v <- bb1r180_par2tail(th1, th2);  return(list(Lr = v[1], Ur = v[2])) }
  if (code == FAM_BB7)      { v <- bb7_par2tail(th1, th2);      return(list(L = v[1], U = v[2])) }
  if (code == FAM_BB7R180)  { v <- bb7r180_par2tail(th1, th2);  return(list(Lr = v[1], Ur = v[2])) }
  if (code == FAM_BB8R180)  { v <- bb8r180_par2tail(th1, th2);  return(list(Lr = v[1])) }
  list()  # Indep/Gauss
}

init_from_tails <- function(new_code, tails) {
  draw_beta <- function() rbeta(1, 2, 2)
  draw_delta_bb8 <- function() runif(1, 1.1, 5.0)
  
  if (new_code == FAM_INDEP)  return(c(0, 0))
  if (new_code == FAM_GAUSS)  return(c(runif(1, -0.99, 0.99), 0))
  
  if (new_code == FAM_BB1) {
    lamL <- clamp01(if (!is.null(tails$L))  tails$L else draw_beta())
    lamU <- clamp01(if (!is.null(tails$U))  tails$U else draw_beta())
    tp   <- bb1_tail2par(lamL, lamU);  return(sanitize_bb1(tp[1], tp[2]))
  }
  if (new_code == FAM_BB1R180) {
    lamLr <- clamp01(if (!is.null(tails$Lr)) tails$Lr else if (!is.null(tails$L)) tails$L else draw_beta())
    lamUr <- clamp01(if (!is.null(tails$Ur)) tails$Ur else if (!is.null(tails$U)) tails$U else draw_beta())
    tp    <- bb1r180_tail2par(lamLr, lamUr); return(sanitize_bb1(tp[1], tp[2]))
  }
  if (new_code == FAM_BB7) {
    lamL <- clamp01(if (!is.null(tails$L))  tails$L else draw_beta())
    lamU <- clamp01(if (!is.null(tails$U))  tails$U else draw_beta())
    tp   <- bb7_tail2par(lamL, lamU);  return(sanitize_bb7(tp[1], tp[2]))
  }
  if (new_code == FAM_BB7R180) {
    lamLr <- clamp01(if (!is.null(tails$Lr)) tails$Lr else if (!is.null(tails$L)) tails$L else draw_beta())
    lamUr <- clamp01(if (!is.null(tails$Ur)) tails$Ur else if (!is.null(tails$U)) tails$U else draw_beta())
    tp    <- bb7r180_tail2par(lamLr, lamUr); return(sanitize_bb7(tp[1], tp[2]))
  }
  if (new_code == FAM_BB8R180) {
    lamLr <- clamp01(if (!is.null(tails$Lr)) tails$Lr else if (!is.null(tails$L)) tails$L else draw_beta())
    delta <- if (!is.null(tails$delta_hint)) tails$delta_hint else draw_delta_bb8()
    tp    <- bb8r180_tail2par(lamLr, delta); return(sanitize_bb8(tp[1], tp[2]))
  }
  stop("Unknown family code in init_from_tails")
}

mh_step <- function(p, data_up_to_t, skeleton, cfg) {
  
  K    <- cfg$K
  prop <- p
  
  ## (a) Family flips — tail-seeded init (Tree 1 only; no MLE)
  end_first_tr <- sum(cfg$edge_tree == 1L)  # number of edges in Tree 1
  if (end_first_tr > 0L) {
    flip_mask <- runif(end_first_tr) < cfg$q_flip
    if (any(flip_mask)) {
      idxs <- which(flip_mask)
      for (idx in idxs) {
        tr_idx <- cfg$edge_tree[idx]
        if (tr_idx != 1L) next  # safety; we only flip on Tree 1
        
        old_code <- prop$fam[idx]
        allowed_names <- if (tr_idx == 1L) cfg$families_first else cfg$families_deep
        allowed_codes <- FAM_INFO$code[FAM_INFO$name %in% allowed_names]
        allowed_codes <- setdiff(allowed_codes, old_code)  # avoid no-op
        if (!length(allowed_codes)) next
        
        new_code <- sample(allowed_codes, 1L)
        prop$fam[idx] <- new_code
        
        if (isTRUE(cfg$tip_init)) {
          ## find the local Tree 1 index and the pair for THIS edge
          local_idx <- sum(cfg$edge_tree[seq_len(idx)] == 1L)  # 1-based within Tree 1
          pair      <- cfg$edge_pair[local_idx, ]              # E1 x 2; aligned to your order
          uv        <- data_up_to_t[, pair, drop = FALSE]
          
          emps <- emp_tails_FRAPO(uv, method = cfg$tip_method, k = cfg$tip_k)
          emps
          mL   <- emps["L"];  mU <- emps["U"]
          
          th_new <- seed_family_from_emp(
            new_code, mL, mU,
            cur_delta   = prop$th2[idx],
            step_sd     = cfg$step_sd,
            tip_sd_logit= cfg$tip_sd_logit
          )
          prop$th1[idx] <- th_new[1]; prop$th2[idx] <- th_new[2]
        } else {
        
          ## seed the proposal using CURRENT tails when available; otherwise prior
          tails_old <- get_tails(old_code, prop$th1[idx], prop$th2[idx])
          tp_new    <- init_from_tails(new_code, tails_old)
          
          prop$th1[idx] <- tp_new[1]
          prop$th2[idx] <- tp_new[2]
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
  
  ### BB8^180 — RW in λ_L^rot only (the only nonzero tail), plus small RW on δ
  bb8r_idx <- which(prop$fam == FAM_BB8R180)
  if (length(bb8r_idx)) {
    lamLr <- vapply(bb8r_idx, function(i) bb8r180_par2tail(prop$th1[i], prop$th2[i])[1], numeric(1))
    lamLr <- rtnorm_vec(lamLr, cfg$step_sd, 1e-3, 0.99)
    pars  <- t(vapply(seq_along(bb8r_idx), function(j) {
      i <- bb8r_idx[j]
      ## update θ from λ_L^rot; δ small RW directly
      delta_prop <- pmin(pmax(prop$th2[i] + rnorm(1, 0, cfg$step_sd), 0.05), 7.0)
      bb8r180_tail2par(lamLr[j], delta_prop)
    }, numeric(2)))
    pars <- t(apply(pars, 1, \(x) sanitize_bb8(x[1], x[2])))
    prop$th1[bb8r_idx] <- pars[,1]; prop$th2[bb8r_idx] <- pars[,2]
  }
  
  ## BB7 — RW in (λL, λU), then map back and sanitize
  bb7_idx <- which(prop$fam == FAM_BB7)
  if (length(bb7_idx)) {
    lam <- t(mapply(bb7_par2tail, prop$th1[bb7_idx], prop$th2[bb7_idx]))  # [n×2]
    lam[,1] <- rtnorm_vec(lam[,1], cfg$step_sd, 1e-3, 0.99)  # λL
    lam[,2] <- rtnorm_vec(lam[,2], cfg$step_sd, 1e-3, 0.99)  # λU
    par <- t(apply(lam, 1, \(x) bb7_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb7(x[1], x[2])))
    prop$th1[bb7_idx] <- par[,1]; prop$th2[bb7_idx] <- par[,2]
  }
  
  ## BB7^180 — RW in (λL^rot, λU^rot), swap map back, sanitize
  bb7r_idx <- which(prop$fam == FAM_BB7R180)
  if (length(bb7r_idx)) {
    lamr <- t(mapply(bb7r180_par2tail, prop$th1[bb7r_idx], prop$th2[bb7r_idx]))
    lamr[,1] <- rtnorm_vec(lamr[,1], cfg$step_sd, 1e-3, 0.99)  # λL^rot
    lamr[,2] <- rtnorm_vec(lamr[,2], cfg$step_sd, 1e-3, 0.99)  # λU^rot
    par <- t(apply(lamr, 1, \(x) bb7r180_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb7(x[1], x[2])))
    prop$th1[bb7r_idx] <- par[,1]; prop$th2[bb7r_idx] <- par[,2]
  }
  
  
  ## (c) MH ratio
  vine_prop <- fast_vine_from_particle(prop, skeleton)
  vine_curr <- fast_vine_from_particle(p,    skeleton)
  
  #ll_prop <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_prop), .Machine$double.eps))))
  #ll_curr <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_curr), .Machine$double.eps))))

  if (isTRUE(cfg$use_weighted_ll)) {
    wobs    <- tail_weights(
      data_up_to_t,
      tau = cfg$tauL,
      k   = cfg$joint_k,
      eps = cfg$tail_eps
    )
    ll_prop <- sum(wobs * safe_logdens(data_up_to_t, vine_prop))
    ll_curr <- sum(wobs * safe_logdens(data_up_to_t, vine_curr))
  } else {
    ll_prop <- sum(safe_logdens(data_up_to_t, vine_prop))
    ll_curr <- sum(safe_logdens(data_up_to_t, vine_curr))
  }
  
  
  # --- PRIOR (toggle-safe) ----------------------------------------------------
  lp_prop <- if (isTRUE(cfg$use_tail_informed_prior)) 
    log_prior_with_tip_time(prop, cfg, data_up_to_t) else 
      log_prior(prop, cfg)
  
  lp_curr <- if (isTRUE(cfg$use_tail_informed_prior)) 
    log_prior_with_tip_time(p,    cfg, data_up_to_t) else 
      log_prior(p,    cfg)
  
  log_acc <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
  
  if (log(runif(1)) < log_acc) {
    prop$last_accept <- TRUE
    return(prop)
  } else {
    p$last_accept <- FALSE
    return(p)
  }
}


## ---- Tail weighting helpers ----------------------------------------------

## Soft weights: 1 for "joint lower-tail" rows, eps for others (renormalized)
tail_weights <- function(U, tau = 0.05, k = 2L, eps = 0.10) {
  jt <- rowSums(U <= tau) >= k     # joint lower-tail indicator per row
  w  <- ifelse(jt, 1, eps)
  w / mean(w)                      # keep average weight = 1 (stability)
}

## Safe log-density vector for a vine on matrix U
## Follows your convention: NA or <= 0  → -1e100
safe_logdens <- function(U, vine) {
  dens <- dvinecop(U, vine)
  ifelse(is.na(dens) | dens <= 0, -1e100, log(dens))
}


###########################################

## --- math helpers ----------------------------------------------------------
logit   <- function(x)  log(x/(1-x))
ilogit  <- function(z)  1/(1+exp(-z))
clamp01 <- function(x, lo = 1e-6, hi = 1-1e-6) pmin(pmax(x, lo), hi)

## Logit-Normal log-density on (0,1)
dlogitnorm <- function(lambda, mean_logit, sd_logit) {
  lam <- clamp01(lambda)
  dnorm(logit(lam), mean_logit, sd_logit, log = TRUE) - log(lam) - log1p(-lam)
}

## --- FRAPO-based empirical tails for a 2-column matrix uv ------------------
emp_tails_FRAPO <- function(uv, method = "EmpTC", k = NULL) {
  stopifnot(is.matrix(uv), ncol(uv) == 2)
  # lower tail
  tL <- FRAPO::tdc(uv, method = method, lower = TRUE,  k = k)
  # upper tail
  tU <- FRAPO::tdc(uv, method = method, lower = FALSE, k = k)
  # robust extraction (FRAPO returns a 2x2 matrix for 2 columns)
  get12 <- function(x) if (is.matrix(x)) as.numeric(x[1,2]) else as.numeric(x)[1]
  c(L = clamp01(get12(tL)), U = clamp01(get12(tU)))
}

## --- Seeding new family params from empirical tails (Tree 1) ---------------
seed_family_from_emp <- function(new_code, mL, mU, cur_delta, step_sd, tip_sd_logit) {
  r_lam <- function(m) {
    z <- rnorm(1, mean = logit(clamp01(m)), sd = tip_sd_logit)
    ilogit(z)
  }
  if (new_code == FAM_INDEP)  return(c(0, 0))
  if (new_code == FAM_GAUSS)  return(c(runif(1, -0.99, 0.99), 0))
  
  if (new_code == FAM_BB1) {
    lamL <- r_lam(mL); lamU <- r_lam(mU)
    th   <- bb1_tail2par(lamL, lamU); th <- sanitize_bb1(th[1], th[2]); return(th)
  }
  if (new_code == FAM_BB7) {
    lamL <- r_lam(mL); lamU <- r_lam(mU)
    th   <- bb7_tail2par(lamL, lamU); th <- sanitize_bb7(th[1], th[2]); return(th)
  }
  if (new_code == FAM_BB1R180) {
    # rotation: lower-rot ↔ upper(orig), upper-rot ↔ lower(orig)
    lamLr <- r_lam(mU); lamUr <- r_lam(mL)
    th    <- bb1r180_tail2par(lamLr, lamUr); th <- sanitize_bb1(th[1], th[2]); return(th)
  }
  if (new_code == FAM_BB7R180) {
    lamLr <- r_lam(mU); lamUr <- r_lam(mL)
    th    <- bb7r180_tail2par(lamLr, lamUr); th <- sanitize_bb7(th[1], th[2]); return(th)
  }
  if (new_code == FAM_BB8R180) {
    lamLr <- r_lam(mL)  # choose L as anchor
    delta <- pmin(pmax(cur_delta + rnorm(1, 0, step_sd), 0.05), 7.0)
    th    <- bb8r180_tail2par(lamLr, delta); th <- sanitize_bb8(th[1], th[2]); return(th)
  }
  stop("Unknown family code in seed_family_from_emp()")
}



log_prior_edge_strong <- function(fam, th1, th2, cfg, tip_means) {
  if (fam == FAM_INDEP) return(0 - cfg$lambda)
  if (fam == FAM_GAUSS) return(log(1/1.98) - cfg$lambda)  # uniform ρ in (-.99,.99) + penalty
  
  s <- cfg$tip_sd_logit
  if (fam == FAM_BB1) {
    lam <- bb1_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mL"]), s) +
      dlogitnorm(lam[2], logit(tip_means["mU"]), s)
    return(tip + bb1_log_jacobian(lam[1], lam[2]) - cfg$lambda)
  }
  if (fam == FAM_BB7) {
    lam <- bb7_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mL"]), s) +
      dlogitnorm(lam[2], logit(tip_means["mU"]), s)
    return(tip + bb7_log_jacobian(lam[1], lam[2]) - cfg$lambda)
  }
  if (fam == FAM_BB1R180) {
    lam <- bb1r180_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mU"]), s) +  # rotated mapping
      dlogitnorm(lam[2], logit(tip_means["mL"]), s)
    return(tip + bb1_log_jacobian(lam[1], lam[2]) - cfg$lambda)
  }
  if (fam == FAM_BB7R180) {
    lam <- bb7r180_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mU"]), s) +
      dlogitnorm(lam[2], logit(tip_means["mL"]), s)
    return(tip + bb7_log_jacobian(lam[1], lam[2]) - cfg$lambda)
  }
  if (fam == FAM_BB8R180) {
    lamLr <- bb8r180_par2tail(th1, th2)[1]
    tip   <- if (is.null(tip_means)) 0 else dlogitnorm(lamLr, logit(tip_means["mL"]), s)
    return(tip + bb8_log_jacobian(lamLr) - cfg$lambda)
  }
  stop("unknown family in log_prior_edge_strong")
}

## empirical means for edge e at time t (Tree 1 only; NULL otherwise)
.tip_means_for_edge_t <- function(e, data_up_to_t, cfg) {
  if (!isTRUE(cfg$use_tail_informed_prior)) return(NULL)
  if (cfg$edge_tree[e] != 1L) return(NULL)
  ## map global edge index e -> local index within Tree 1, then pick its (i,j) pair
  local_idx <- sum(cfg$edge_tree[seq_len(e)] == 1L)
  pair      <- cfg$edge_pair[local_idx, ]
  uv        <- data_up_to_t[, pair, drop = FALSE]
  emps      <- emp_tails_FRAPO(uv, method = cfg$tip_method, k = cfg$tip_k)
  c(mL = emps["L"], mU = emps["U"])
}

## TIP-aware prior (falls back to your base prior when TIP is OFF)
log_prior_with_tip_time <- function(p, cfg, data_up_to_t) {
  if (!isTRUE(cfg$use_tail_informed_prior)) return(log_prior(p, cfg))
  lp <- 0
  for (e in seq_len(cfg$K)) {
    fam  <- p$fam[e]; th1 <- p$th1[e]; th2 <- p$th2[e]
    tipm <- .tip_means_for_edge_t(e, data_up_to_t, cfg)  # NULL if not Tree 1
    lp   <- lp + log_prior_edge_strong(fam, th1, th2, cfg, tipm)
  }
  lp
}



############################################3












# mh_step <- function(p, data_up_to_t, skeleton, cfg) {
#   K    <- cfg$K
#   prop <- p
#   
#   ## (a) Family flips — prior-only init (no MLE, no data-length gating)
#   end_first_tr <- sum(cfg$edge_tree == 1L)  # #edges in Tree 1 (by your ordering)
#   if (end_first_tr > 0L) {
#     flip_mask <- runif(end_first_tr) < cfg$q_flip
#     if (any(flip_mask)) {
#       idxs <- which(flip_mask)
#       for (idx in idxs) {
#         tr_idx <- cfg$edge_tree[idx]
#         if (tr_idx != 1L) next  # keep flips on Tree 1 only
#         
#         old_code <- prop$fam[idx]
#         allowed_names <- if (tr_idx == 1L) cfg$families_first else cfg$families_deep
#         allowed_codes <- FAM_INFO$code[FAM_INFO$name %in% allowed_names]
#         allowed_codes <- setdiff(allowed_codes, old_code)  # avoid no-op flip
#         if (!length(allowed_codes)) next
#         
#         new_code <- sample(allowed_codes, 1L)
#         prop$fam[idx] <- new_code
#         
#         ## Initialize from priors (no MLE)
#         if (new_code == FAM_INDEP) {
#           prop$th1[idx] <- 0;                       prop$th2[idx] <- 0
#           
#         } else if (new_code == FAM_GAUSS) {
#           prop$th1[idx] <- runif(1, -0.99, 0.99);   prop$th2[idx] <- 0
#           
#         } else if (new_code == FAM_BB1) {
#           lamL <- rbeta(1, 2, 2); lamU <- rbeta(1, 2, 2)
#           tp   <- bb1_tail2par(lamL, lamU)
#           tp   <- sanitize_bb1(tp[1], tp[2])
#           prop$th1[idx] <- tp[1];                  prop$th2[idx] <- tp[2]
#           
#         } else if (new_code == FAM_BB1R180) {
#           lamLr <- rbeta(1, 2, 2); lamUr <- rbeta(1, 2, 2)
#           tp    <- bb1r180_tail2par(lamLr, lamUr)
#           tp    <- sanitize_bb1(tp[1], tp[2])
#           prop$th1[idx] <- tp[1];                  prop$th2[idx] <- tp[2]
#           
#         } else if (new_code == FAM_BB8R180) {
#           lamLr <- rbeta(1, 2, 2)
#           tp    <- bb8r180_tail2par(lamLr, delta_free = runif(1, 1.1, 5.0))
#           tp    <- sanitize_bb8(tp[1], tp[2])
#           prop$th1[idx] <- tp[1]; prop$th2[idx] <- tp[2]
#           
#         } else if (new_code == FAM_BB7) {
#           lamL <- rbeta(1, 2, 2); lamU <- rbeta(1, 2, 2)
#           tp   <- bb7_tail2par(lamL, lamU)
#           tp   <- sanitize_bb7(tp[1], tp[2])
#           prop$th1[idx] <- tp[1];  prop$th2[idx] <- tp[2]
#           
#         } else if (new_code == FAM_BB7R180) {
#           lamLr <- rbeta(1, 2, 2); lamUr <- rbeta(1, 2, 2)
#           tp    <- bb7r180_tail2par(lamLr, lamUr)
#           tp    <- sanitize_bb7(tp[1], tp[2])
#           prop$th1[idx] <- tp[1];  prop$th2[idx] <- tp[2]
#         }
#       }
#     }
#   }
#   
#   ## (b) Random-walk on parameters of current families
#   # Gaussian — truncated normal on ρ
#   ga_idx <- which(prop$fam == FAM_GAUSS)
#   if (length(ga_idx))
#     prop$th1[ga_idx] <- rtnorm_vec(prop$th1[ga_idx], cfg$step_sd, -0.99, 0.99)
#   
#   # BB1 (0°) — RW in (λL, λU), then map back and sanitize
#   bb1_idx <- which(prop$fam == FAM_BB1)
#   if (length(bb1_idx)) {
#     lam <- t(mapply(bb1_par2tail, prop$th1[bb1_idx], prop$th2[bb1_idx]))
#     lam[, 1] <- rtnorm_vec(lam[, 1], cfg$step_sd, 1e-3, 0.99)
#     lam[, 2] <- rtnorm_vec(lam[, 2], cfg$step_sd, 1e-3, 0.99)
#     par <- t(apply(lam, 1, \(x) bb1_tail2par(x[1], x[2])))
#     par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
#     prop$th1[bb1_idx] <- par[, 1]; prop$th2[bb1_idx] <- par[, 2]
#   }
#   
#   # BB1^180 — RW in (λL_rot, λU_rot), map back with rotated converter
#   bb1r_idx <- which(prop$fam == FAM_BB1R180)
#   if (length(bb1r_idx)) {
#     lamr <- t(mapply(bb1r180_par2tail, prop$th1[bb1r_idx], prop$th2[bb1r_idx]))
#     lamr[, 1] <- rtnorm_vec(lamr[, 1], cfg$step_sd, 1e-3, 0.99)
#     lamr[, 2] <- rtnorm_vec(lamr[, 2], cfg$step_sd, 1e-3, 0.99)
#     par <- t(apply(lamr, 1, \(x) bb1r180_tail2par(x[1], x[2])))
#     par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
#     prop$th1[bb1r_idx] <- par[, 1]; prop$th2[bb1r_idx] <- par[, 2]
#   }
#   
#   ### BB8^180 — RW in λ_L^rot only (the only nonzero tail), plus small RW on δ
#   bb8r_idx <- which(prop$fam == FAM_BB8R180)
#   if (length(bb8r_idx)) {
#     lamLr <- vapply(bb8r_idx, function(i) bb8r180_par2tail(prop$th1[i], prop$th2[i])[1], numeric(1))
#     lamLr <- rtnorm_vec(lamLr, cfg$step_sd, 1e-3, 0.99)
#     pars  <- t(vapply(seq_along(bb8r_idx), function(j) {
#       i <- bb8r_idx[j]
#       ## update θ from λ_L^rot; δ small RW directly
#       delta_prop <- pmin(pmax(prop$th2[i] + rnorm(1, 0, cfg$step_sd), 0.05), 7.0)
#       bb8r180_tail2par(lamLr[j], delta_prop)
#     }, numeric(2)))
#     pars <- t(apply(pars, 1, \(x) sanitize_bb8(x[1], x[2])))
#     prop$th1[bb8r_idx] <- pars[,1]; prop$th2[bb8r_idx] <- pars[,2]
#   }
#   
#   ## BB7 — RW in (λL, λU), then map back and sanitize
#   bb7_idx <- which(prop$fam == FAM_BB7)
#   if (length(bb7_idx)) {
#     lam <- t(mapply(bb7_par2tail, prop$th1[bb7_idx], prop$th2[bb7_idx]))  # [n×2]
#     lam[,1] <- rtnorm_vec(lam[,1], cfg$step_sd, 1e-3, 0.99)  # λL
#     lam[,2] <- rtnorm_vec(lam[,2], cfg$step_sd, 1e-3, 0.99)  # λU
#     par <- t(apply(lam, 1, \(x) bb7_tail2par(x[1], x[2])))
#     par <- t(apply(par, 1, \(x) sanitize_bb7(x[1], x[2])))
#     prop$th1[bb7_idx] <- par[,1]; prop$th2[bb7_idx] <- par[,2]
#   }
#   
#   ## BB7^180 — RW in (λL^rot, λU^rot), swap map back, sanitize
#   bb7r_idx <- which(prop$fam == FAM_BB7R180)
#   if (length(bb7r_idx)) {
#     lamr <- t(mapply(bb7r180_par2tail, prop$th1[bb7r_idx], prop$th2[bb7r_idx]))
#     lamr[,1] <- rtnorm_vec(lamr[,1], cfg$step_sd, 1e-3, 0.99)  # λL^rot
#     lamr[,2] <- rtnorm_vec(lamr[,2], cfg$step_sd, 1e-3, 0.99)  # λU^rot
#     par <- t(apply(lamr, 1, \(x) bb7r180_tail2par(x[1], x[2])))
#     par <- t(apply(par, 1, \(x) sanitize_bb7(x[1], x[2])))
#     prop$th1[bb7r_idx] <- par[,1]; prop$th2[bb7r_idx] <- par[,2]
#   }
#   
#   
#   ## (c) MH ratio
#   vine_prop <- fast_vine_from_particle(prop, skeleton)
#   vine_curr <- fast_vine_from_particle(p,    skeleton)
#   
#   ll_prop <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_prop), .Machine$double.eps))))
#   ll_curr <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_curr), .Machine$double.eps))))
#   
#   log_acc <- (ll_prop + log_prior(prop, cfg)) - (ll_curr + log_prior(p, cfg))
#   
#   if (log(runif(1)) < log_acc) {
#     prop$last_accept <- TRUE
#     return(prop)
#   } else {
#     p$last_accept <- FALSE
#     return(p)
#   }
# }



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
  BB7   <- FAM_BB7
  BB7s  <- FAM_BB7R180
  
  
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
    p_bb7   <- sum(w_new[fam_mat[, e] == BB7])
    p_bb7r  <- sum(w_new[fam_mat[, e] == BB7s])
    
    res <- list(edge   = e,
                p_indep = p_indep,
                p_gauss = p_gauss,
                p_bb1   = p_bb1,
                p_bb1r  = p_bb1r,
                p_bb8r  = p_bb8r)
    
    if (print_flag)
      if (p_indep > 0){
      cat(sprintf("     Edge %2d  Indep   : P=%.3f\n", e, p_indep))
      }
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
    
    ## BB7 (Joe–Clayton) parameter stats
    if (p_bb7 > 0) {
      mask  <- fam_mat[, e] == BB7
      th1_e <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      wloc  <- w_new[mask] / p_bb7
      stats <- function(x){
        mu <- sum(wloc * x)
        sd <- sqrt(sum(wloc * (x - mu)^2))
        qs <- w_quantile(x, wloc, q_probs)
        c(mu, sd, qs)
      }
      s1 <- stats(th1_e); s2 <- stats(th2_e)
      res <- c(res,
               list(mu_bb7_theta = s1[1], sd_bb7_theta = s1[2],
                    bb7_theta_q025 = s1[3], bb7_theta_q975 = s1[4],
                    mu_bb7_delta = s2[1], sd_bb7_delta = s2[2],
                    bb7_delta_q025 = s2[3], bb7_delta_q975 = s2[4]))
      
      if (print_flag)
        cat(sprintf(
          "                BB7     : P=%.3f | θ=%.2f±%.2f CI=[%.2f,%.2f] | δ=%.2f±%.2f CI=[%.2f,%.2f]\n",
          p_bb7,
          s1[1], s1[2], s1[3], s1[4],
          s2[1], s2[2], s2[3], s2[4]
        ))
    }
    
    ## BB7^180 (rotated) parameter stats
    if (p_bb7r > 0) {
      mask  <- fam_mat[, e] == BB7s
      th1_e <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      wloc  <- w_new[mask] / p_bb7r
      stats <- function(x){
        mu <- sum(wloc * x)
        sd <- sqrt(sum(wloc * (x - mu)^2))
        qs <- w_quantile(x, wloc, q_probs)
        c(mu, sd, qs)
      }
      s1 <- stats(th1_e); s2 <- stats(th2_e)
      res <- c(res,
               list(mu_bb7r_theta = s1[1], sd_bb7r_theta = s1[2],
                    bb7r_theta_q025 = s1[3], bb7r_theta_q975 = s1[4],
                    mu_bb7r_delta = s2[1], sd_bb7r_delta = s2[2],
                    bb7r_delta_q025 = s2[3], bb7r_delta_q975 = s2[4]))
      
      if (print_flag)
        cat(sprintf(
          "                BB7^180 : P=%.3f | θ=%.2f±%.2f CI=[%.2f,%.2f] | δ=%.2f±%.2f CI=[%.2f,%.2f]\n",
          p_bb7r,
          s1[1], s1[2], s1[3], s1[4],
          s2[1], s2[2], s2[3], s2[4]
        ))
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
















