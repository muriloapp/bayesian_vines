library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(matrixStats)


FAM_INFO <- data.frame(
  name      = c("indep","gaussian","bb1","bb1r180","bb8r180","bb7","bb7r180","t"),
  rv_name   = c("indep","gaussian","bb1","bb1","bb8","bb7","bb7","t"),
  rotation  = c(0L, 0L, 0L, 180L, 180L, 0L, 180L, 0L),
  code      = c(0L, 1L, 3L, 13L, 16L, 7L, 17L, 2L),  # any unique ints are fine internally
  npar      = c(0L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
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
FAM_T <- FAM_INFO$code[FAM_INFO$name == "t"]

## Minimal templates (we’ll build bicops from rv_name/rotation on the fly too)
T_INDEP     <- bicop_dist("indep")
T_GAUSS     <- bicop_dist("gaussian", parameters = 0)

T_BB1       <- bicop_dist("bb1", rotation = 0,   parameters = c(1, 2))
T_BB1R180   <- bicop_dist("bb1", rotation = 180, parameters = c(1, 2))

T_BB8R180   <- bicop_dist("bb8", rotation = 180, parameters = c(1, 1))  # dummy two-par

T_BB7       <- bicop_dist("bb7", rotation = 0,   parameters = c(1, 2))   # NEW
T_BB7R180   <- bicop_dist("bb7", rotation = 180, parameters = c(1, 2))   # NEW

T_T <- bicop_dist("t", rotation = 0, parameters = c(0, 5))  # rho=0, nu=5 dummy



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
  
  if (fam == FAM_T) {
    # map to tail, add λ prior + Jacobian + df prior
    lam <- t_par2tail(th1, th2)
    lp_tail <- dbeta(lam, 2, 2, log = TRUE)           # diffuse
    lp_jac  <- t_log_jacobian(lam, th2)
    lp_df <- if (th2 >= 2 && th2 <= 30) -log(th2) - log(log(30)) else -Inf
    return(fillna_neg(lp_tail + lp_jac + lp_df))      # - cfg$lambda if you re-enable penalties
  }
  
  
  stop("Unknown family code in log_prior_edge: ", fam)
}



log_prior <- function(p, cfg) {
  sum(mapply(log_prior_edge,
             p$fam, p$th1, p$th2,
             MoreArgs = list(cfg = cfg)))
}


safe_sample <- function(x, size, replace = FALSE, prob = NULL) {
  # same call shape as base::sample
  if (size != 1L) {
    # fall back to base for any size != 1
    return(base::sample(x, size = size, replace = replace, prob = prob))
  }
  n <- length(x)
  if (n == 0L) stop("safe_sample: no candidates to sample from")
  if (n == 1L) return(x[[1L]])  # short-circuit RNG
  
  # sanitize probabilities (accept NULL, wrong length, non-finite, all zeros)
  p <- prob
  if (!is.null(p)) {
    p <- as.numeric(p)
    if (length(p) != n || any(!is.finite(p))) p <- rep(1, n)
    p <- pmax(p, 0)
    s <- sum(p)
    if (s <= 0) p <- rep(1, n) else p <- p / s
  }
  base::sample(x, size = 1L, replace = replace, prob = p)
}



# returns a single particle as 1×K row matrices (so can rbind into big mats)
new_particle_row_mat <- function(cfg, U_init = NULL) {
  K   <- cfg$K
  fam <- integer(K); th1 <- numeric(K); th2 <- numeric(K)
  
  tip_on  <- isTRUE(cfg$use_tail_informed_prior) && !is.null(U_init)
  if (tip_on && is.null(cfg$edge_pair))
    stop("TIP init requested but cfg$edge_pair is NULL. Provide Tree-1 pairs aligned with edge order.")
  
  for (k in seq_len(K)) {
    tr_k <- cfg$edge_tree[k]
    allowed_names <- if (tr_k == 1L) cfg$families_first else cfg$families_deep
    fam_tbl <- FAM_INFO[FAM_INFO$name %in% allowed_names, , drop = FALSE]
    
    if (tip_on && tr_k == 1L) {
      # ---- TIP-aware family selection on Tree 1 ----
      local_idx <- sum(cfg$edge_tree[seq_len(k)] == 1L)
      pair      <- cfg$edge_pair[local_idx, ]
      uv        <- U_init[, pair, drop = FALSE]
      
      emps <- emp_tails_FRAPO(uv, method = cfg$tip_method, k = cfg$tip_k)
      mL   <- as.numeric(emps["L"]);  mU <- as.numeric(emps["U"])
      tail_strength <- max(mL, mU)
      
      #tail_codes <- c(FAM_BB1, FAM_BB1R180, FAM_BB7, FAM_BB7R180, FAM_BB8R180, FAM_T)
      #base_w     <- rep(1.0, nrow(fam_tbl))
      #is_tailfam <- fam_tbl$code %in% tail_codes
      #w          <- base_w * ifelse(is_tailfam, 1 + 3*tail_strength, pmax(1 - 0.9*tail_strength, 0.1))
      #w          <- w / sum(w)
      w          <- rep(1 / length(cfg$families_first), length(cfg$families_first))
      
      code_k  <- sample(fam_tbl$code, 1L, prob = w)
      fam[k]  <- code_k
      
      if (code_k == FAM_INDEP) {
        th1[k] <- 0; th2[k] <- 0
      } else if (code_k == FAM_GAUSS) {
        th1[k] <- runif(1, -0.99, 0.99); th2[k] <- 0
      } else {
        th_new <- seed_family_from_emp(
          new_code     = code_k,
          mL           = mL,
          mU           = mU,
          cur_delta    = runif(1, 1.1, 5.0),
          step_sd      = cfg$step_sd,
          tip_sd_logit = cfg$tip_sd_logit
        )
        th1[k] <- th_new[1]; th2[k] <- th_new[2]
      }
      
    } else {
      # ---- original selection / seeding ----
      weights <- rep(1 / nrow(fam_tbl), nrow(fam_tbl))
      code_k  <- safe_sample(fam_tbl$code, 1L, prob = weights)
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
        
      } else if (code_k == FAM_T) {
        nu  <- exp(runif(1, log(2), log(30)))   # log-uniform
        lam <- rbeta(1, 2, 2)
        rho <- t_tail2rho(lam, nu)
        pr  <- sanitize_t(rho, nu)
        th1[k] <- pr[1]; th2[k] <- pr[2]
      }
    }
  }
  
  list(
    fam = matrix(fam, nrow = 1L),
    th1 = matrix(th1, nrow = 1L),
    th2 = matrix(th2, nrow = 1L),
    w   = 1 / cfg$M,
    last_accept = FALSE
  )
}

# Build M particles and stack row-wise: fast to serialize, easy to slice by row
new_particles_mats <- function(cfg, U_init = NULL) {
  M <- as.integer(cfg$M)
  rows <- vector("list", M)
  for (m in seq_len(M)) {
    rows[[m]] <- new_particle_row_mat(cfg, U_init)  # must return 1×K matrices
  }
  
  fam_mat <- do.call(rbind, lapply(rows, `[[`, "fam"))
  th1_mat <- do.call(rbind, lapply(rows, `[[`, "th1"))
  th2_mat <- do.call(rbind, lapply(rows, `[[`, "th2"))
  
  storage.mode(fam_mat) <- "integer"  # <-- fix: mutate in place
  
  w_vec   <- rep(1 / M, M)
  acc_vec <- rep(FALSE, M)
  
  list(
    fam_mat     = fam_mat,
    th1_mat     = th1_mat,
    th2_mat     = th2_mat,
    w           = w_vec,
    last_accept = acc_vec
  )
}

# 
# 
# fast_vine_from_particle <- function(p, skel) {
#   pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)
#   idx <- 1L
#   
#   for (tr in seq_along(pcs)) {
#     for (ed in seq_along(pcs[[tr]])) {
#       code <- p$fam[idx]
#       spec <- fam_spec(code)
#       
#       if (code == FAM_INDEP) {
#         pcs[[tr]][[ed]] <- T_INDEP
#         
#       } else if (code == FAM_GAUSS) {
#         bc <- rlang::duplicate(T_GAUSS, shallow = TRUE)
#         bc$parameters[1] <- p$th1[idx]
#         pcs[[tr]][[ed]] <- bc
#         
#       } else if (code %in% c(FAM_BB1, FAM_BB1R180)) {
#         pars <- sanitize_bb1(p$th1[idx], p$th2[idx])
#         p$th1[idx] <- pars[1]; p$th2[idx] <- pars[2]
#         
#         tmpl <- if (code == FAM_BB1) T_BB1 else T_BB1R180
#         bc <- rlang::duplicate(tmpl, shallow = TRUE)
#         bc$parameters[1:2] <- pars
#         pcs[[tr]][[ed]] <- bc
#         
#       } else if (code == FAM_BB8R180) {
#         bc <- rlang::duplicate(T_BB8R180, shallow = TRUE)
#         bc$parameters[1:2] <- c(p$th1[idx], p$th2[idx])
#         pcs[[tr]][[ed]] <- bc
#         
#       } else if (code == FAM_BB7) {                       # NEW
#         bc <- rlang::duplicate(T_BB7, shallow = TRUE)
#         bc$parameters[1:2] <- c(p$th1[idx], p$th2[idx])
#         pcs[[tr]][[ed]] <- bc
#         
#       } else if (code == FAM_BB7R180) {                   # NEW
#         bc <- rlang::duplicate(T_BB7R180, shallow = TRUE)
#         bc$parameters[1:2] <- c(p$th1[idx], p$th2[idx])
#         pcs[[tr]][[ed]] <- bc
#         
#       } else if (code == FAM_T) {
#         bc <- rlang::duplicate(T_T, shallow = TRUE)
#         bc$parameters[1:2] <- c(p$th1[idx], p$th2[idx])  # rho, nu
#         pcs[[tr]][[ed]] <- bc
#         
#       } else {
#         stop("fast_vine_from_particle: unsupported family code: ", code)
#       }
#       
#       idx <- idx + 1L
#     }
#   }
#   
#   vinecop_dist(pcs, skel$structure)
# }

## ---- 0) Small helper: count edges in skeleton (K must match) ----
K_of_skeleton <- function(skel) {
  # number of pair-copulas across all trees
  sum(vapply(skel$pair_copulas, length, integer(1)))
}

## ---- 1) Coerce any particle shape into simple vectors ----
.as_particle_vectors <- function(p) {
  # Supports:
  #   - old: p$fam, p$th1, p$th2 as atomic vectors (length K)
  #   - new: p$fam, p$th1, p$th2 as 1×K matrices (row)
  fam <- if (is.matrix(p$fam)) as.integer(p$fam[1, ]) else as.integer(p$fam)
  th1 <- if (is.matrix(p$th1)) as.numeric(p$th1[1, ]) else as.numeric(p$th1)
  th2 <- if (is.matrix(p$th2)) as.numeric(p$th2[1, ]) else as.numeric(p$th2)
  list(fam = fam, th1 = th1, th2 = th2)
}

## ---- 2) Core builder: from fam/th1/th2 vectors to a vine ----
.build_vine_from_vectors <- function(fam, th1, th2, skel) {
  K <- K_of_skeleton(skel)
  if (length(fam) != K || length(th1) != K || length(th2) != K) {
    stop(sprintf("fast vine build: K mismatch (fam/th1/th2 lengths %d/%d/%d vs K=%d)",
                 length(fam), length(th1), length(th2), K))
  }
  
  pcs <- rlang::duplicate(skel$pair_copulas, shallow = TRUE)
  idx <- 1L
  
  for (tr in seq_along(pcs)) {
    for (ed in seq_along(pcs[[tr]])) {
      
      code <- fam[idx]
      
      if (code == FAM_INDEP) {
        pcs[[tr]][[ed]] <- T_INDEP
        
      } else if (code == FAM_GAUSS) {
        bc <- rlang::duplicate(T_GAUSS, shallow = TRUE)
        bc$parameters[1] <- th1[idx]
        pcs[[tr]][[ed]] <- bc
        
      } else if (code %in% c(FAM_BB1, FAM_BB1R180)) {
        pars <- sanitize_bb1(th1[idx], th2[idx])
        tmpl <- if (code == FAM_BB1) T_BB1 else T_BB1R180
        bc <- rlang::duplicate(tmpl, shallow = TRUE)
        bc$parameters[1:2] <- pars
        pcs[[tr]][[ed]] <- bc
        
      } else if (code == FAM_BB8R180) {
        bc <- rlang::duplicate(T_BB8R180, shallow = TRUE)
        bc$parameters[1:2] <- c(th1[idx], th2[idx])
        pcs[[tr]][[ed]] <- bc
        
      } else if (code == FAM_BB7) {
        bc <- rlang::duplicate(T_BB7, shallow = TRUE)
        bc$parameters[1:2] <- c(th1[idx], th2[idx])
        pcs[[tr]][[ed]] <- bc
        
      } else if (code == FAM_BB7R180) {
        bc <- rlang::duplicate(T_BB7R180, shallow = TRUE)
        bc$parameters[1:2] <- c(th1[idx], th2[idx])
        pcs[[tr]][[ed]] <- bc
        
      } else if (code == FAM_T) {
        bc <- rlang::duplicate(T_T, shallow = TRUE)
        bc$parameters[1:2] <- c(th1[idx], th2[idx])  # rho, nu
        pcs[[tr]][[ed]] <- bc
        
      } else {
        stop("fast_vine_from_particle: unsupported family code: ", code)
      }
      
      idx <- idx + 1L
    }
  }
  
  vinecop_dist(pcs, skel$structure)
}


## ---- 3B) Fast path for matrix storage: take row m directly ----
fast_vine_from_row <- function(mats, m, skel) {
  fam <- as.integer(mats$fam_mat[m, ])
  th1 <- as.numeric(mats$th1_mat[m, ])
  th2 <- as.numeric(mats$th2_mat[m, ])
  .build_vine_from_vectors(fam, th1, th2, skel)
}


# Guarded single-row density → log
.safe_logdens1 <- function(u_row, vine) {
  U1 <- if (is.matrix(u_row)) u_row else matrix(u_row, nrow = 1L)
  d  <- dvinecop(U1, vine, cores = 1)
  if (!is.finite(d) || d <= 0) -1e100 else log(d)
}

compute_log_incr <- function(particles, u_row, skeleton) {
  M <- nrow(particles$fam_mat)
  log_incr <- numeric(M)
  
  for (j in seq_len(M)) {
    vine_j      <- fast_vine_from_row(particles, j, skeleton)
    log_incr[j] <- .safe_logdens1(u_row, vine_j)
  }
  
  log_incr
}


update_weights <- function(particles, log_lik) {
  w_prev <- particles$w
  if (length(w_prev) != length(log_lik))
    stop(sprintf("length mismatch: w=%d, log_lik=%d",
                 length(w_prev), length(log_lik)))
  
  # log-sum-exp normalization for numerical stability
  log_w <- log(pmax(w_prev, .Machine$double.xmin)) + log_lik
  log_w <- log_w - max(log_w)
  w_new <- exp(log_w)
  s     <- sum(w_new)
  if (!is.finite(s) || s <= 0) stop("weights under/overflow")
  w_new <- w_new / s
  
  particles$w <- w_new
  attr(particles, "ESS") <- ESS(w_new)   # handy to read right after
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
  if (code == FAM_T) { v <- t_par2tail(th1, th2); return(list(L = v, U = v)) }
  
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
  if (new_code == FAM_T) {
    lam <- clamp01(if (!is.null(tails$L)) tails$L else if (!is.null(tails$U)) tails$U else rbeta(1,2,2))
    nu <- exp(runif(1, log(2), log(30)))             # log-uniform
    rho <- t_tail2rho(lam, nu)
    pr  <- sanitize_t(rho, nu)
    return(pr)
  }
  
  stop("Unknown family code in mh_step")
}
# 
# mh_step <- function(p, data_up_to_t, skeleton, cfg, tip_means_cache = NULL,   # <-- NEW
#                     wobs = NULL)              # (optional, unchanged logic)) 
#   {
#   
#   K    <- cfg$K
#   prop <- p
#   
#   ## (a) Family flips — tail-seeded init (Tree 1 only; no MLE)
#   end_first_tr <- sum(cfg$edge_tree == 1L)  # number of edges in Tree 1
#   if (end_first_tr > 0L) {
#     flip_mask <- runif(end_first_tr) < cfg$q_flip
#     if (any(flip_mask)) {
#       idxs <- which(flip_mask)
#       for (idx in idxs) {
#         tr_idx <- cfg$edge_tree[idx]
#         if (tr_idx != 1L) next  # safety; we only flip on Tree 1
#         
#         old_code <- prop$fam[idx]
#         allowed_names <- if (tr_idx == 1L) cfg$families_first else cfg$families_deep
#         allowed_codes <- setdiff(FAM_INFO$code[FAM_INFO$name %in% allowed_names], old_code)
#         if (!length(allowed_codes)) next
#         
#         new_code <- sample(allowed_codes, 1L)
#         prop$fam[idx] <- new_code
#         
#         if (isTRUE(cfg$use_tail_informed_prior)) {
#           tipm <- if (is.null(tip_means_cache)) NULL else tip_means_cache[[idx]]
#           if (is.null(tipm)) {
#             tails_old <- get_tails(old_code, prop$th1[idx], prop$th2[idx])
#             tp_new    <- init_from_tails(new_code, tails_old)
#             prop$th1[idx] <- tp_new[1]; prop$th2[idx] <- tp_new[2]
#           } else {
#             th_new <- seed_family_from_emp(
#               new_code,
#               mL           = tipm["mL"],
#               mU           = tipm["mU"],
#               cur_delta    = prop$th2[idx],
#               step_sd      = cfg$step_sd,
#               tip_sd_logit = cfg$tip_sd_logit
#             )
#             prop$th1[idx] <- th_new[1]; prop$th2[idx] <- th_new[2]
#           }
#         } else {
#           tails_old <- get_tails(old_code, prop$th1[idx], prop$th2[idx])
#           tp_new    <- init_from_tails(new_code, tails_old)
#           prop$th1[idx] <- tp_new[1]; prop$th2[idx] <- tp_new[2]
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
#   # t-copula — RW in λ (tail) and z=log(nu), then map back to (rho, nu)
#   t_idx <- which(prop$fam == FAM_T)
#   if (length(t_idx)) {
#     lam <- vapply(t_idx, function(i) t_par2tail(prop$th1[i], prop$th2[i]), numeric(1))
#     lam <- rtnorm_vec(lam, cfg$step_sd, 1e-3, 0.99)
#     
#     z  <- log(pmax(prop$th2[t_idx], 2))                 # enforce ≥2
#     z  <- rtnorm_vec(z, cfg$step_sd, log(2), log(30))   # truncated normal in [log 2, log 30]
#     nu <- exp(z)
#     
#     pars <- t(vapply(seq_along(t_idx), function(j) {
#       i   <- t_idx[j]
#       rho <- t_tail2rho(lam[j], nu[j])
#       sanitize_t(rho, nu[j])
#     }, numeric(2)))
#     prop$th1[t_idx] <- pars[,1]; prop$th2[t_idx] <- pars[,2]
#   }
#   
#   
#   
#   
#   
#   ## (c) MH ratio
#   vine_prop <- fast_vine_from_particle(prop, skeleton)
#   vine_curr <- fast_vine_from_particle(p,    skeleton)
#   
#   #ll_prop <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_prop), .Machine$double.eps))))
#   #ll_curr <- sum(fillna_neg(log(pmax(dvinecop(data_up_to_t, vine_curr), .Machine$double.eps))))
# 
#   if (!is.null(wobs)) {
#     ll_prop <- sum(wobs * safe_logdens(data_up_to_t, vine_prop))
#     ll_curr <- sum(wobs * safe_logdens(data_up_to_t, vine_curr))
#   } else {
#     ll_prop <- sum(safe_logdens(data_up_to_t, vine_prop))
#     ll_curr <- sum(safe_logdens(data_up_to_t, vine_curr))
#   }
#   
#   lp_prop <- if (isTRUE(cfg$use_tail_informed_prior))
#     log_prior_with_tip_cached(prop, cfg, tip_means_cache) else
#       log_prior(prop, cfg)
#   
#   lp_curr <- if (isTRUE(cfg$use_tail_informed_prior))
#     log_prior_with_tip_cached(p,    cfg, tip_means_cache) else
#       log_prior(p,    cfg)
#   
#   log_acc <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
#   
#   if (log(runif(1)) < log_acc) { prop$last_accept <- TRUE; return(prop) }
#   p$last_accept <- FALSE; return(p)
# }
# 






resample_move_old <- function(particles, newAncestors,
                              data_up_to_t, cl, type, cfg,
                              skeleton = NULL, tr = NULL, temp_skel = NULL) {
  
  # ----- 1) RESAMPLE (row reindex) -----
  M <- nrow(particles$fam_mat)
  particles$fam_mat <- particles$fam_mat[newAncestors, , drop = FALSE]
  particles$th1_mat <- particles$th1_mat[newAncestors, , drop = FALSE]
  particles$th2_mat <- particles$th2_mat[newAncestors, , drop = FALSE]
  particles$w       <- rep(1 / M, M)
  if (!is.null(particles$last_accept)) particles$last_accept[] <- FALSE
  
  # ----- 2) PRECOMPUTE (once) -----
  tip_cache <- if (isTRUE(cfg$use_tail_informed_prior))
    precompute_tip_means(data_up_to_t, cfg) else NULL
  wobs <- if (isTRUE(cfg$use_weighted_ll))
    tail_weights(data_up_to_t, tau = cfg$tauL,
                 k = cfg$joint_k, eps = cfg$tail_eps) else NULL
  
  # ship big read-only objects once
  parallel::clusterExport(cl, c("particles","tip_cache","wobs"), envir = environment())
  
  mh_n_prop <- M * cfg$n_mh
  tic("MH move")
  
  if (identical(type, "standard")) {
    mh_results <- parallel::parLapply(
      cl, seq_len(M),
      function(i, data_up_to_t, skeleton, cfg, tip_cache, wobs) {
        fam <- as.integer(particles$fam_mat[i, ])
        th1 <- as.numeric(particles$th1_mat[i, ])
        th2 <- as.numeric(particles$th2_mat[i, ])
        
        acc_local <- 0L
        for (k in seq_len(cfg$n_mh)) {
          res <- mh_step(fam, th1, th2, data_up_to_t, skeleton, cfg,
                             tip_means_cache = tip_cache, wobs = wobs)
          fam <- res$fam; th1 <- res$th1; th2 <- res$th2
          if (isTRUE(res$last_accept)) acc_local <- acc_local + 1L
        }
        list(fam = fam, th1 = th1, th2 = th2, acc = acc_local)
      },
      data_up_to_t, skeleton, cfg, tip_cache, wobs
    )
  } else {
    stop("Only type='standard' is implemented for the matrix-based path.")
  }
  
  toc()
  
  # ----- 3) COLLECT -----
  mh_n_acc <- sum(vapply(mh_results, `[[`, integer(1), "acc"))
  for (i in seq_len(M)) {
    particles$fam_mat[i, ] <- as.integer(mh_results[[i]]$fam)
    particles$th1_mat[i, ] <- mh_results[[i]]$th1
    particles$th2_mat[i, ] <- mh_results[[i]]$th2
  }
  
  acc_pct <- 100 * mh_n_acc / mh_n_prop
  cat(sprintf("MH acceptance = %4d / %4d  =  %.2f%%\n\n",
              mh_n_acc, mh_n_prop, acc_pct))
  
  list(particles = particles, acc_pct = acc_pct)
}


mh_step <- function(fam, th1, th2,
                        data_up_to_t, skeleton, cfg,
                        tip_means_cache = NULL,   # list length K, each NULL or c(mL=..,mU=..)
                        wobs = NULL) {           # optional obs weights
  K <- length(fam)
  
  # --- propose ---
  fam_p  <- fam
  th1_p  <- th1
  th2_p  <- th2
  
  ## (a) Family flips — Tree 1 only (tail-seeded init if TIP)
  end_first_tr <- sum(cfg$edge_tree == 1L)
  if (end_first_tr > 0L) {
    flip_mask <- runif(end_first_tr) < cfg$q_flip
    if (any(flip_mask)) {
      idxs <- which(flip_mask)
      for (idx in idxs) {
        tr_idx <- cfg$edge_tree[idx]
        if (tr_idx != 1L) next
        
        old_code <- fam_p[idx]
        allowed_names <- if (tr_idx == 1L) cfg$families_first else cfg$families_deep
        allowed_codes <- setdiff(FAM_INFO$code[FAM_INFO$name %in% allowed_names], old_code)
        if (!length(allowed_codes)) next
        
        new_code <- safe_sample(allowed_codes, 1L)
        fam_p[idx] <- new_code
        
        if (isTRUE(cfg$use_tail_informed_prior)) {
          tipm <- if (is.null(tip_means_cache)) NULL else tip_means_cache[[idx]]
          if (is.null(tipm)) {
            tails_old <- get_tails(old_code, th1_p[idx], th2_p[idx])
            tp_new    <- init_from_tails(new_code, tails_old)
            th1_p[idx] <- tp_new[1]; th2_p[idx] <- tp_new[2]
          } else {
            th_new <- seed_family_from_emp(
              new_code,
              mL           = tipm["mL"],
              mU           = tipm["mU"],
              cur_delta    = th2_p[idx],
              step_sd      = cfg$step_sd,
              tip_sd_logit = cfg$tip_sd_logit
            )
            th1_p[idx] <- th_new[1]; th2_p[idx] <- th_new[2]
          }
        } else {
          tails_old <- get_tails(old_code, th1_p[idx], th2_p[idx])
          tp_new    <- init_from_tails(new_code, tails_old)
          th1_p[idx] <- tp_new[1]; th2_p[idx] <- tp_new[2]
        }
      }
    }
  }
  
  ## (b) Random-walk in tail space, map back, sanitize
  # Gaussian
  ga_idx <- which(fam_p == FAM_GAUSS)
  if (length(ga_idx))
    th1_p[ga_idx] <- rtnorm_vec(th1_p[ga_idx], cfg$step_sd, -0.99, 0.99)
  
  # BB1
  bb1_idx <- which(fam_p == FAM_BB1)
  if (length(bb1_idx)) {
    lam <- t(mapply(bb1_par2tail, th1_p[bb1_idx], th2_p[bb1_idx]))
    lam[,1] <- rtnorm_vec(lam[,1], cfg$step_sd, 1e-3, 0.99)
    lam[,2] <- rtnorm_vec(lam[,2], cfg$step_sd, 1e-3, 0.99)
    par <- t(apply(lam, 1, \(x) bb1_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
    th1_p[bb1_idx] <- par[,1]; th2_p[bb1_idx] <- par[,2]
  }
  
  # BB1^180
  bb1r_idx <- which(fam_p == FAM_BB1R180)
  if (length(bb1r_idx)) {
    lamr <- t(mapply(bb1r180_par2tail, th1_p[bb1r_idx], th2_p[bb1r_idx]))
    lamr[,1] <- rtnorm_vec(lamr[,1], cfg$step_sd, 1e-3, 0.99)
    lamr[,2] <- rtnorm_vec(lamr[,2], cfg$step_sd, 1e-3, 0.99)
    par <- t(apply(lamr, 1, \(x) bb1r180_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb1(x[1], x[2])))
    th1_p[bb1r_idx] <- par[,1]; th2_p[bb1r_idx] <- par[,2]
  }
  
  # BB8^180
  bb8r_idx <- which(fam_p == FAM_BB8R180)
  if (length(bb8r_idx)) {
    lamLr <- vapply(bb8r_idx, function(i) bb8r180_par2tail(th1_p[i], th2_p[i])[1], numeric(1))
    lamLr <- rtnorm_vec(lamLr, cfg$step_sd, 1e-3, 0.99)
    pars  <- t(vapply(seq_along(bb8r_idx), function(j) {
      i <- bb8r_idx[j]
      delta_prop <- pmin(pmax(th2_p[i] + rnorm(1, 0, cfg$step_sd), 0.05), 7.0)
      bb8r180_tail2par(lamLr[j], delta_prop)
    }, numeric(2)))
    pars <- t(apply(pars, 1, \(x) sanitize_bb8(x[1], x[2])))
    th1_p[bb8r_idx] <- pars[,1]; th2_p[bb8r_idx] <- pars[,2]
  }
  
  # BB7
  bb7_idx <- which(fam_p == FAM_BB7)
  if (length(bb7_idx)) {
    lam <- t(mapply(bb7_par2tail, th1_p[bb7_idx], th2_p[bb7_idx]))
    lam[,1] <- rtnorm_vec(lam[,1], cfg$step_sd, 1e-3, 0.99)
    lam[,2] <- rtnorm_vec(lam[,2], cfg$step_sd, 1e-3, 0.99)
    par <- t(apply(lam, 1, \(x) bb7_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb7(x[1], x[2])))
    th1_p[bb7_idx] <- par[,1]; th2_p[bb7_idx] <- par[,2]
  }
  
  # BB7^180
  bb7r_idx <- which(fam_p == FAM_BB7R180)
  if (length(bb7r_idx)) {
    lamr <- t(mapply(bb7r180_par2tail, th1_p[bb7r_idx], th2_p[bb7r_idx]))
    lamr[,1] <- rtnorm_vec(lamr[,1], cfg$step_sd, 1e-3, 0.99)
    lamr[,2] <- rtnorm_vec(lamr[,2], cfg$step_sd, 1e-3, 0.99)
    par <- t(apply(lamr, 1, \(x) bb7r180_tail2par(x[1], x[2])))
    par <- t(apply(par, 1, \(x) sanitize_bb7(x[1], x[2])))
    th1_p[bb7r_idx] <- par[,1]; th2_p[bb7r_idx] <- par[,2]
  }
  
  # t-copula
  t_idx <- which(fam_p == FAM_T)
  if (length(t_idx)) {
    lam <- vapply(t_idx, function(i) t_par2tail(th1_p[i], th2_p[i]), numeric(1))
    lam <- rtnorm_vec(lam, cfg$step_sd, 1e-3, 0.99)
    z   <- log(pmax(th2_p[t_idx], 2))
    z   <- rtnorm_vec(z, cfg$step_sd, log(2), log(30))
    nu  <- exp(z)
    pars <- t(vapply(seq_along(t_idx), function(j) {
      i   <- t_idx[j]
      rho <- t_tail2rho(lam[j], nu[j])
      sanitize_t(rho, nu[j])
    }, numeric(2)))
    th1_p[t_idx] <- pars[,1]; th2_p[t_idx] <- pars[,2]
  }
  
  ## (c) MH ratio
  vine_prop <- .build_vine_from_vectors(fam_p, th1_p, th2_p, skeleton)
  vine_curr <- .build_vine_from_vectors(fam,   th1,   th2,   skeleton)
  
  if (is.null(wobs)) {
    ll_prop <- sum(safe_logdens(data_up_to_t, vine_prop))
    ll_curr <- sum(safe_logdens(data_up_to_t, vine_curr))
  } else {
    ll_prop <- sum(wobs * safe_logdens(data_up_to_t, vine_prop))
    ll_curr <- sum(wobs * safe_logdens(data_up_to_t, vine_curr))
  }
  
  # reuse existing prior fns via tiny wrappers
  p_prop <- list(fam = fam_p, th1 = th1_p, th2 = th2_p)
  p_curr <- list(fam = fam,   th1 = th1,   th2 = th2)
  
  lp_prop <- if (isTRUE(cfg$use_tail_informed_prior))
    log_prior_with_tip_cached(p_prop, cfg, tip_means_cache) else
      log_prior(p_prop, cfg)
  lp_curr <- if (isTRUE(cfg$use_tail_informed_prior))
    log_prior_with_tip_cached(p_curr, cfg, tip_means_cache) else
      log_prior(p_curr, cfg)
  
  log_acc <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
  
  if (log(runif(1)) < log_acc) {
    return(list(fam = fam_p, th1 = th1_p, th2 = th2_p, last_accept = TRUE))
  } else {
    return(list(fam = fam,   th1 = th1,   th2 = th2,   last_accept = FALSE))
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
    z <- rnorm(1, mean = logit(clamp01(m)), sd = step_sd) #tip_sd_logit
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
  #if (new_code == FAM_BB8R180) {
  #  lamLr <- r_lam(mL)  # choose L as anchor
  #  delta <- pmin(pmax(cur_delta + rnorm(1, 0, step_sd), 0.05), 7.0)
  #  th    <- bb8r180_tail2par(lamLr, delta); th <- sanitize_bb8(th[1], th[2]); return(th)
  #}
  if (new_code == FAM_T) {
    lam <- {
      z <- rnorm(1, mean = logit(clamp01(max(mL, mU))), sd = step_sd) #tip_sd_logit
      ilogit(z)
    }
    nu  <- exp(runif(1, 0, log(30)))
    rho <- t_tail2rho(lam, nu)
    pr  <- sanitize_t(rho, nu)
    return(pr)
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
      dlogitnorm(lam[1], logit(tip_means["mL.L"]), s) +
      dlogitnorm(lam[2], logit(tip_means["mU.U"]), s)
    return(tip + bb1_log_jacobian(lam[1], lam[2]))
  }
  if (fam == FAM_BB7) {
    lam <- bb7_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mL.L"]), s) +
      dlogitnorm(lam[2], logit(tip_means["mU.U"]), s)
    return(tip + bb7_log_jacobian(lam[1], lam[2]) )
  }
  if (fam == FAM_BB1R180) {
    lam <- bb1r180_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mU.U"]), s) +  # rotated mapping
      dlogitnorm(lam[2], logit(tip_means["mL.L"]), s)
    return(tip + bb1_log_jacobian(lam[1], lam[2]))
  }
  if (fam == FAM_BB7R180) {
    lam <- bb7r180_par2tail(th1, th2)
    tip <- if (is.null(tip_means)) 0 else
      dlogitnorm(lam[1], logit(tip_means["mU.U"]), s) +
      dlogitnorm(lam[2], logit(tip_means["mL.L"]), s)
    return(tip + bb7_log_jacobian(lam[1], lam[2]))
  }
  if (fam == FAM_BB8R180) {
    lamLr <- bb8r180_par2tail(th1, th2)[1]
    tip   <- if (is.null(tip_means)) 0 else dlogitnorm(lamLr, logit(tip_means["mL.L"]), s)
    return(tip + bb8_log_jacobian(lamLr))
  }
  stop("unknown family in log_prior_edge_strong")
}

## TIP-aware prior using precomputed tip_means_cache (list length K).
## If tip_means_cache[[e]] is NULL (not Tree 1 or no TIP), fall back to diffuse Beta tail priors.
log_prior_with_tip_cached <- function(p, cfg, tip_means_cache) {
  lp <- 0
  s  <- cfg$tip_sd_logit
  for (e in seq_len(cfg$K)) {
    fam <- p$fam[e]; th1 <- p$th1[e]; th2 <- p$th2[e]
    tipm <- tip_means_cache[[e]]  # NULL if not Tree 1 (or TIP off)
    
    ## match your constants (adjust these two lines if you don't want the -cfg$lambda penalty)
    if (fam == FAM_INDEP) { lp <- lp + (0 - cfg$lambda); next }
    if (fam == FAM_GAUSS) { lp <- lp + (log(1/1.98) - cfg$lambda); next }
    
    if (fam == FAM_BB1) {
      lam <- bb1_par2tail(th1, th2)   # (λL, λU)
      tip_or_diffuse <-
        if (is.null(tipm)) {
          dbeta(lam[1], 2, 2, log=TRUE) + dbeta(lam[2], 2, 2, log=TRUE)
        } else {
          dlogitnorm(lam[1], logit(tipm["mL"]), s) +
            dlogitnorm(lam[2], logit(tipm["mU"]), s)
        }
      lp <- lp + tip_or_diffuse + bb1_log_jacobian(lam[1], lam[2]); next
    }
    
    if (fam == FAM_BB1R180) {
      lamr <- bb1r180_par2tail(th1, th2)   # (λL^rot, λU^rot)
      tip_or_diffuse <-
        if (is.null(tipm)) {
          dbeta(lamr[1], 2, 2, log=TRUE) + dbeta(lamr[2], 2, 2, log=TRUE)
        } else {
          ## rotation: lower-rot ↔ upper(orig), upper-rot ↔ lower(orig)
          dlogitnorm(lamr[1], logit(tipm["mU"]), s) +
            dlogitnorm(lamr[2], logit(tipm["mL"]), s)
        }
      lp <- lp + tip_or_diffuse + bb1r180_log_jacobian(lamr[1], lamr[2]); next
    }
    
    if (fam == FAM_BB7) {
      lam <- bb7_par2tail(th1, th2)   # (λL, λU)
      tip_or_diffuse <-
        if (is.null(tipm)) {
          dbeta(lam[1], 2, 2, log=TRUE) + dbeta(lam[2], 2, 2, log=TRUE)
        } else {
          dlogitnorm(lam[1], logit(tipm["mL"]), s) +
            dlogitnorm(lam[2], logit(tipm["mU"]), s)
        }
      lp <- lp + tip_or_diffuse + bb7_log_jacobian(lam[1], lam[2]); next
    }
    
    if (fam == FAM_BB7R180) {
      lamr <- bb7r180_par2tail(th1, th2)   # (λL^rot, λU^rot)
      tip_or_diffuse <-
        if (is.null(tipm)) {
          dbeta(lamr[1], 2, 2, log=TRUE) + dbeta(lamr[2], 2, 2, log=TRUE)
        } else {
          dlogitnorm(lamr[1], logit(tipm["mU"]), s) +
            dlogitnorm(lamr[2], logit(tipm["mL"]), s)
        }
      lp <- lp + tip_or_diffuse + bb7r180_log_jacobian(lamr[1], lamr[2]); next
    }
    
    if (fam == FAM_BB8R180) {
      lamLr <- bb8r180_par2tail(th1, th2)[1]  # only λ_L^rot is nonzero
      tip_or_diffuse <-
        if (is.null(tipm)) {
          dbeta(lamLr, 2, 2, log=TRUE)
        } else {
          dlogitnorm(lamLr, logit(tipm["mL"]), s)
        }
      lp <- lp + tip_or_diffuse + bb8r180_log_jacobian_1d(lamLr); next
    }
    
    if (fam == FAM_T) {
    lam <- t_par2tail(th1, th2)
    tip_or_diffuse <-
      if (is.null(tipm)) {
        dbeta(lam, 2, 2, log = TRUE)
      } else {
        m <- max(tipm["mL"], tipm["mU"])
        dlogitnorm(lam, logit(m), s)  # symmetric: use max/mU; see seeding below
      }
    lp_jac <- t_log_jacobian(lam, th2)
    lp_df <- if (th2 >= 2 && th2 <= 30) -log(th2) - log(log(30)) else -Inf
    lp <- tip_or_diffuse + lp_jac + lp_df
    return(fillna_neg(lp))
  }

    
    stop("unknown family in log_prior_with_tip_cached")
  }
  lp
}



## empirical means for edge e at time t (Tree 1 only; NULL otherwise)
.tip_means_for_edge_t <- function(e, data_up_to_t, cfg) {
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
  lp <- 0
  for (e in seq_len(cfg$K)) {
    fam  <- p$fam[e]; th1 <- p$th1[e]; th2 <- p$th2[e]
    tipm <- .tip_means_for_edge_t(e, data_up_to_t, cfg)  # NULL if not Tree 1
    lp   <- lp + log_prior_edge_strong(fam, th1, th2, cfg, tipm)
  }
  lp
}


## Build a cache of empirical tails for all edges at time t (Tree 1 only).
## Returns a list length K; each Tree-1 edge has c(mL=..., mU=...), else NULL.
precompute_tip_means <- function(data_up_to_t, cfg) {
  K <- cfg$K
  out <- vector("list", K)
  if (is.null(cfg$edge_pair)) return(out)
  for (e in seq_len(K)) {
    if (cfg$edge_tree[e] != 1L) {
      next
    } else {
      local_idx <- sum(cfg$edge_tree[seq_len(e)] == 1L)
      pair      <- cfg$edge_pair[local_idx, ]
      uv        <- data_up_to_t[, pair, drop = FALSE]
      emps      <- emp_tails_FRAPO(uv, method = cfg$tip_method, k = cfg$tip_k)
      out[[e]]  <- c(mL = as.numeric(emps["L"]), mU = as.numeric(emps["U"]))
    }
  }
  out
}


## Precompute tail weights for weighted log-likelihood (if used)
precompute_tail_weights <- function(data_up_to_t, cfg) {
  if (isTRUE(cfg$use_weighted_ll)) {
    tail_weights(
      data_up_to_t,
      tau = cfg$tauL,
      k   = cfg$joint_k,
      eps = cfg$tail_eps
    )
  } else {
    NULL
  }
}



######################################

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




diagnostic_report <- function(t, tr, U, particles, w_new,
                              cfg, q_probs = c(0.025, 0.975)) {
  # ==== inputs & basic flags ====
  fam_mat <- particles$fam_mat            # M × K (integer)
  th1_mat <- particles$th1_mat            # M × K (numeric)
  th2_mat <- particles$th2_mat            # M × K (numeric)
  M       <- nrow(fam_mat)
  K       <- ncol(fam_mat)
  
  # pick weights: prefer explicit w_new; fall back to particles$w
  if (missing(w_new) || is.null(w_new) || length(w_new) != M) {
    w_new <- particles$w
  }
  w_new <- w_new / sum(w_new)
  
  k_step     <- cfg$k_step
  print_flag <- (t %% k_step == 0L) || (t == nrow(U))
  
  # ==== family codes (short aliases) ====
  IND  <- FAM_INDEP
  GAU  <- FAM_GAUSS
  BB1c <- FAM_BB1
  BB1s <- FAM_BB1R180
  BB8s <- FAM_BB8R180
  BB7  <- FAM_BB7
  BB7s <- FAM_BB7R180
  TT   <- FAM_T
  
  # ==== weights & ESS ====
  ess_t <- 1 / sum(w_new^2)
  
  # ==== sparsity & diversity ====
  dep_mask   <- fam_mat != IND
  edge_ct    <- rowSums(dep_mask)
  mean_edges <- w_mean(edge_ct, w_new)
  se_edges   <- mc_se(edge_ct, w_new, ess_t)
  sparsity   <- 1 - mean_edges / K
  
  key_vec  <- apply(th1_mat, 1, \(row) paste(row, collapse = ","))
  n_unique <- length(unique(key_vec))
  
  if (M >= 2) {
    dists    <- as.matrix(stats::dist(th1_mat))
    avg_dist <- mean(dists[lower.tri(dists)])
  } else {
    avg_dist <- 0
  }
  
  if (print_flag) {
    cat(sprintf(
      "t=%4d | tr=%2d | ESS/M=%.3f | dep.edges=%.2f ± %.2f | sparsity=%.3f | unique=%d | L2=%.4f\n",
      t, tr, ess_t / M, mean_edges, se_edges, sparsity, n_unique, avg_dist))
  }
  
  # ==== family proportions (by cfg$families order) ====
  fam_codes <- FAM_INFO$code[match(cfg$families, FAM_INFO$name)]
  prop_line <- vapply(seq_along(fam_codes), function(j) {
    code <- fam_codes[j]
    if (is.na(code)) {
      sprintf("%s n/a", cfg$families[j])
    } else {
      p_j <- w_mean(rowSums(fam_mat == code), w_new) / K
      sprintf("%s %.2f", cfg$families[j], p_j)
    }
  }, character(1))
  
  if (print_flag)
    cat("   Family proportions | ", paste(prop_line, collapse = " | "), "\n")
  
  # ==== per-edge summaries ====
  edge_summ <- vector("list", K)
  
  for (e in seq_len(K)) {
    p_indep <- sum(w_new[fam_mat[, e] == IND ])
    p_gauss <- sum(w_new[fam_mat[, e] == GAU ])
    p_bb1   <- sum(w_new[fam_mat[, e] == BB1c])
    p_bb1r  <- sum(w_new[fam_mat[, e] == BB1s])
    p_bb8r  <- sum(w_new[fam_mat[, e] == BB8s])
    p_bb7   <- sum(w_new[fam_mat[, e] == BB7 ])
    p_bb7r  <- sum(w_new[fam_mat[, e] == BB7s])
    p_t     <- sum(w_new[fam_mat[, e] == TT  ])
    
    res <- list(edge   = e,
                p_indep = p_indep,
                p_gauss = p_gauss,
                p_bb1   = p_bb1,
                p_bb1r  = p_bb1r,
                p_bb8r  = p_bb8r)
    # (kept same columns as your original; add p_bb7/p_bb7r/p_t if you want)
    
    if (print_flag && p_indep > 0)
      cat(sprintf("     Edge %2d  Indep   : P=%.3f\n", e, p_indep))
    
    # ---- Gaussian
    if (p_gauss > 0) {
      mask   <- fam_mat[, e] == GAU
      rho_e  <- th1_mat[mask, e]
      w_g    <- w_new[mask] / p_gauss
      mu_rho <- sum(w_g * rho_e)
      sd_rho <- sqrt(sum(w_g * (rho_e - mu_rho)^2))
      qs     <- w_quantile(rho_e, w_g, q_probs)
      
      res <- c(res, list(mu_rho = mu_rho, sd_rho = sd_rho,
                         rho_q025 = qs[1], rho_q975 = qs[2]))
      if (print_flag)
        cat(sprintf("                Gaussian: P=%.3f | ρ=%.3f (SD %.3f) CI=[%.3f,%.3f]\n",
                    p_gauss, mu_rho, sd_rho, qs[1], qs[2]))
    }
    
    # ---- BB1
    if (p_bb1 > 0) {
      mask   <- fam_mat[, e] == BB1c
      th1_e  <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      lam    <- t(mapply(bb1_par2tail, th1_e, th2_e))
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
    
    # ---- BB1^180
    if (p_bb1r > 0) {
      mask    <- fam_mat[, e] == BB1s
      th1_e   <- th1_mat[mask, e]; th2_e <- th2_mat[mask, e]
      lam_rot <- t(mapply(bb1r180_par2tail, th1_e, th2_e))
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
    
    # ---- BB8^180
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
    
    # ---- BB7
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
    
    # ---- BB7^180
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
    
    # ---- Student-t
    if (p_t > 0) {
      mask <- fam_mat[, e] == TT
      rho_e <- th1_mat[mask, e]
      nu_e  <- th2_mat[mask, e]
      wloc  <- w_new[mask] / p_t
      stats <- function(x){
        mu <- sum(wloc * x); sd <- sqrt(sum(wloc * (x - mu)^2)); qs <- w_quantile(x, wloc, q_probs)
        c(mu, sd, qs)
      }
      s_rho <- stats(rho_e); s_nu <- stats(nu_e)
      if (print_flag)
        cat(sprintf("                t      : P=%.3f | rho=%.3f±%.3f CI=[%.3f,%.3f] | nu=%.2f±%.2f CI=[%.2f,%.2f]\n",
                    p_t, s_rho[1], s_rho[2], s_rho[3], s_rho[4], s_nu[1], s_nu[2], s_nu[3], s_nu[4]))
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
  lik <- vapply(seq_len(cfg$M), function(i) {
    vine_i <- .build_vine_from_vectors(
      fam = particles$fam_mat[i, ],
      th1 = particles$th1_mat[i, ],
      th2 = particles$th2_mat[i, ],
      skel    = skel
    )
    as.numeric(dvinecop(u_obs, vine_i))
  }, numeric(1))
  log_pred_density <- log(sum(w_prev_for_prediction * lik) + 1e-30)
  
  list(
    log_pred_density = log_pred_density
  )
}

# Empirical Rosenblatt PIT from predictive copula draws (mixture)
# u: numeric length d (observed copula row at time t)
# U_draws: L x d matrix of mixture draws in copula space
# K: number of neighbors used for conditional CDFs (default ~ sqrt(L))
# Fast empirical Rosenblatt using mixture draws
# u: numeric length d (observed copula row)
# U_draws: L x d numeric matrix (mixture draws)
# K: neighbors (default ~ sqrt(L), min 30)
empirical_rosenblatt_from_draws <- function(u, U_draws, K = NULL) {
  # --- defensive coercions ---
  if (is.matrix(u)) u <- as.numeric(u[1, ])
  u <- as.numeric(u)
  U_draws <- as.matrix(U_draws)
  storage.mode(U_draws) <- "double"
  
  L <- nrow(U_draws); d <- ncol(U_draws)
  stopifnot(length(u) == d, L >= 1, d >= 1)
  
  if (is.null(K)) K <- max(30L, floor(sqrt(L)))
  K <- as.integer(min(K, L))
  
  z <- numeric(d)
  
  # Z1: marginal empirical CDF (with add-1 smoothing)
  c1 <- sum(U_draws[, 1] <= u[1])
  z[1] <- (c1 + 1) / (L + 2)
  
  if (d == 1L) return(z)
  
  # We’ll grow the conditioning space incrementally
  Xprev <- U_draws[, 1, drop = FALSE]
  u_prev <- u[1]
  
  use_fnn <- requireNamespace("FNN", quietly = TRUE)
  
  for (j in 2:d) {
    if (use_fnn) {
      # Fast kNN indices in (j-1)-dim space
      idx <- FNN::get.knnx(Xprev, t(as.matrix(u_prev)), k = K)$nn.index[1, ]
    } else {
      # Fallback: full sort (robust across R versions)
      if (ncol(Xprev) == 1L) {
        dist2 <- abs(Xprev[, 1] - u_prev)
      } else {
        diff  <- sweep(Xprev, 2, u_prev, FUN = "-")
        dist2 <- sqrt(rowSums(diff * diff))
      }
      idx <- head(order(as.numeric(dist2)), K)
    }
    
    # Conditional CDF estimate with add-1 smoothing
    cj <- sum(U_draws[idx, j] <= u[j])
    z[j] <- (cj + 1) / (K + 2)
    
    # Expand conditioning set for next step
    if (j < d) {
      Xprev <- cbind(Xprev, U_draws[, j, drop = FALSE])
      u_prev <- c(u_prev, u[j])
    }
  }
  z
}




# compute_predictive_metrics <- function(u_obs, particles, skel,
#                                        w_prev_for_prediction, cfg,
#                                        q_probs = c(0.025, 0.975)) {
#   
# M <- nrow(particles$fam_mat)
# d <- ncol(skel$structure)
# vines <- lapply(seq_len(M), function(i) {
#   .build_vine_from_vectors(
#     fam = particles$fam_mat[i, ],
#     th1 = particles$th1_mat[i, ],
#     th2 = particles$th2_mat[i, ],
#     skel    = skel
#   )
# })
#   
# 
# lik <- vapply(seq_len(cfg$M), function(i) {
#   as.numeric(dvinecop(u_obs, vines[[i]]))
#   }, numeric(1))
# log_pred_density <- log(sum(w_prev_for_prediction * lik) + 1e-30)
# 
# Ci <- vapply(seq_len(cfg$M), function(i) {
#   as.numeric(pvinecop(u_obs, vines[[i]]))
# }, numeric(1))
# pit_joint <- sum(w_prev_for_prediction * Ci)
# 
# list(
#   log_pred_density = log_pred_density,
#   pit_joint = pit_joint
# )
# }


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







residual_counts <- function(w, L) {
  w <- pmax(w, 0); w <- w/sum(w)
  base <- floor(L * w)
  R <- L - sum(base)
  if (R > 0) {
    r <- L*w - base
    if (sum(r) > 0) {
      r <- r / sum(r)
      base <- base + as.vector(rmultinom(1, size = R, prob = r))
    }
  }
  base
}





smc_predictive_sample_fast2_scoped <- function(particles, skel, w, L = 10000, cl = NULL) {
  # --- normalize weights robustly ---
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0
  sw <- sum(w)
  if (sw <= 0) stop("All weights are zero/invalid.")
  w <- w / sw
  
  # --- sanity & sizes ---
  if (!is.matrix(particles$fam_mat) || !is.matrix(particles$th1_mat) || !is.matrix(particles$th2_mat))
    stop("particles must contain fam_mat, th1_mat, th2_mat matrices.")
  M <- nrow(particles$fam_mat)
  if (length(w) != M) stop("length(w) must equal number of particles (nrow(fam_mat)).")
  if (L < 1L) stop("L must be >= 1.")
  
  # --- counts per particle (sum == L) ---
  cnt  <- residual_counts(w, L)
  used <- which(cnt > 0L)
  if (length(used) == 0L) {
    d0 <- ncol(skel$structure)
    return(matrix(numeric(0), nrow = 0, ncol = d0))
  }
  
  # --- build each used vine ON MASTER (once) ---
  vines_used <- lapply(used, function(i) {
    .build_vine_from_vectors(
      fam = as.integer(particles$fam_mat[i, ]),
      th1 = as.numeric(particles$th1_mat[i, ]),
      th2 = as.numeric(particles$th2_mat[i, ]),
      skel = skel
    )
  })
  
  # --- SERIAL PATH ---
  if (is.null(cl)) {
    d_ref <- NULL
    sims  <- vector("list", length(used))
    for (k in seq_along(used)) {
      ni  <- cnt[used[k]]
      sim <- rvinecop(ni, vines_used[[k]], cores = 1)
      vec <- as.numeric(sim)
      if (length(vec) %% ni != 0L)
        stop(sprintf("rvinecop length mismatch: len=%d, ni=%d", length(vec), ni))
      d_k <- length(vec) / ni
      dim(vec) <- c(ni, d_k)
      if (is.null(d_ref)) d_ref <- d_k else if (d_k != d_ref)
        stop(sprintf("Inconsistent simulated dimension: %d vs %d", d_k, d_ref))
      sims[[k]] <- vec
    }
    out <- do.call(rbind, sims)
    dimnames(out) <- NULL
    return(out)
  }
  
  ## ================= PARALLEL: export-once + batched indices =================
  W <- length(cl)
  parallel::clusterExport(cl, varlist = c("vines_used", "cnt", "used"), envir = environment())
  
  # greedy bin packing by workload
  ord  <- order(cnt[used], decreasing = TRUE)
  idxs <- seq_along(used)[ord]
  bins <- vector("list", W)
  load <- integer(W)
  for (j in idxs) {
    k <- which.min(load)
    bins[[k]] <- c(bins[[k]], j)           # index into 'used' / 'vines_used'
    load[k]   <- load[k] + cnt[ used[j] ]
  }
  
  worker_sim <- function(local_idx) {
    if (length(local_idx) == 0L) return(NULL)
    
    # First job: infer d and init parts
    first <- local_idx[[1]]
    ni1   <- cnt[ used[first] ]
    sim1  <- rvinecop(ni1, vines_used[[first]], cores = 1)
    v1    <- as.numeric(sim1)
    if (length(v1) %% ni1 != 0L)
      stop(sprintf("rvinecop length mismatch on worker: len=%d, ni=%d", length(v1), ni1))
    d <- length(v1) / ni1
    dim(v1) <- c(ni1, d)
    out_parts <- list(v1)
    
    # Remaining jobs
    if (length(local_idx) > 1L) {
      for (t in local_idx[-1]) {
        ni  <- cnt[ used[t] ]
        sim <- rvinecop(ni, vines_used[[t]], cores = 1)
        vv  <- as.numeric(sim)
        if (length(vv) %% ni != 0L)
          stop(sprintf("rvinecop length mismatch on worker: len=%d, ni=%d", length(vv), ni))
        dim(vv) <- c(ni, d)
        out_parts[[length(out_parts) + 1L]] <- vv
      }
    }
    do.call(rbind, out_parts)
  }
  
  sims_list <- parallel::parLapply(cl, bins, worker_sim)
  
  sims_list <- Filter(Negate(is.null), sims_list)
  d_vec <- vapply(sims_list, ncol, integer(1))
  if (length(unique(d_vec)) != 1L)
    stop(sprintf("Inconsistent simulated dimension across workers: %s",
                 paste(unique(d_vec), collapse = ", ")))
  out <- do.call(rbind, sims_list)
  dimnames(out) <- NULL
  out
}






smc_predictive_sample_fast2_scoped2 <- function(particles, skel, w, L = 10000, cl = NULL) {
  # normalize weights
  w[!is.finite(w)] <- 0; w[w < 0] <- 0
  sw <- sum(w); if (sw <= 0) stop("All weights are zero/invalid.")
  w <- w / sw
  
  # sanity
  stopifnot(is.matrix(particles$fam_mat), is.matrix(particles$th1_mat), is.matrix(particles$th2_mat))
  M <- nrow(particles$fam_mat)
  if (length(w) != M) stop("length(w) must equal number of particles (nrow(fam_mat)).")
  if (L < 1L) stop("L must be >= 1.")
  
  # counts per particle
  cnt  <- residual_counts(w, L)
  used <- which(cnt > 0L)
  if (!length(used)) {
    d0 <- ncol(skel$structure)
    return(matrix(numeric(0), nrow = 0, ncol = d0))
  }
  
  # ---------- SERIAL PATH ----------
  if (is.null(cl)) {
    sims <- vector("list", length(used))
    d_ref <- NULL
    for (k in seq_along(used)) {
      i   <- used[k]
      ni  <- cnt[i]
      vine_i <- .build_vine_from_vectors(
        fam  = as.integer(particles$fam_mat[i, ]),
        th1  = as.numeric(particles$th1_mat[i, ]),
        th2  = as.numeric(particles$th2_mat[i, ]),
        skel = skel
      )
      sim <- rvinecop(ni, vine_i, cores = 1)  # already ni × d
      if (!is.matrix(sim)) sim <- as.matrix(sim)
      if (is.null(d_ref)) d_ref <- ncol(sim)
      else if (ncol(sim) != d_ref)
        stop(sprintf("Inconsistent simulated dimension: %d vs %d", ncol(sim), d_ref))
      sims[[k]] <- sim
    }
    out <- do.call(rbind, sims)
    rownames(out) <- NULL
    return(out)
  }
  
  # ---------- PARALLEL PATH (export specs once; pass indices only) ----------
  fam_mat <- particles$fam_mat
  th1_mat <- particles$th1_mat
  th2_mat <- particles$th2_mat
  
  parallel::clusterExport(
    cl,
    varlist = c("fam_mat","th1_mat","th2_mat","skel","used","cnt"),
    envir   = environment()
  )
  
  # greedy bin packing by workload
  W <- length(cl)
  ord  <- order(cnt[used], decreasing = TRUE)
  idxs <- seq_along(used)[ord]
  bins <- vector("list", W)
  load <- integer(W)
  for (j in idxs) {
    k <- which.min(load)
    bins[[k]] <- c(bins[[k]], j)      # index into 'used'
    load[k]   <- load[k] + cnt[ used[j] ]
  }
  
  worker_sim <- function(local_idx) {
    if (!length(local_idx)) return(NULL)
    parts <- vector("list", length(local_idx))
    d_loc <- NULL
    for (m in seq_along(local_idx)) {
      u_pos <- used[ local_idx[m] ]   # actual particle row
      ni    <- cnt[ u_pos ]
      
      vine_i <- .build_vine_from_vectors(
        fam  = as.integer(fam_mat[u_pos, ]),
        th1  = as.numeric(th1_mat[u_pos, ]),
        th2  = as.numeric(th2_mat[u_pos, ]),
        skel = skel
      )
      
      sim <- rvinecop(ni, vine_i, cores = 1)
      if (!is.matrix(sim)) sim <- as.matrix(sim)
      if (is.null(d_loc)) d_loc <- ncol(sim)
      else if (ncol(sim) != d_loc)
        stop(sprintf("Worker dimension mismatch: %d vs %d", ncol(sim), d_loc))
      parts[[m]] <- sim
    }
    do.call(rbind, parts)
  }
  
  sims_list <- parallel::parLapply(cl, bins, worker_sim)
  sims_list <- Filter(Negate(is.null), sims_list)
  
  # final consistency + combine
  d_vec <- vapply(sims_list, ncol, integer(1))
  if (length(unique(d_vec)) != 1L)
    stop(sprintf("Inconsistent simulated dimension across workers: %s",
                 paste(unique(d_vec), collapse = ", ")))
  out <- do.call(rbind, sims_list)
  rownames(out) <- NULL
  out
}


smc_predictive_sample_fast2_grouped_epoch <- function(
    w,                    # numeric length M (you can pass particles$w)
    L = 10000,
    cl,                   # PSOCK cluster (REQUIRED for this epoch version)
    round_digits = 6,     # grouping granularity
    top_k = NULL,         # optionally keep only top K groups by count
    min_count = 0,        # optionally drop tiny groups
    nc = 1                # threads per worker for rvinecop; keep small
) {
  # --- normalize weights robustly ---
  w[!is.finite(w)] <- 0; w[w < 0] <- 0
  sw <- sum(w); if (sw <= 0) stop("All weights are zero/invalid.")
  w <- w / sw
  
  # discover M and d from workers (state must already be pushed)
  Ms <- parallel::clusterEvalQ(cl, nrow(.vine_epoch$fam_mat))
  if (!all(unlist(Ms) == Ms[[1]])) stop("Inconsistent M across workers.")
  M <- Ms[[1]]
  if (length(w) != M) stop("length(w) must equal number of particles.")
  
  # counts per particle
  cnt  <- residual_counts(w, L)
  used <- which(cnt > 0L)
  if (!length(used)) {
    d0 <- parallel::clusterEvalQ(cl, ncol(.vine_epoch$skel$structure))[[1]]
    return(matrix(numeric(0), nrow = 0, ncol = d0))
  }
  
  # compute grouping keys ON A WORKER (no big exports)
  parallel::clusterExport(cl, c("used","round_digits"), envir = environment())
  keys_list <- parallel::clusterEvalQ(cl, {
    fm <- .vine_epoch$fam_mat; t1 <- .vine_epoch$th1_mat; t2 <- .vine_epoch$th2_mat
    kfun <- function(i) paste(c(fm[i,], round(t1[i,], round_digits), round(t2[i,], round_digits)), collapse="|")
    vapply(used, kfun, "", USE.NAMES = FALSE)
  })
  keys <- keys_list[[1]]                 # they’re identical on each worker
  
  ukeys <- unique(keys)
  gix   <- match(keys, ukeys)
  rep_row     <- as.integer(tapply(used, gix, function(ix) ix[1]))
  group_count <- as.numeric(tapply(cnt[used], gix, sum))
  
  # optional pruning
  if (!is.null(top_k)) {
    ord <- order(group_count, decreasing = TRUE)
    keep <- ord[seq_len(min(top_k, length(ord)))]
    rep_row     <- rep_row[keep]
    group_count <- group_count[keep]
  }
  if (min_count > 0) {
    keep <- which(group_count >= min_count)
    if (!length(keep)) {
      d0 <- parallel::clusterEvalQ(cl, ncol(.vine_epoch$skel$structure))[[1]]
      return(matrix(numeric(0), nrow = 0, ncol = d0))
    }
    rep_row     <- rep_row[keep]
    group_count <- group_count[keep]
  }
  
  # bin packing (balanced)
  W <- length(cl)
  ord <- order(group_count, decreasing = TRUE)
  idxs <- seq_along(rep_row)[ord]
  bins <- vector("list", W)
  load <- numeric(W)
  for (j in idxs) {
    k <- which.min(load)
    bins[[k]] <- c(bins[[k]], j)
    load[k]   <- load[k] + group_count[j]
  }
  
  # worker function: uses cached vines via .worker_get_vine()
  worker_sim_epoch <- function(local_idx, rep_row, group_count, nc = 1) {
    if (!length(local_idx)) return(NULL)
    parts <- vector("list", length(local_idx))
    d_loc <- NULL
    for (m in seq_along(local_idx)) {
      j  <- local_idx[m]
      i  <- rep_row[j]
      ni <- group_count[j]
      vine_i <- .worker_get_vine(i)                 # <<— cached
      sim <- rvinecop(ni, vine_i, cores = nc)
      sim <- as.matrix(sim)
      if (is.null(d_loc)) d_loc <- ncol(sim)
      else if (ncol(sim) != d_loc) stop("Worker dimension mismatch")
      parts[[m]] <- sim
    }
    do.call(rbind, parts)
  }
  
  parallel::clusterExport(cl, c("worker_sim_epoch","rep_row","group_count","bins","nc"), envir = environment())
  sims_list <- parallel::parLapplyLB(cl, bins, worker_sim_epoch, rep_row, group_count, nc)
  sims_list <- Filter(Negate(is.null), sims_list)
  d_vec <- vapply(sims_list, ncol, integer(1))
  if (length(unique(d_vec)) != 1L) stop("Inconsistent simulated dimension across workers")
  out <- do.call(rbind, sims_list)
  rownames(out) <- NULL
  out
}



# Push the big read-only state to workers ONCE per time step (epoch).
push_epoch_state <- function(cl, particles, skel, epoch_id = NULL) {
  stopifnot(!is.null(cl))
  fam_mat <- particles$fam_mat
  th1_mat <- particles$th1_mat
  th2_mat <- particles$th2_mat
  
  parallel::clusterExport(cl, c("fam_mat","th1_mat","th2_mat","skel","epoch_id"), envir = environment())
  
  parallel::clusterEvalQ(cl, {
    if (!exists(".vine_epoch", inherits = FALSE)) .vine_epoch <- new.env(parent = emptyenv())
    .vine_epoch$fam_mat  <- fam_mat
    .vine_epoch$th1_mat  <- th1_mat
    .vine_epoch$th2_mat  <- th2_mat
    .vine_epoch$skel     <- skel
    .vine_epoch$epoch_id <- epoch_id
    
    # Simple per-epoch vine cache
    if (!exists(".vine_cache", inherits = FALSE)) .vine_cache <- new.env(parent = emptyenv())
    .vine_cache$epoch_id <- epoch_id
    .vine_cache$vines    <- new.env(parent = emptyenv())
    
    rm(fam_mat, th1_mat, th2_mat, skel, epoch_id)
    invisible(TRUE)
  })
}

# Optional tidy-up when you are completely done (not needed each step)
clear_epoch_state <- function(cl) {
  parallel::clusterEvalQ(cl, { 
    if (exists(".vine_epoch", inherits = FALSE)) rm(.vine_epoch, inherits = FALSE)
    if (exists(".vine_cache", inherits = FALSE)) rm(.vine_cache, inherits = FALSE)
    invisible(TRUE)
  })
}

# Worker helper: build or fetch cached vine for particle row i
.worker_get_vine <- function(i) {
  env   <- .vine_epoch
  cache <- .vine_cache
  if (!identical(cache$epoch_id, env$epoch_id)) {  # epoch changed → reset cache
    cache$vines    <- new.env(parent = emptyenv())
    cache$epoch_id <- env$epoch_id
  }
  key <- as.character(i)
  v <- get0(key, envir = cache$vines, inherits = FALSE)
  if (!is.null(v)) return(v)
  
  v <- .build_vine_from_vectors(
    fam  = as.integer(env$fam_mat[i, ]),
    th1  = as.numeric(env$th1_mat[i, ]),
    th2  = as.numeric(env$th2_mat[i, ]),
    skel = env$skel
  )
  assign(key, v, envir = cache$vines)
  v
}

# Register the worker helper once after you create the cluster
register_worker_helpers <- function(cl) {
  parallel::clusterExport(cl, c(".worker_get_vine"), envir = environment())
}


