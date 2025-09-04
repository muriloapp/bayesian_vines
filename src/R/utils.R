
import_data <- function(path = "data", 
                        n_days = NULL, 
                        n_assets = NULL, 
                        var = NULL, 
                        drop_first_col = FALSE) {
  
  # subset one object; return matrix if drop_first_col=TRUE, else data.table
  subset_data <- function(obj, n_days, n_assets, drop_first_col) {
    if (inherits(obj, "data.table")) {
      dt <- obj
    } else if (is.matrix(obj)) {
      dt <- data.table::as.data.table(obj)
    } else {
      dt <- data.table::as.data.table(obj)
    }
    cn <- colnames(dt)
    if (length(cn) == 0L) return(if (drop_first_col) as.matrix(dt) else dt)
    
    # detect if first column is a date column (by name OR class)
    first_is_date <- identical(cn[1L], "Date") || inherits(dt[[1L]], "Date")
    
    # row subset
    if (!is.null(n_days)) {
      dt <- dt[seq_len(min(n_days, nrow(dt)))]
    }
    
    # figure asset names (exclude first if it's a date column)
    asset_names <- if (first_is_date) cn[-1L] else cn
    
    # limit to n_assets if requested
    if (!is.null(n_assets)) {
      asset_names <- asset_names[seq_len(min(n_assets, length(asset_names)))]
    }
    
    # decide kept columns
    if (drop_first_col) {
      # drop the first date-like column if present; keep only assets
      keep_names <- asset_names
    } else {
      # keep the date column (if present) + assets
      keep_names <- if (first_is_date) c(cn[1L], asset_names) else asset_names
    }
    
    # column subset (data.table-friendly)
    dt <- dt[, ..keep_names]
    
    # return type to match your requirement:
    if (drop_first_col) {
      return(as.matrix(dt))   # pure numeric matrix (no date column)
    } else {
      return(dt)              # data.table with the date + assets
    }
  }
  
  datasets <- list(
    U        = readRDS(file.path(path, "PIT.rds")),
    mu_fc    = readRDS(file.path(path, "returns_mean_forecast.rds")),
    sig_fc   = readRDS(file.path(path, "returns_vol_forecast.rds")),
    df_fc    = readRDS(file.path(path, "df_fc.rds")),
    shape_fc = readRDS(file.path(path, "shape_fc.rds")),
    y_real   = readRDS(file.path(path, "returns_actual.rds"))
  )
  
  if (!is.null(var)) {
    if (!var %in% names(datasets)) {
      stop("Invalid var. Must be one of: ", paste(names(datasets), collapse = ", "))
    }
    return(subset_data(datasets[[var]], n_days, n_assets, drop_first_col))
  }
  
  lapply(datasets, subset_data, n_days = n_days, n_assets = n_assets, drop_first_col = drop_first_col)
}





## RW step

rtnorm1 <- function(mu, sd, a = -0.99, b = 0.99) {
  repeat {
    x <- rnorm(1, mu, sd)
    if (x > a && x < b) return(x)
  }
}

rtnorm_vec <- function(mu_vec, sd, a = -0.99, b = 0.99) {
  x <- mu_vec + sd * rnorm(length(mu_vec))
  pmin(pmax(x, a), b)                     # hard clip (reflection is slower)
}


## BB1 family utils

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

bb1_log_jacobian <- function(lambdaL, lambdaU) {
  log2 <- log(2)
  denom <- (log(1 / lambdaL))^2 * lambdaL * log(2 - lambdaU) * (2 - lambdaU)
  log(log2) - log(denom)
}

sanitize_bb1 <- function(theta, delta,
                         eps   = 0.01,   # lower-bound cushion
                         upper = 7 - 0.01) {
  
  theta <- pmin(pmax(theta, eps),   upper)   # (0 , 7]
  delta <- pmin(pmax(delta, 1+eps), upper)   # (1 , 7]
  c(theta, delta)
}

bb1r180_tail2par <- function(lambdaL_rot, lambdaU_rot) {
  # survival: (λL_rot, λU_rot) correspond to (λU_orig, λL_orig)
  bb1_tail2par(lambdaL = lambdaU_rot, lambdaU = lambdaL_rot)
}

bb1r180_par2tail <- function(theta, delta) {
  lam_orig <- bb1_par2tail(theta, delta)
  # swap for rotated copula’s own tails
  c(lambdaL = lam_orig[2], lambdaU = lam_orig[1])
}

bb1r180_log_jacobian <- function(lambdaL_rot, lambdaU_rot) {
  # J_rot = J_orig evaluated at swapped lambdas
  bb1_log_jacobian(lambdaL = lambdaU_rot, lambdaU = lambdaL_rot)
}

## ---------- BB7 (Joe–Clayton): both tails nonzero ----------
## Tail formulas:
##   λ_L = 2^{-1/δ},         δ > 0
##   λ_U = 2 - 2^{1/θ},      θ ≥ 1
bb7_tail2par <- function(lambdaL, lambdaU) {
  c_log2 <- log(2)
  theta  <- c_log2 / log(2 - lambdaU)      # ≥ 1 if lambdaU ∈ (0,1)
  delta  <- c_log2 / (-log(lambdaL))       # > 0
  c(theta, delta)
}
bb7_par2tail <- function(theta, delta) {
  lambdaL <- 2^(-1 / delta)
  lambdaU <- 2 - 2^(1 / theta)
  c(lambdaL, lambdaU)
}
bb7_log_jacobian <- function(lambdaL, lambdaU) {
  ## J = diag(dδ/dλ_L, dθ/dλ_U); log|J| = log(dδ/dλ_L) + log(dθ/dλ_U)
  c_log2 <- log(2)
  d_delta <- c_log2 / (lambdaL * (log(lambdaL))^2)
  d_theta <- c_log2 / ((log(2 - lambdaU))^2 * (2 - lambdaU))
  log(abs(d_delta)) + log(abs(d_theta))
}

## ---------- BB7^180 (survival rotation): tails swap ----------
## Rotated tails: (λ_L^rot, λ_U^rot) = (λ_U^orig, λ_L^orig)
bb7r180_tail2par <- function(lambdaL_rot, lambdaU_rot) {
  ## map back via the original:
  ##   λ_U^orig = λ_L^rot,  λ_L^orig = λ_U^rot
  bb7_tail2par(lambdaL = lambdaU_rot, lambdaU = lambdaL_rot)
}
bb7r180_par2tail <- function(theta, delta) {
  lam <- bb7_par2tail(theta, delta)
  c(lambdaL_rot = lam[2], lambdaU_rot = lam[1])  # swap
}
bb7r180_log_jacobian <- function(lambdaL_rot, lambdaU_rot) {
  ## same Jacobian as BB7 but with swapped arguments
  bb7_log_jacobian(lambdaL = lambdaU_rot, lambdaU = lambdaL_rot)
}

## ---------- BB8^180 (Joe–Frank survival): single nonzero tail ----------
## Rotated lower tail equals Joe's upper tail; upper tail is 0.
## Use a tail prior on λ_L^rot only (δ from Frank is tail-free).
bb8r180_tail2par <- function(lambdaL_rot, delta_free) {
  c_log2 <- log(2)
  theta  <- c_log2 / log(2 - lambdaL_rot) # from Joe's λ_U formula
  c(theta, delta_free)
}
bb8r180_par2tail <- function(theta, delta) {
  lambdaL_rot <- 2 - 2^(1 / theta)
  c(lambdaL_rot = lambdaL_rot, lambdaU_rot = 0)
}
bb8r180_log_jacobian_1d <- function(lambdaL_rot) {
  c_log2 <- log(2)
  d_theta <- c_log2 / ((log(2 - lambdaL_rot))^2 * (2 - lambdaL_rot))
  log(abs(d_theta))
}

## ---------- Simple sanitizers (mirror your BB1 style) ----------
sanitize_bb7 <- function(theta, delta, eps = 0.01, upper_theta = 6 - 0.01, upper_delta= 25 - 0.1) {
  theta <- pmin(pmax(theta, 1 + eps), upper_theta)   # Joe typically ≥ 1
  delta <- pmin(pmax(delta, 0 + eps), upper_delta)   # keep > 1 like BB1
  c(theta, delta)
}
sanitize_bb8 <- function(theta, delta, eps = 0.01, upper = 7 - 0.01) {
  theta <- pmin(pmax(theta, 1 + eps), upper)   # Joe part
  delta <- pmin(pmax(delta, eps),     upper)   # Frank: > 0
  c(theta, delta)
}

## ---- Student-t tail helpers ----------------------------------------------

t_par2tail <- function(rho, nu) {
  # symmetric: λL = λU = λ
  nu1 <- nu + 1
  s   <- sqrt( nu1 * (1 - rho) / (1 + rho) )
  # λ = 2 * CDF_t(nu+1)(-s)
  2 * stats::pt(-s, df = nu1)
}

t_tail2rho <- function(lambda, nu) {
  lam <- clamp01(lambda, 1e-6, 1 - 1e-6)
  nu1 <- nu + 1
  x   <- stats::qt(lam / 2, df = nu1)               # x <= 0
  s2  <- x^2
  ((nu1) - s2) / ((nu1) + s2)                       # in (-1,1)
}

t_log_jacobian <- function(lambda, nu) {
  lam <- clamp01(lambda, 1e-6, 1 - 1e-6)
  nu1 <- nu + 1
  x   <- stats::qt(lam / 2, df = nu1)               # inverse-CDF
  ft  <- stats::dt(x, df = nu1)                     # pdf at x
  # |dr/dλ| = [ 2 (nu+1) |x| ] / [ ((nu+1)+x^2)^2 * f_t(x) ]
  log( 2 * nu1 * abs(x) ) - 2*log(nu1 + x^2) - log(ft)
}

sanitize_t <- function(rho, nu, rho_max = 0.99, nu_lo = 2, nu_hi = 30) {
  rho <- pmin(pmax(rho, -rho_max), rho_max)
  nu  <- pmin(pmax(nu, nu_lo), nu_hi)
  c(rho, nu)
}



fam_spec <- function(code) {
  row <- FAM_INFO[FAM_INFO$code == code, ]
  list(name = row$name, rv_name = row$rv_name, rotation = row$rotation, npar = row$npar)
}


# Monte-Carlo 

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

## Resample

ESS <- function(w) 1 / sum(w^2)

need_resample <- function(w, cfg, t, T) {
  ess(w) < cfg$ess_thr * cfg$M && t < T
}

systematic_resample <- function(w) {
  M <- length(w)
  u0 <- runif(1)/M
  cumw <- cumsum(w)
  (u0 + 0:(M-1)/M) |> 
    (\(u) findInterval(u, cumw)+1L)()
}

stratified_resample <- function(w) {
  M <- length(w)
  cumw <- cumsum(w)                     # cumulative weights
  u    <- ((seq_len(M) - 1) + runif(M)) / M   # one uniform draw per stratum
  findInterval(u, cumw) + 1L            # ancestor indices (1-based)
}

## Save

save_result <- function(res, dir_out) {
  
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  
  ## Safe fall-backs -----------------------------------------------------
  alg   <- if (length(res$alg)      && nzchar(res$alg[1]))      res$alg[1]      else "unnamedAlg"
  label <- if (!is.null(res$cfg$label) && nzchar(res$cfg$label[1]))
    res$cfg$label[1]
  else format(Sys.time(), "%Y%m%d_%H%M%S")
  
  fn <- file.path(dir_out, sprintf("%s_%s.rds", alg, label))
  
  saveRDS(res, fn)
  message("✓ saved → ", fn)
}


# make_skeleton_CVM <- function(U_train) {
#   CVM <- RVineStructureSelect(U_train, familyset = c(1), type=1, indeptest = TRUE, level = 0.1)
#   old_M  <- CVM$Matrix          
#   order  <- old_M[, 1]
#   skeleton <- vinecop(U_train, family_set = c("gaussian"), structure = cvine_structure(order))
#   skeleton
# }

make_skeleton_CVM <- function(U_train, trunc_tree = 2) {            # NEW arg
  CVM <- RVineStructureSelect(
    U_train, familyset = c(1), type = 1,
    indeptest = TRUE, level = 0.1,
    trunclevel = trunc_tree                                   # NEW
  )
  old_M  <- CVM$Matrix
  order  <- old_M[, 1]
  skeleton <- vinecop(
    U_train,
    family_set = c("gaussian"),
    structure = cvine_structure(order),
    trunc_lvl = trunc_tree                                 # NEW
  )
  skeleton
}


make_cluster <- function(n_cores, seed, exports) {
  cl <- makeCluster(n_cores)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, exports , envir = parent.frame())
  cl
}

active_fams <- function(cfg) FAM_INFO[FAM_INFO$name %in% cfg$families, ]

# edge_tree_map <- function(d) {
#   K   <- d * (d - 1) / 2
#   map <- integer(K)
#   idx <- 1L
#   for (tr in 1:(d - 1)) {          # tree index
#     n_edges <- d - tr              # edges in that tree
#     map[idx:(idx + n_edges - 1)] <- tr
#     idx <- idx + n_edges
#   }
#   map
# }

# edge_tree_map <- function(d, trunc_tree = Inf) {                    # NEW arg
#   K   <- d * (d - 1) / 2
#   map <- integer(K)
#   idx <- 1L
#   max_tr <- min(trunc_tree, d - 1)                                  # NEW
#   for (tr in 1:max_tr) {                                            # CHANGED
#     n_edges <- d - tr
#     map[idx:(idx + n_edges - 1)] <- tr
#     idx <- idx + n_edges
#   }
#   map[map > 0L]                                                     # drop unused tail
# }

K_from_trunc <- function(d, trunc_tree) {
  trunc_tree <- max(1L, min(as.integer(trunc_tree), d - 1L))
  as.integer(d * trunc_tree - trunc_tree * (trunc_tree + 1L) / 2L)
}

edge_tree_map <- function(d, trunc_tree = d - 1L) {
  Tmax <- max(1L, min(as.integer(trunc_tree), d - 1L))
  levels <- integer(0)
  for (t in 1:Tmax) {
    levels <- c(levels, rep.int(t, d - t))
  }
  levels  # length == K_from_trunc(d, Tmax)
}


add_first_tree_map <- function(cfg, skeleton) {
  ord <- skeleton$structure$order      # e.g. c(1, 3, 2)
  d   <- length(ord)
  
  root   <- ord[d]                     # last element = root
  leaves <- ord[seq_len(d - 1)]        # everything before it
  
  # one row per edge, same order as pair_copulas[[1]]
  cfg$edge_pair <- cbind(u = rep(root, d - 1),
                         v = leaves)
  
  cfg
}





# st_inv <- function(U_dt, shape, df_row_dt) {
#   U  <- as.matrix(U_dt)            # L × d numeric matrix
#   ν  <- as.numeric(df_row_dt)      # length-d numeric vector
#   shape  <- as.numeric(shape)
#   sweep(U, 2, shape, ν,  function(u, shape,  df) sn::qst(u, shape, df))
# }

# st_inv <- function(U_dt, shape, df_row_dt) {
#   U      <- as.matrix(U_dt)            # L × d
#   ν      <- as.numeric(df_row_dt)      # length-d
#   shape  <- as.numeric(shape)          # length-d
#   
#   stopifnot(ncol(U) == length(shape), ncol(U) == length(ν))
#   
#   # Apply sn::qst(u, shape, df) to each column
#   Z <- mapply(function(u, s, df) sn::qst(u, xi = 0, omega = 1, alpha = s, nu = df),
#               as.data.frame(U), shape, ν, SIMPLIFY = TRUE)
#   return(Z)
# }

# st_inv <- function(U_dt, shape, df_row_dt) {
#   U     <- as.matrix(U_dt)
#   shape <- as.numeric(shape)
#   nu    <- as.numeric(df_row_dt)
#   
#   out <- matrix(NA_real_, nrow = nrow(U), ncol = ncol(U))
#   for (j in seq_along(shape)) {
#     out[, j] <- sn::qst(U[, j], shape = shape[j], df = nu[j])
#   }
#   colnames(out) <- colnames(U_dt)
#   out
# }

# st_inv_fast <- function(U_dt, shape, df_row_dt) {
#   U          <- as.matrix(U_dt)
#   n_obs      <- nrow(U)
#   n_dim      <- ncol(U)
#   
#   # vectorise everything ----------------------------------------------------
#   p_long     <- as.vector(U)                       # length = n_obs * n_dim
#   shape_long <- rep(as.numeric(shape   ), each = n_obs)
#   df_long    <- rep(as.numeric(df_row_dt), each = n_obs)
#   
#   # single call to qst ------------------------------------------------------
#   out <- sn::qst(p_long, shape = shape_long, df = df_long)
#   
#   # reshape back to L × d ---------------------------------------------------
#   dim(out) <- c(n_obs, n_dim)
#   dimnames(out) <- list(NULL, colnames(U_dt))
#   out
# }


st_inv_fast <- function(U_dt, xi, nu, eps = 1e-12) {
  stopifnot(is.matrix(U_dt))
  U <- pmin(pmax(as.matrix(U_dt), eps), 1 - eps)  # guard 0/1
  n_obs <- nrow(U); n_dim <- ncol(U)
  stopifnot(length(xi) == n_dim, length(nu) == n_dim)
  
  # vectorized single call
  p_long   <- as.vector(U)
  xi_long  <- rep(as.numeric(xi), each = n_obs)
  nu_long  <- rep(as.numeric(nu), each = n_obs)
  out <- fGarch::qsstd(p_long,
                       mean = 0, sd = 1,
                       nu   = nu_long,
                       xi   = xi_long)
  dim(out) <- c(n_obs, n_dim)
  dimnames(out) <- list(NULL, colnames(U_dt))
  out
}



t_inv <- function(U_dt, df_row_dt) {
  U  <- as.matrix(U_dt)            # L × d numeric matrix
  ν  <- as.numeric(df_row_dt)      # length-d numeric vector
  sweep(U, 2, ν, function(u, df) stats::qt(u, df))
}

fillna_neg <- function(x, neg_val = -1e10) {
  
  # handle vectors
  if (is.atomic(x) && is.numeric(x)) {
    x[is.na(x) | is.nan(x)] <- neg_val
    return(x)
  }
  
  # handle matrices / data frames column-wise
  if (is.matrix(x) || is.data.frame(x)) {
    for (j in seq_len(ncol(x))) {
      if (is.numeric(x[[j]]))
        x[[j]][is.na(x[[j]]) | is.nan(x[[j]])] <- neg_val
    }
    return(x)
  }
  
  stop("`x` must be numeric, a numeric matrix, or a data.frame with numeric columns.")
}
