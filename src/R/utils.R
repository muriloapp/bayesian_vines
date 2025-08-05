


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


make_skeleton_CVM <- function(U_train) {
  CVM <- RVineStructureSelect(U_train, familyset = c(1), type=1, indeptest = TRUE, level = 0.1)
  old_M  <- CVM$Matrix          
  order  <- old_M[, 1]
  skeleton <- vinecop(U_train, family_set = c("gaussian"), structure = cvine_structure(order))
  skeleton
}

make_cluster <- function(n_cores, seed, exports) {
  cl <- makeCluster(n_cores)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, exports , envir = parent.frame())
  cl
}

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

st_inv <- function(U_dt, shape, df_row_dt) {
  U     <- as.matrix(U_dt)
  shape <- as.numeric(shape)
  nu    <- as.numeric(df_row_dt)
  
  out <- matrix(NA_real_, nrow = nrow(U), ncol = ncol(U))
  for (j in seq_along(shape)) {
    out[, j] <- sn::qst(U[, j], shape = shape[j], df = nu[j])
  }
  colnames(out) <- colnames(U_dt)
  out
}

st_inv_fast <- function(U_dt, shape, df_row_dt) {
  U          <- as.matrix(U_dt)
  n_obs      <- nrow(U)
  n_dim      <- ncol(U)
  
  # vectorise everything ----------------------------------------------------
  p_long     <- as.vector(U)                       # length = n_obs * n_dim
  shape_long <- rep(as.numeric(shape   ), each = n_obs)
  df_long    <- rep(as.numeric(df_row_dt), each = n_obs)
  
  # single call to qst ------------------------------------------------------
  out <- sn::qst(p_long, shape = shape_long, df = df_long)
  
  # reshape back to L × d ---------------------------------------------------
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
