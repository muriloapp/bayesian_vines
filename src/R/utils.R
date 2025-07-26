


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
                         eps   = 1e-6,   # lower-bound cushion
                         upper = 7 - 1e-6) {
  
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
  fn <- file.path(dir_out,
                  sprintf("%s_%s.rds", res$alg, res$cfg$label))
  saveRDS(res, fn)
  message("✓ saved → ", fn)
}


make_skeleton_CVM <- function(U_train) {
  CVM <- RVineStructureSelect(U_train, familyset = c(1), type=1, indeptest = TRUE, level = 0.1)
  old_M  <- CVM$Matrix          
  order  <- old_M[, 1]
  skeleton <- vinecop(U_train, family_set = "gaussian", structure = cvine_structure(order))
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

