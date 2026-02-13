

port_stats <- function(r_vec, alphas) {

  # downside risk
  losses <- r_vec                  # “loss = −return”
  VaR <- vapply(alphas, function(a)  quantile(losses, a), numeric(1))
  ES  <- vapply(seq_along(alphas), function(k) {
    thr <- VaR[k]
    mean(losses[losses <= thr])
  }, numeric(1))
  
  list(VaR = VaR, ES = ES)
}


risk_stats_full <- function(R_t,           # L × d matrix  (future PIT‑draws ⇒ returns)
                            #mu_fc,            # length‑d numeric
                            #sig_fc,           # length‑d numeric
                            alphas)
{

  
  ## 2.  VaR & ES  -------------------------------------------------------------
  losses <- R_t
  VaR <- sapply(alphas, function(a)
    apply(losses, 2, quantile, probs = a))  
  ES  <- sapply(seq_along(alphas), function(k) {       # loop over α’s
    vapply(seq_len(ncol(losses)), function(j) {  # loop over assets
      thr <- VaR[j, k]                           # ← correct index
      mean(losses[losses[, j] <= thr, j])        # tail average
    }, numeric(1))
  })
  
  ## 3.  Return as a clean list  ----------------------------------------------
  list(
    VaR  = VaR,                          # d × length(alphas)
    ES   = ES                            # d × length(alphas)
  )
}

# ------------ PINBALL (Quantile) LOSS ------------
# y: scalar realized value; q: scalar forecast quantile at level alpha
pinball_loss <- function(y, q, alpha) {
  u <- y - q
  (alpha - (u < 0)) * u
}

# Vectorized over assets & alphas:
# y_vec: length-d realized vector; VaR: d x A matrix of VaR forecasts
pinball_matrix <- function(y_vec, VaR, alphas) {
  d <- length(y_vec); A <- length(alphas)
  out <- matrix(NA_real_, d, A)
  for (j in seq_len(d)) {
    for (k in seq_len(A)) 
      out[j, k] <- pinball_loss(y_vec[1,j], VaR[j, k], alphas[k])
  }
  out
}

# ------------ FISSLER–ZIEGEL JOINT LOSS (VaR, ES) ------------


fzl_pzc_matrix <- function(y_vec, VaR, ES, alphas) {
  d <- length(y_vec); A <- length(alphas)
  if (!all(dim(VaR) == c(d, A)) || !all(dim(ES) == c(d, A)))
    stop("Dimensions of VaR/ES must be d x A.")
  if (any(ES >= 0)) stop("ES entries must be < 0 for this score.")
  out <- matrix(NA_real_, d, A)
  for (k in seq_len(A)) {
    I <- as.numeric(y_vec <= VaR[, k])
    out[, k] <- as.numeric((I * (y_vec - VaR[, k])) / (alphas[k] * ES[, k]) +
                             (VaR[, k] / ES[, k]) + log(-ES[, k]) - 1)
  }
  out
}

fzl_pzc_scalar <- function(y, v, e, alpha) {
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0,1).")
  if (!is.finite(y) || !is.finite(v) || !is.finite(e)) stop("Non-finite input.")
  if (e >= 0) stop("This FZL form assumes ES < 0 (left tail).")
  I <- as.numeric(y <= v)
  (I * (y - v)) / (alpha * e) + (v / e) + log(-e) - 1
}




wcrps_gr <- function(draws, y, ngrid = 1000) {
  # draws: numeric vector of predictive draws
  # y: realized value (scalar)
  # ngrid: number of grid points
  rng <- range(c(draws, y))
  grid <- seq(rng[1], rng[2], length.out = ngrid)
  
  # empirical CDF at each grid point
  Fhat <- ecdf(draws)
  Fz   <- Fhat(grid)
  
  # weight function omega(z) = 1 - Phi(z)
  w <- 1 - pnorm(grid)
  
  # indicator 1{y <= z}
  Iy <- as.numeric(y <= grid)
  
  integrand <- w * (Fz - Iy)^2
  
  # numeric integration by trapezoid
  trapz <- sum(0.5 * (integrand[-1] + integrand[-length(integrand)]) *
                 diff(grid))
  
  trapz
}

# vectorized across d assets:
wcrps_gr_matrix <- function(draws, y_vec, ngrid = 1000) {
  d <- ncol(draws); out <- numeric(d)
  for (j in seq_len(d)) {
    out[j] <- wcrps_gr(draws[, j], y_vec[1,j], ngrid = ngrid)
  }
  out
}

# Weighted CRPS (left-tail) for a single forecast distribution and realization
# draws: numeric vector of portfolio predictive draws
# y: scalar realized portfolio return (or loss)
# ngrid: grid size for numerical integration
wcrps_gr_scalar <- function(draws, y, ngrid = 2000) {
  rng <- range(c(draws, y))
  grid <- seq(rng[1], rng[2], length.out = ngrid)
  
  # empirical CDF at grid points
  Fhat <- ecdf(draws)
  Fz   <- Fhat(grid)
  
  # weight function: omega(z) = 1 - Phi(z)
  w <- 1 - pnorm(grid)
  
  Iy <- as.numeric(y <= grid)
  integrand <- w * (Fz - Iy)^2
  
  # trapezoidal integration
  sum(0.5 * (integrand[-1] + integrand[-length(integrand)]) * diff(grid))
}




# Tail-CoVaR:  VaR_alpha of portfolio conditional on stock j being in its left tail (<= VaR_j^alpha)
# R_draws: L x d matrix (predictive draws of asset returns at time t)
# r_p:     length-L vector (EW portfolio draws at time t)
# VaRj:    length-d vector of each asset's VaR_j^alpha at time t (same sign convention as R_draws)
# alpha:   quantile for the portfolio VaR (usually same as 0.05 used for the conditioning)
# minN:    minimum number of tail draws to accept before falling back to nearest-neighbor band
# covar_tail_vec <- function(R_draws, r_p, VaRj, alpha = 0.05, minN = 200) {
#   L <- nrow(R_draws); d <- ncol(R_draws)
#   out <- numeric(d)
#   for (j in seq_len(d)) {
#     idx <- which(R_draws[, j] <= VaRj[j])
#     if (length(idx) < minN) {
#       # fallback: take nearest neighbors around VaR_j if the tail set is too small
#       ord <- order(abs(R_draws[, j] - VaRj[j]))
#       idx <- ord[1:minN]
#     }
#     out[j] <- as.numeric(quantile(r_p[idx], probs = alpha, names = FALSE))
#   }
#   out
# }


# Tail-CoVaR for ALL stocks at once:
# R_draws: L×d draws; r_p: length-L portfolio draws; VaRj: length-d vector (at cond_alpha)
# port_alpha: portfolio VaR level; cond_alpha used to build VaRj upstream
covar_tail_vec <- function(R_draws, r_p, VaRj, port_alpha = 0.05, minN = 200) {
  d <- ncol(R_draws)
  out <- numeric(d)
  for (j in seq_len(d)) {
    idx <- which(R_draws[, j] <= VaRj[j])
    if (length(idx) < minN) {
      ord <- order(abs(R_draws[, j] - VaRj[j])); idx <- ord[1:minN]
    }
    out[j] <- as.numeric(quantile(r_p[idx], probs = port_alpha, names = FALSE))
  }
  out
}



covar_tail_vec_asset <- function(R_draws, r_p, VaRj, port_alpha = 0.05, minN = 200) {
  
  d <- ncol(R_draws)
  out <- numeric(d)
  for (j in seq_len(d)) {
    xx  <- 3L - j ### ONLY WORK FOR 2 assets
    idx <- which(R_draws[, j] <= VaRj[j])
    out[j] <- as.numeric(quantile(R_draws[idx, xx], probs = port_alpha, names = FALSE))
  }
  out
}

coes_tail_vec_asset <- function(R_draws, r_p, VaRj, port_alpha = 0.05, minN = 200) {
  
  # ---- avoid anomalies: returns cannot be less than -100% ----------------------
  R_draws <- pmax(R_draws, -1)
  
  d <- ncol(R_draws)
  out <- rep(NA_real_, d)
  
  for (j in seq_len(d)) {
    xx  <- 3L - j  # ONLY WORK FOR 2 assets (other asset index)
    
    idx <- which(R_draws[, j] <= VaRj[j])
    
    # optional: enforce minimum number of distress observations
    if (length(idx) < minN) {
      out[j] <- NA_real_
      next
    }
    
    x_other <- R_draws[idx, xx]
    
    q <- as.numeric(stats::quantile(x_other, probs = port_alpha, names = FALSE, na.rm = TRUE))
    
    out[j] <- mean(x_other[x_other <= q], na.rm = TRUE)
  }
  
  out
}



