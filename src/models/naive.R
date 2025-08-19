library(rvinecopulib)
library(VineCopula)
library(here)

source(here("src/R", "utils.R"))
source(here("src/R", "metrics.R"))


build_cfg <- function(d,
                      K = d * (d - 1) / 2) {
  list(
    d       = d,
    K       = K,
    W_predict    = 756L,                                 
    alphas     = c(.1, .05, .025, .01)
  )
}

map_families <- function(fams) {
  mapping <- c(indep = 0, gaussian = 1, bb1 = 3, bb8 = 6)
  return(unname(mapping[fams]))
}

extract_params <- function(model, fill = 0) {
  params <- get_all_parameters(model, trees = NA)   # list (by tree) of lists (by edge) of matrices
  mats   <- do.call(c, params)                      # flatten to a single list of matrices (one per edge)
  
  p1 <- sapply(mats, function(m) as.numeric(m)[1])  # first parameter
  p2 <- sapply(mats, function(m) {
    v <- as.numeric(m)
    if (length(v) >= 2) v[2] else fill              # second parameter or 0
  })
  
  list(par1 = unname(p1), par2 = unname(p2))
}

extract_rotations <- function(model, as_quadrant = FALSE) {
  pcs  <- model$pair_copulas        # nested list [[tree]][[edge]] of bicop_dist
  cops <- do.call(c, pcs)            # flatten to a single list (length = d(d-1)/2)
  
  rot <- vapply(cops, function(c) c$rotation, numeric(1))
  unname(rot)
}

extend_with_gaussian <- function(vinecop_t1, structure_full) {
  # get dimension from structure (or from T1 length + 1)
  d <- dim(structure_full)[1]
  
  obj <- vinecop_t1                       # keep class 'vinecop'
  obj$structure <- structure_full         # ensure full structure set
  
  # ensure 'pair_copulas' has d-1 trees; keep Tree 1 as-is
  pcs <- obj$pair_copulas
  pcs <- c(pcs, vector("list", (d - 1) - length(pcs)))  # pad if needed
  
  # fill Trees 2..(d-1) with Gaussian placeholders
  for (t in 2:(d - 1)) {
    pcs[[t]] <- replicate(d - t, bicop_dist("gaussian", parameters = 0), simplify = FALSE)
  }
  obj$pair_copulas <- pcs
  obj
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
covar_tail_vec <- function(R_draws, r_p, VaRj, alpha = 0.05, minN = 200) {
  L <- nrow(R_draws); d <- ncol(R_draws)
  out <- numeric(d)
  for (j in seq_len(d)) {
    idx <- which(R_draws[, j] <= VaRj[j])
    if (length(idx) < minN) {
      # fallback: take nearest neighbors around VaR_j if the tail set is too small
      ord <- order(abs(R_draws[, j] - VaRj[j]))
      idx <- ord[1:minN]
    }
    out[j] <- as.numeric(quantile(r_p[idx], probs = alpha, names = FALSE))
  }
  out
}



#####################################################################################




n_assets <- 1:5

data <- list(
  U      = readRDS("data/PIT.rds")[,n_assets],
  mu_fc  = readRDS("data/returns_mean_forecast.rds")[,n_assets+1],# [,-1],  # drop date col
  sig_fc = readRDS("data/returns_vol_forecast.rds")[n_assets+1],  #[,-1],
  df_fc = readRDS("data/df_fc.rds")[,n_assets+1],#[,-1]
  shape_fc = readRDS("data/shape_fc.rds")[,n_assets+1]#[,-1]
)


U         <- data$U
mu_fc     <- data$mu_fc
sig_fc    <- data$sig_fc
df_fc     <- data$df_fc 
shape_fc     <- data$shape_fc 

y_real = readRDS("data/returns_actual.rds")[,n_assets+1]

cfg <- build_cfg(d = ncol(U))
N <- nrow(U); K <- cfg$K; tickers <- colnames(U); A <- length(cfg$alphas)
n_oos <- N - cfg$W_predict; d <- cfg$d; t_train <- cfg$W_predict


out <- list(
  log_pred    = numeric(n_oos),
  fam_hist  = matrix(NA_integer_, n_oos,K),
  par1_hist = matrix(NA_real_, n_oos,K),
  par2_hist = matrix(NA_real_, n_oos,K),
  rotation_hist = matrix(NA_real_, n_oos,K),
  risk    =  list(
    dates = integer(N),                                        
    mean  = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    var   = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    ci_lo = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    ci_hi = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    VaR   = array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),  
  out$QL   <- array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
  out$FZL  <- array(NA_real_, dim = c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
  out$wCRPS <- matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
  port = list(
    dates  = integer(N),
    mean   = numeric(N),
    VaR    = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    ES     = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    QL    <- matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    FZL   <- matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    wCRPS <- numeric(n_oos)
  ),

  
  out$CoVaR_tail <- matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers))
)



for (t in seq_len(n_oos)) {
  test_idx <- t_train + t
  u_train <- U[(test_idx - t_train):(test_idx - 1), , drop = FALSE]
  if (t == 1)  skel <- make_skeleton_CVM(u_train)
  fit_t1 <- vinecop(u_train,
                   family_set = c("indep", "gaussian", "bb1", "bb8"),
                   structure   = skel$structure,
                   allow_rotations = TRUE,
                   trunc_lvl       = 1)
  
  template_vinecop  <- extend_with_gaussian(fit_t1, skel$structure)
  
  model <- vinecop(
    u_train,
    vinecop_object  = template_vinecop,  # <-- must be class 'vinecop'
    par_method      = "mle",
    allow_rotations = FALSE              # rotations irrelevant for Gaussian
  )
  
  out$fam_hist[t, ] <- map_families(unlist(get_all_families(model, trees = NA)))
  out$par1_hist[t, ] <- extract_params(model)$par1
  out$par2_hist[t, ] <- extract_params(model)$par2
  out$rotation_hist[t, ] <- extract_rotations(model)
  
  
  
  ## Predictive density
  u_test            <- U[test_idx, , drop = FALSE]
  out$log_pred[t]   <- dvinecop(u_test, model, cores = 7)
  
  ## Draws & risk metrics 
  draws <- rvinecop(n=5000, model, qrng = FALSE, cores = 7)
  Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])  
  R_t <- sweep(Z_pred, 2, as.numeric(sig_fc[t, ]), `*`) + as.numeric(mu_fc[t, ])
  
  # Risk metrics 
  rs <- risk_stats_full(R_t, cfg$alphas)
    
  out$risk$dates[t]   <- t
  out$risk$VaR [t, , ] <- rs$VaR         
  out$risk$ES  [t, , ] <- rs$ES
  
  
  
  
  y_real_t <- y_real[t,]
  
  # EW-portfolio metrics 
  r_p  <- rowMeans(R_t)                                            
  ps   <- port_stats(r_p, cfg$alphas)                             
  
  out$port$dates[t]   <- t                                     
  out$port$VaR [t, ]  <- ps$VaR
  out$port$ES  [t, ]  <- ps$ES  
  
  r_p_real <- mean(as.numeric(y_real_t))
  out$port$QL[t, ]   <- vapply(seq_along(cfg$alphas), function(k) pinball_loss(r_p_real, ps$VaR[k], cfg$alphas[k]), numeric(1))
  out$port$FZL[t, ]  <- vapply(seq_along(cfg$alphas), function(k) fzl_pzc_scalar(r_p_real, ps$VaR[k], ps$ES[k], cfg$alphas[k]), numeric(1))
  out$port$wCRPS[t]  <- wcrps_gr_scalar(r_p, r_p_real)
  
  
  
  # Quantile loss per asset & alpha using the VaR you already computed
  out$QL[t, , ]  <- pinball_matrix(y_real_t, rs$VaR, cfg$alphas)
  # FZL joint loss for (VaR, ES)
  out$FZL[t, , ] <- fzl_pzc_matrix(y_real_t, rs$VaR, rs$ES, cfg$alphas)
  
  # Weighted CRPS from predictive draws 'R_t' and realization
  # (choose tau=0.10 for left 10% focus; change if you want a different tail focus)
  out$wCRPS[t, ] <- wcrps_gr_matrix(R_t, y_real_t)
  
  
  
  k5 <- which.min(abs(cfg$alphas - 0.05))  # index of alpha = 0.05
  VaRj_5 <- rs$VaR[, k5]     
  out$CoVaR_tail[t, ] <- covar_tail_vec(R_t, r_p, VaRj_5, alpha = 0.05, minN = 20)
  print(t)
  
}
  
  
saveRDS(out, file = file.path("empirical_results", "out_naive_4.rds"))















###########################################################################################


sum(na.omit(out$fam_hist[,4] == 4))/length(na.omit(out$fam_hist[,1]))


  
actual = readRDS("data/returns_actual.rds")[,2:4]
  
dim(actual)
  
port_actual <- rowMeans(actual)#[757:1501] 

length(port_actual)
  
  
out$port$VaR[,1]
  
  
  
length(out$port$VaR[,1])

sum(port_actual<out$port$VaR[,1])/length(port_actual)



sum(actual[757:1501,3]<out$risk$VaR[,3,4])/(1501-757)

actual[757:1501,stock]  
  
# Extract series
actual_series <- port_actual[2757:3500]
var_series <- out$port$VaR[, 1]  # make sure dimensions match

# Plot actual portfolio values
plot(actual_series, type = "l", col = "black", lwd = 2,
     ylab = "Value", xlab = "Time", main = "Actual vs VaR",
     ylim = range(c(actual_series, var_series)))

# Add the VaR line
lines(var_series, col = "red", lwd = 2, lty = 2)

# Optional: Add legend
legend("topright", legend = c("Actual", "VaR"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)











out2 <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/out_naive_3.rds")


order(out$log_pred)


sort(out$log_pred, decreasing = TRUE)
