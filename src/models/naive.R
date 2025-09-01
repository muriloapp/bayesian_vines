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
    alphas     = c(.1, .05, .025, .01),
    nc       = max(parallel::detectCores()-1, 1)
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





#####################################################################################




n_assets <- 1:3
n_days <- 1:3000

# data <- list(
#   U      = readRDS("data/PIT.rds")[1:length(n_days),n_assets, with=FALSE],
#   mu_fc  = readRDS("data/returns_mean_forecast.rds")[n_days,n_assets+1, with=FALSE],# [,-1],  # drop date col
#   sig_fc = readRDS("data/returns_vol_forecast.rds")[n_days,n_assets+1, with=FALSE],  #[,-1],
#   df_fc = readRDS("data/df_fc.rds")[n_days,n_assets+1, with=FALSE],#[,-1]
#   shape_fc = readRDS("data/shape_fc.rds")[n_days,n_assets+1, with=FALSE]#[,-1]
# )
# y_real = readRDS("data/returns_actual.rds")[n_days,n_assets+1, with=FALSE]


data <- list(
  U      = readRDS("data/PIT.rds")[,n_assets],
  mu_fc  = readRDS("data/returns_mean_forecast.rds")[,n_assets+1, with=FALSE],# [,-1],  # drop date col
  sig_fc = readRDS("data/returns_vol_forecast.rds")[,n_assets+1, with=FALSE],  #[,-1],
  df_fc = readRDS("data/df_fc.rds")[,n_assets+1, with=FALSE],#[,-1]
  shape_fc = readRDS("data/shape_fc.rds")[,n_assets+1, with=FALSE]#[,-1]
)
y_real = readRDS("data/returns_actual.rds")[,n_assets+1, with=FALSE]


cfg <- build_cfg(d = ncol(U))
N <- nrow(U); K <- cfg$K; tickers <- colnames(U); A <- length(cfg$alphas)
n_oos <- N - cfg$W_predict; d <- cfg$d; t_train <- cfg$W_predict



U         <- data$U
mu_fc     <- data$mu_fc[(.N - n_oos + 1):.N]
sig_fc    <- data$sig_fc[(.N - n_oos + 1):.N]
df_fc     <- data$df_fc[(.N - n_oos + 1):.N]
shape_fc     <- data$shape_fc[(.N - n_oos + 1):.N]
y_real <- y_real[(.N - n_oos + 1):.N]


out <- list(
  log_pred    = numeric(n_oos),
  fam_hist    = matrix(NA_integer_, n_oos, K),
  par1_hist   = matrix(NA_real_, n_oos, K),
  par2_hist   = matrix(NA_real_, n_oos, K),
  rotation_hist = matrix(NA_real_, n_oos, K),
  
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  QL    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
  FZL   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
  wCRPS = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
  
  port = list(
    dates = integer(n_oos),
    VaR   = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    ES    = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    QL    = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    FZL   = matrix(NA_real_, n_oos, A, dimnames = list(NULL, paste0("a", cfg$alphas))),
    wCRPS = numeric(n_oos)
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 4), dimnames = list(NULL, tickers, c("a0.05b0.05","a0.05b0.1","a0.1b0.1","a0.1b0.05")))
)





for (t in seq_len(n_oos)) {
  test_idx <- t_train + t
  u_train <- U[(test_idx - t_train):(test_idx - 1), , drop = FALSE]
  y_real_t <- y_real[t,]
  
  if (t == 1)  skel <- make_skeleton_CVM(u_train)
  fit_t1 <- vinecop(u_train,
                   family_set = c("bb1", "bb7", "t"),
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
  out$log_pred[t]   <- dvinecop(u_test, model, cores = cfg$nc)
  
  ## Draws & risk metrics 
  draws <- rvinecop(n=20000, model, qrng = FALSE, cores = cfg$nc)
  Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])  
  R_t <- sweep(Z_pred, 2, as.numeric(sig_fc[t, ]), `*`) + as.numeric(mu_fc[t, ])
  
  # Risk metrics 
  rs <- risk_stats_full(R_t, cfg$alphas)
    
  out$risk$dates[t]   <- t
  out$risk$VaR [t, , ] <- rs$VaR         
  out$risk$ES  [t, , ] <- rs$ES
  
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
  
  
  
  out$QL[t, , ]  <- pinball_matrix(as.matrix(y_real_t), rs$VaR, cfg$alphas) # Quantile loss per asset & alpha using the VaR you already computed
  out$FZL[t, , ] <- fzl_pzc_matrix(y_real_t, rs$VaR, rs$ES, cfg$alphas) # FZL joint loss for (VaR, ES)
  out$wCRPS[t, ] <- wcrps_gr_matrix(R_t, as.matrix(y_real_t)) # Weighted CRPS from predictive draws 'R_t' and realization
  
  
  # CoVaR
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  
  VaRj_5  <- rs$VaR[, k5]   # d-vector
  VaRj_10 <- rs$VaR[, k10]
  covar5  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar5b10  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.1, minN = 50)
  covar10 <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.10, minN = 50)
  covar10b5 <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.05, minN = 50)
  
  out$CoVaR_tail[t, , "a0.05b0.05"] <- covar5
  out$CoVaR_tail[t, , "a0.05b0.1"] <- covar5b10
  out$CoVaR_tail[t, , "a0.1b0.1"] <- covar10
  out$CoVaR_tail[t, , "a0.1b0.05"] <- covar10b5
  
  
  print(t)
  
}
  
  
saveRDS(out, file = file.path("empirical_results", "out_naive_10.rds"))




###########################################################################################





# Hits for VaR (unconditional): y <= q  (vectors of same length)
var_hits <- function(y, q) as.numeric(y <= q)

# Kupiec (1995) unconditional coverage
kupiec_test <- function(h, alpha) {
  n <- length(h); n1 <- sum(h); n0 <- n - n1
  pi_hat <- n1 / n
  lr_uc <- -2 * ( n0*log1p(-alpha) + n1*log(alpha) - ( n0*log1p(-pi_hat) + n1*log(pi_hat) ) )
  pval  <- 1 - pchisq(lr_uc, df = 1)
  list(stat = lr_uc, pval = pval, rate = pi_hat)
}

# Christoffersen (1998) independence & conditional coverage
christoffersen_ind_test <- function(h) {
  # 2-state Markov chain counts
  hlag <- c(NA, h[-length(h)])
  valid <- !is.na(hlag)
  n00 <- sum(hlag[valid] == 0 & h[valid] == 0)
  n01 <- sum(hlag[valid] == 0 & h[valid] == 1)
  n10 <- sum(hlag[valid] == 1 & h[valid] == 0)
  n11 <- sum(hlag[valid] == 1 & h[valid] == 1)
  n0  <- n00 + n01; n1 <- n10 + n11
  # MLEs
  pi   <- (n01 + n11) / (n0 + n1)
  pi0  <- if (n0 > 0) n01 / n0 else 0
  pi1  <- if (n1 > 0) n11 / n1 else 0
  # LR_ind
  ll_ind <- (n00*log1p(-pi0) + n01*log(pi0) + n10*log1p(-pi1) + n11*log(pi1))
  ll_iid <- ((n0 + n1 - (n01 + n11))*log1p(-pi) + (n01 + n11)*log(pi))
  lr_ind <- -2 * (ll_iid - ll_ind)
  pval   <- 1 - pchisq(lr_ind, df = 1)
  list(stat = lr_ind, pval = pval)
}

christoffersen_cc_test <- function(h, alpha) {
  uc <- kupiec_test(h, alpha); ind <- christoffersen_ind_test(h)
  list(stat = uc$stat + ind$stat, pval = 1 - pchisq(uc$stat + ind$stat, df = 2))
}

# Build conditional hit series for CoVaR^{p|j} at given alpha (same for cond/port here)
# Inputs are OOS arrays/matrices aligned by time:
# r_p_real: n_oos vector (realized portfolio)
# y_real_oos: n_oos×d realized returns
# VaRj_oos:  n_oos×d asset VaR forecasts at alpha (same sign convention)
# CoVaR_oos: n_oos×d CoVaR forecasts at alpha (from your array)
covar_hits_by_j <- function(r_p_real, y_real_oos, VaRj_oos, CoVaR_oos, alpha) {
  d <- ncol(y_real_oos)
  hits_list <- vector("list", d)
  n_list    <- integer(d)
  for (j in seq_len(d)) {
    mask <- y_real_oos[, j] <= VaRj_oos[, j]  # days when j is distressed
    n_list[j] <- sum(mask)
    if (n_list[j] > 0) {
      hits_list[[j]] <- as.numeric(r_p_real[mask] <= CoVaR_oos[mask, j])
    } else {
      hits_list[[j]] <- numeric(0)
    }
  }
  list(hits = hits_list, n = n_list)
}





out$port$VaR <- na.omit(out$port$VaR)
out$CoVaR_tail <- na.omit(out$CoVaR_tail)


y_real_oos = readRDS("data/returns_actual.rds")[1:(length(n_days)-t_train),n_assets+1]


y_real_oos <- y_real[1:3061,]
rp_real_oos <- rowMeans(y_real_oos)



# ==== A) Portfolio VaR backtests for all alphas ====
eval_port_var <- lapply(seq_along(cfg$alphas), function(k) {
  a <- cfg$alphas[k]
  q <- out$port$VaR[, k]
  h <- var_hits(rp_real_oos, q)
  list(
    alpha = a,
    n     = length(h),
    hits  = sum(h),
    rate  = mean(h),
    kupiec = kupiec_test(h, a),
    chr_ind = christoffersen_ind_test(h),
    chr_cc  = christoffersen_cc_test(h, a)
  )
})




# ==== B) Asset VaR backtests (per asset, chosen alphas e.g. 5% & 10%) ====
alphas_eval <- c(0.10, 0.05, 0.025, 0.01)
eval_asset_var <- do.call(rbind, lapply(seq_len(d), function(j) {
  do.call(rbind, lapply(alphas_eval, function(a) {
    k <- which.min(abs(cfg$alphas - a))
    qj <- out$risk$VaR[, j, k]
    hj <- var_hits(y_real_oos[, j, with=FALSE], qj)
    data.frame(
      asset = tickers[j], alpha = a,
      n = length(hj), hits = sum(hj), rate = mean(hj),
      kupiec_p = kupiec_test(hj, a)$pval,
      ind_p    = christoffersen_ind_test(hj)$pval,
      cc_p     = christoffersen_cc_test(hj, a)$pval,
      row.names = NULL
    )
  }))
}))


alphas_eval <- c(0.10, 0.05)
# ==== C) CoVaR backtests (per conditioning stock j, α = 5% and 10%) ====
eval_covar <- do.call(rbind, lapply(alphas_eval, function(a) {
  lab <- if (a == 0.05) "a0.05" else "a0.10"
  k   <- which.min(abs(cfg$alphas - a))
  
  VaRj_oos   <- out$risk$VaR[, , k, drop = FALSE][, , 1]   # n_oos×d
  CoVaR_oos  <- out$CoVaR_tail[, , lab]                    # n_oos×d
  cond_hits  <- covar_hits_by_j(rp_real_oos, y_real_oos, VaRj_oos, CoVaR_oos, alpha = a)
  
  do.call(rbind, lapply(seq_len(d), function(j) {
    hj <- cond_hits$hits[[j]]
    if (length(hj) == 0) {
      data.frame(asset = tickers[j], alpha = a, T_event = 0,
                 rate = NA, kupiec_p = NA, ind_p = NA, cc_p = NA, row.names = NULL)
    } else {
      data.frame(
        asset    = tickers[j],
        alpha    = a,
        T_event  = length(hj),
        rate     = mean(hj),
        kupiec_p = kupiec_test(hj, a)$pval,
        ind_p    = christoffersen_ind_test(hj)$pval,
        cc_p     = christoffersen_cc_test(hj, a)$pval,
        row.names = NULL
      )
    }
  }))
}))



alphas_eval <- c(0.10, 0.05)  # order doesn't matter

# helper to match your dimnames exactly: 0.05 -> "0.05", 0.10 -> "0.1"
fmt_ab <- function(x) ifelse(abs(x - 0.05) < 1e-12, "0.05", 
                             ifelse(abs(x - 0.10) < 1e-12, "0.1", as.character(x)))

grid_ab <- expand.grid(alpha_j = alphas_eval,
                       alpha_port = alphas_eval,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# ==== C) CoVaR backtests (per conditioning stock j, α_j ∈ {5%,10%} & portfolio α_p ∈ {5%,10%}) ====
eval_covar <- do.call(rbind, lapply(seq_len(nrow(grid_ab)), function(i) {
  a <- grid_ab$alpha_j[i]      # VaR threshold for conditioning asset j
  b <- grid_ab$alpha_port[i]   # portfolio CoVaR threshold
  k <- which.min(abs(alphas_eval - a))  # pick matching VaR slice
  
  # labels like "a0.05b0.1", "a0.1b0.05", etc.
  lab <- paste0("a", fmt_ab(a), "b", fmt_ab(b))
  
  # n_oos × d
  VaRj_oos  <- out$risk$VaR[1:3061, , k, drop = FALSE][, , 1]
  CoVaR_oos <- out$CoVaR_tail[1:3061, , lab]  # n_oos × d
  
  # hits when r_p ≤ CoVaR(j; a,b) conditional on asset j being in VaR event at level a
  cond_hits <- covar_hits_by_j(rp_real_oos, as.matrix(y_real_oos), VaRj_oos, CoVaR_oos, alpha = a)
  
  do.call(rbind, lapply(seq_len(d), function(j) {
    hj <- cond_hits$hits[[j]]
    if (length(hj) == 0) {
      data.frame(
        asset      = tickers[j],
        alpha_j    = a,
        alpha_port = b,
        T_event    = 0,
        rate       = NA_real_,
        kupiec_p   = NA_real_,
        ind_p      = NA_real_,
        cc_p       = NA_real_,
        row.names  = NULL
      )
    } else {
      data.frame(
        asset      = tickers[j],
        alpha_j    = a,
        alpha_port = b,
        T_event    = length(hj),
        rate       = mean(hj),
        kupiec_p   = kupiec_test(hj, b)$pval,                 # test vs portfolio level b
        ind_p      = christoffersen_ind_test(hj)$pval,
        cc_p       = christoffersen_cc_test(hj, b)$pval,      # cc vs portfolio level b
        row.names  = NULL
      )
    }
  }))
}))




# Tidy summaries ready to print/plot:
print(eval_port_var)     # list (one per alpha) with stats
print(eval_asset_var)    # data.frame per asset per alpha
print(eval_covar)        # data.frame per conditioning asset per alpha







apply(out$port$FZL, 2, sum)


apply(out$port$QL, 2, sum)


sum(out$port$wCRPS)



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











out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/out_naive_8.rds")


order(out$log_pred)


sort(out$log_pred, decreasing = TRUE)
