library(here)
source(here("src/R", "config.R"))     

tip_means <- function(data_up_to_t, cfg) {
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


tip_means_evt <- function(data_up_to_t, cfg) {
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
      emps      <- emp_tails_FRAPO(uv, method = "EVT", k = cfg$tip_k)
      out[[e]]  <- c(mL = as.numeric(emps["L"]), mU = as.numeric(emps["U"]))
    }
  }
  out
}

n_assets <- 1:2

data <- list(
  U      = readRDS("data/PIT.rds")[,n_assets],
  mu_fc  = readRDS("data/returns_mean_forecast.rds")[,(n_assets+1), with = FALSE],# [,-1],  # drop date col
  sig_fc = readRDS("data/returns_vol_forecast.rds")[,(n_assets+1), with = FALSE],  #[,-1],
  df_fc = readRDS("data/df_fc.rds")[,(n_assets+1), with = FALSE],#[,-1]
  shape_fc = readRDS("data/shape_fc.rds")[,(n_assets+1), with = FALSE],
  y_real = readRDS("data/returns_actual.rds")[,n_assets+1, with = FALSE]
  
)
dat = data

U      <- data$U
mu_fc  <- data$mu_fc
sig_fc <- data$sig_fc
df_fc     <- data$df_fc 
shape_fc     <- data$shape_fc 
y_real <- data$y_real

cfg_variants <- list(
list(label = "tip",  use_tail_informed_prior = TRUE, tip_method = "EmpTC",          
     tip_k = NULL, tip_sd_logit = 0.025, q_flip = 0.2)
)
v=cfg_variants[[1]]
cfg <- modifyList(build_cfg(ncol(dat$U)), v[ setdiff(names(v),"label") ])


t_train <- 756
skeleton <- make_skeleton_CVM(U[1:t_train, ])
cfg <- add_first_tree_map(cfg, skeleton)
N <- nrow(U); K <- cfg$K; tickers <- colnames(U); A <- length(cfg$alphas)
n_oos <- N - cfg$W_predict; d <- cfg$d;


out <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

out_tip <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

out_tip3 <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

out_tip4 <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

out_tip5 <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

out_tip6 <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

out_tip7 <- list(
  risk = list(
    dates = integer(n_oos),
    VaR   = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, c(n_oos, d, A), dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),
  
  CoVaR_tail = array(NA_real_, c(n_oos, d, 2), dimnames = list(NULL, tickers, c("a0.05","a0.10")))
)

for (t in seq_len(n_oos)) {
  idx = t
  
  test_idx <- t_train + t
  u_train <- U[(test_idx - t_train):(test_idx - 1), , drop = FALSE]
  y_real_t <- y_real[t,]
  
  data_up_to_t <- u_train#U[1512:1756, , drop = FALSE]
  
  fit_t1 <- vinecop(data_up_to_t,
                    family_set = c("bb1"),
                    structure   = skeleton$structure,
                    allow_rotations = FALSE,
                    trunc_lvl       = 1)
  params <- fit_t1$pair_copulas[[1]][[1]]$parameters
  #bb1_par2tail(params[1], params[2])
  
  draws_f <- rvinecop(10000, fit_t1, cores = 7)
  Z_pred_f <- st_inv_fast(draws_f, shape_fc[idx, ], df_fc[idx, ])  
  Rt_f  <- sweep(Z_pred_f, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_f  <- rowMeans(Rt_f)  
  rs_f <- risk_stats_full(Rt_f, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5_f  <- rs_f$VaR[, k5]   # d-vector
  VaRj_10_f <- rs_f$VaR[, k10]
  covar5_f  <- covar_tail_vec(Rt_f, r_p_f, VaRj_5_f,  port_alpha = 0.05, minN = 50)
  covar10_f <- covar_tail_vec(Rt_f, r_p_f, VaRj_10_f, port_alpha = 0.10, minN = 50)
  
  
  out$risk$dates[t]   <- t
  out$risk$VaR [t, , ] <- rs_f$VaR         
  out$risk$ES  [t, , ] <- rs_f$ES
  out$CoVaR_tail[t, , "a0.05"] <- covar5_f
  out$CoVaR_tail[t, , "a0.10"] <- covar10_f
  
  
  
  #option 2
  
  tip <- tip_means(data_up_to_t, cfg)
  #tip[[1]]["mL"] <- 0.90
  #tip[[1]]["mU"] <- 0.90
  
  pars <- bb1_tail2par(as.numeric(tip[[1]]["mL"]), as.numeric(tip[[1]]["mU"]))
  pars <- sanitize_bb1(pars[1], pars[2])
  bc <- bicop_dist("bb1", parameters = c(pars[1], pars[2]), rotation = 0)
  
  draws_tip <- rbicop(n = 10000, bc) 
  Z_pred_tip <- st_inv_fast(draws_tip, shape_fc[idx, ], df_fc[idx, ])  
  Rt_tip  <- sweep(Z_pred_tip, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_tip  <- rowMeans(Rt_tip)  
  rs_tip <- risk_stats_full(Rt_tip, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5  <- rs_tip$VaR[, k5]   # d-vector
  VaRj_10 <- rs_tip$VaR[, k10]
  covar5  <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar10 <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_10, port_alpha = 0.10, minN = 50)
  
  out_tip$risk$dates[t]   <- t
  out_tip$risk$VaR [t, , ] <- rs_tip$VaR         
  out_tip$risk$ES  [t, , ] <- rs_tip$ES
  out_tip$CoVaR_tail[t, , "a0.05"] <- covar5
  out_tip$CoVaR_tail[t, , "a0.10"] <- covar10
  
  #option 3
  
  tip <- tip_means(data_up_to_t[630:756,], cfg)
  #tip[[1]]["mL"] <- 0.90
  #tip[[1]]["mU"] <- 0.90
  
  pars <- bb1_tail2par(as.numeric(tip[[1]]["mL"]), as.numeric(tip[[1]]["mU"]))
  pars <- sanitize_bb1(pars[1], pars[2])
  bc <- bicop_dist("bb1", parameters = c(pars[1], pars[2]), rotation = 0)
  
  draws_tip <- rbicop(n = 10000, bc) 
  Z_pred_tip <- st_inv_fast(draws_tip, shape_fc[idx, ], df_fc[idx, ])  
  Rt_tip  <- sweep(Z_pred_tip, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_tip  <- rowMeans(Rt_tip)  
  rs_tip <- risk_stats_full(Rt_tip, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5  <- rs_tip$VaR[, k5]   # d-vector
  VaRj_10 <- rs_tip$VaR[, k10]
  covar5  <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar10 <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_10, port_alpha = 0.10, minN = 50)
  
  out_tip3$risk$dates[t]   <- t
  out_tip3$risk$VaR [t, , ] <- rs_tip$VaR         
  out_tip3$risk$ES  [t, , ] <- rs_tip$ES
  out_tip3$CoVaR_tail[t, , "a0.05"] <- covar5
  out_tip3$CoVaR_tail[t, , "a0.10"] <- covar10
  
  #option 4
  
  mU <- taildep(data_up_to_t[,1], data_up_to_t[,2], u = 0.95)["chi"]  # c(chi, chibar)
  mL <- taildep(-data_up_to_t[,1], -data_up_to_t[,2], u = 0.95)["chi"]  # c(chi, chibar)
  
  #tip[[1]]["mL"] <- 0.90
  #tip[[1]]["mU"] <- 0.90
  
  pars <- bb1_tail2par(as.numeric(mL), as.numeric(mU))
  pars <- sanitize_bb1(pars[1], pars[2])
  bc <- bicop_dist("bb1", parameters = c(pars[1], pars[2]), rotation = 0)
  
  draws_tip <- rbicop(n = 10000, bc) 
  Z_pred_tip <- st_inv_fast(draws_tip, shape_fc[idx, ], df_fc[idx, ])  
  Rt_tip  <- sweep(Z_pred_tip, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_tip  <- rowMeans(Rt_tip)  
  rs_tip <- risk_stats_full(Rt_tip, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5  <- rs_tip$VaR[, k5]   # d-vector
  VaRj_10 <- rs_tip$VaR[, k10]
  covar5  <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar10 <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_10, port_alpha = 0.10, minN = 50)
  
  out_tip4$risk$dates[t]   <- t
  out_tip4$risk$VaR [t, , ] <- rs_tip$VaR         
  out_tip4$risk$ES  [t, , ] <- rs_tip$ES
  out_tip4$CoVaR_tail[t, , "a0.05"] <- covar5
  out_tip4$CoVaR_tail[t, , "a0.10"] <- covar10
  
  #option 5
  
  mU <- taildep(data_up_to_t[504:756,1], data_up_to_t[504:756,2], u = 0.95)["chi"]  # c(chi, chibar)
  mL <- taildep(-data_up_to_t[504:756,1], -data_up_to_t[504:756,2], u = 0.95)["chi"]  # c(chi, chibar)
  
  #tip[[1]]["mL"] <- 0.90
  #tip[[1]]["mU"] <- 0.90
  
  pars <- bb1_tail2par(as.numeric(mL), as.numeric(mU))
  pars <- sanitize_bb1(pars[1], pars[2])
  bc <- bicop_dist("bb1", parameters = c(pars[1], pars[2]), rotation = 0)
  
  draws_tip <- rbicop(n = 10000, bc) 
  Z_pred_tip <- st_inv_fast(draws_tip, shape_fc[idx, ], df_fc[idx, ])  
  Rt_tip  <- sweep(Z_pred_tip, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_tip  <- rowMeans(Rt_tip)  
  rs_tip <- risk_stats_full(Rt_tip, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5  <- rs_tip$VaR[, k5]   # d-vector
  VaRj_10 <- rs_tip$VaR[, k10]
  covar5  <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar10 <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_10, port_alpha = 0.10, minN = 50)
  
  out_tip5$risk$dates[t]   <- t
  out_tip5$risk$VaR [t, , ] <- rs_tip$VaR         
  out_tip5$risk$ES  [t, , ] <- rs_tip$ES
  out_tip5$CoVaR_tail[t, , "a0.05"] <- covar5
  out_tip5$CoVaR_tail[t, , "a0.10"] <- covar10
  
  
  #option 5
  
  mU <- taildep(data_up_to_t[,1], data_up_to_t[,2], u = 0.99)["chi"]  # c(chi, chibar)
  mL <- taildep(-data_up_to_t[,1], -data_up_to_t[,2], u = 0.99)["chi"]  # c(chi, chibar)
  
  #tip[[1]]["mL"] <- 0.90
  #tip[[1]]["mU"] <- 0.90
  
  pars <- bb1_tail2par(as.numeric(mL), as.numeric(mU))
  pars <- sanitize_bb1(pars[1], pars[2])
  bc <- bicop_dist("bb1", parameters = c(pars[1], pars[2]), rotation = 0)
  
  draws_tip <- rbicop(n = 10000, bc) 
  Z_pred_tip <- st_inv_fast(draws_tip, shape_fc[idx, ], df_fc[idx, ])  
  Rt_tip  <- sweep(Z_pred_tip, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_tip  <- rowMeans(Rt_tip)  
  rs_tip <- risk_stats_full(Rt_tip, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5  <- rs_tip$VaR[, k5]   # d-vector
  VaRj_10 <- rs_tip$VaR[, k10]
  covar5  <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar10 <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_10, port_alpha = 0.10, minN = 50)
  
  out_tip6$risk$dates[t]   <- t
  out_tip6$risk$VaR [t, , ] <- rs_tip$VaR         
  out_tip6$risk$ES  [t, , ] <- rs_tip$ES
  out_tip6$CoVaR_tail[t, , "a0.05"] <- covar5
  out_tip6$CoVaR_tail[t, , "a0.10"] <- covar10
  
  #option 6
  
  mU <- taildep(data_up_to_t[504:756,1], data_up_to_t[504:756,2], u = 0.99)["chi"]  # c(chi, chibar)
  mL <- taildep(-data_up_to_t[504:756,1], -data_up_to_t[504:756,2], u = 0.99)["chi"]  # c(chi, chibar)
  
  #tip[[1]]["mL"] <- 0.90
  #tip[[1]]["mU"] <- 0.90
  
  pars <- bb1_tail2par(as.numeric(mL), as.numeric(mU))
  pars <- sanitize_bb1(pars[1], pars[2])
  bc <- bicop_dist("bb1", parameters = c(pars[1], pars[2]), rotation = 0)
  
  draws_tip <- rbicop(n = 10000, bc) 
  Z_pred_tip <- st_inv_fast(draws_tip, shape_fc[idx, ], df_fc[idx, ])  
  Rt_tip  <- sweep(Z_pred_tip, 2, as.numeric(sig_fc[idx, ]), `*`) + as.numeric(mu_fc[idx, ])
  r_p_tip  <- rowMeans(Rt_tip)  
  rs_tip <- risk_stats_full(Rt_tip, cfg$alphas)
  k5  <- which.min(abs(cfg$alphas - 0.05))
  k10 <- which.min(abs(cfg$alphas - 0.10))
  VaRj_5  <- rs_tip$VaR[, k5]   # d-vector
  VaRj_10 <- rs_tip$VaR[, k10]
  covar5  <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_5,  port_alpha = 0.05, minN = 50)
  covar10 <- covar_tail_vec(Rt_tip, r_p_tip, VaRj_10, port_alpha = 0.10, minN = 50)
  
  out_tip7$risk$dates[t]   <- t
  out_tip7$risk$VaR [t, , ] <- rs_tip$VaR         
  out_tip7$risk$ES  [t, , ] <- rs_tip$ES
  out_tip7$CoVaR_tail[t, , "a0.05"] <- covar5
  out_tip7$CoVaR_tail[t, , "a0.10"] <- covar10
  print(t)
}















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
    mask <- y_real_oos[, j, with=FALSE] <= VaRj_oos[, j]  # days when j is distressed
    n_list[j] <- sum(mask)
    if (n_list[j] > 0) {
      hits_list[[j]] <- as.numeric(r_p_real[mask] <= CoVaR_oos[mask, j])
    } else {
      hits_list[[j]] <- numeric(0)
    }
  }
  list(hits = hits_list, n = n_list)
}





evaluate <- function(out){

y_real_oos = readRDS("data/returns_actual.rds")[,n_assets+1, with=FALSE]
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





# Tidy summaries ready to print/plot:
#print(eval_port_var)     # list (one per alpha) with stats
#print(eval_asset_var)    # data.frame per asset per alpha
print(eval_covar)      

}



out <- readRDS("empirical_results/explore2_out.rds")
out <- readRDS("empirical_results/explore2_out_tip.rds")


evaluate(out)
evaluate(out_tip3)





#############



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



y_real = readRDS("data/returns_actual.rds")[,1:3, with = FALSE]

df <- cbind(y_real, hj)



#df$cum_hits <- cumsum(df$hj) / seq_along(df$hj)






alphas_eval <- c(0.10, 0.05)
# ==== C) CoVaR backtests (per conditioning stock j, α = 5% and 10%) ====
eval_covar <- do.call(rbind, lapply(alphas_eval, function(a) {
  lab <- if (a == 0.05) "a0.05" else "a0.10"
  k   <- which.min(abs(cfg$alphas - a))
  
  VaRj_oos   <- out$risk$VaR[, , k, drop = FALSE][, , 1]   # n_oos×d
  CoVaR_oos  <- out$CoVaR_tail[, , lab]                    # n_oos×d
  cond_hits  <- covar_hits_by_j(rp_real_oos, y_real_oos, VaRj_oos, CoVaR_oos, alpha = a)
  
  do.call(rbind, lapply(seq_len(d), function(j) {
    co_hj <- cond_hits$hits[[j]]
    if (length(co_hj) == 0) {
      data.frame(asset = tickers[j], alpha = a, T_event = 0,
                 rate = NA, kupiec_p = NA, ind_p = NA, cc_p = NA, row.names = NULL)
    } else {
      data.frame(
        asset    = tickers[j],
        alpha    = a,
        T_event  = length(co_hj),
        rate     = mean(co_hj),
        kupiec_p = kupiec_test(co_hj)$pval,
        ind_p    = christoffersen_ind_test(co_hj)$pval,
        cc_p     = christoffersen_cc_test(co_hj, a)$pval,
        row.names = NULL
      )
    }
  }))
}))




length(co_hj) == sum(df$hj)
df$co_hj <- NA              # create empty column
df$co_hj[df$hj == 1] <- co_hj

df <- df[Date >= as.Date("2021-01-01")]


## indices where co_hj is observed (0/1)
idx <- which(!is.na(df$co_hj))

## cumulative proportion among observed entries only
cum_prop <- cumsum(df$co_hj[idx]) / seq_along(idx)


df$cum_co_prop <- NA_real_
df$cum_co_prop[idx] <- cum_prop

## plot (step line updates only at observed points)
plot(df$Date, df$cum_co_prop, type = "s", lwd = 2,
     xlab = "Date", ylab = "Cumulative proportion (co_hj)",
     main = "Cumulative Proportion of co_hj (conditional on observed)")





