
# 3.0 transform to 2 dimensions

library(here)
library(fGarch)
source(here("src/R", "config.R"))   

n_sim   <- 50
d       <- 2
n_train <- 252
n_test  <- 3000
n       <- n_train+n_test

fam_names <- c("bb1", "bb1r180", "bb7", "bb7r180", "t")

cvine_struct <- function(d) cvine_structure(order = 1:d)  # C-vine 1-2-3


scale_by_time <- function(mat, v) {
  stopifnot(nrow(mat) == length(v))
  mat * matrix(v, nrow = nrow(mat), ncol = ncol(mat))
}

make_piecewise_dgp_d2 <- function(T, L_switch, draw_fun) {
  n_reg <- ceiling(T / L_switch)
  
  U <- matrix(NA_real_, nrow = T, ncol = 2)
  regime_id <- integer(T)
  regimes <- vector("list", n_reg)
  
  t0 <- 1L
  for (r in seq_len(n_reg)) {
    t1 <- min(T, t0 + L_switch - 1L)
    
    dgp_r <- draw_fun()  # must return list with $vc and $true_bases
    U[t0:t1, ] <- rvinecop(t1 - t0 + 1L, dgp_r$vc)
    
    # store regime-specific truth you will need later
    u_star_r <- c(
      a0.05b0.05 = solve_u_star(dgp_r$vc, alpha = 0.05, beta = 0.05),
      a0.05b0.10 = solve_u_star(dgp_r$vc, alpha = 0.05, beta = 0.10),
      a0.10b0.10 = solve_u_star(dgp_r$vc, alpha = 0.10, beta = 0.10),
      a0.10b0.05 = solve_u_star(dgp_r$vc, alpha = 0.10, beta = 0.05)
    )
    
    regimes[[r]] <- list(
      r = r,
      t_start = t0,
      t_end   = t1,
      idx     = t0:t1,
      vc      = dgp_r$vc,
      true_bases = dgp_r$true_bases,  # for d=2 should be length 1
      u_star  = u_star_r
    )
    
    regime_id[t0:t1] <- r
    t0 <- t1 + 1L
  }
  
  true_base_by_regime <- vapply(regimes, function(g) g$true_bases[1], character(1))
  true_base_t <- true_base_by_regime[regime_id]
  
  list(U = U, regimes = regimes, regime_id = regime_id, true_base_t = true_base_t)
}

# draw one bicop_dist from the allowed set
draw_bicop <- function() {
  f <- sample(fam_names, 1)
  if (f == "t") {
    # target moderate |tau|; map to rho; pick df
    tau <- runif(1, 0.2, 0.6) * 1#sample(c(-1, 1), 1)
    rho <- sin(pi * tau / 2)
    nu  <- runif(1, 2, 7)
    list(name_base = "t",
         bic       = bicop_dist("t", 0, c(rho, nu)))
  } else if (f %in% c("bb1")) {
    # BB1: theta>0, delta>=1  (keep in stable ranges)
    theta <- runif(1, 0.6, 3.0)
    delta <- runif(1, 1.1, 2.5)
    rot   <- 0
    list(name_base = "bb1",
         bic       = bicop_dist("bb1", rot, c(theta, delta)))  
  } else if (f %in% c("bb1r180")) {
    # BB1: theta>0, delta>=1  (keep in stable ranges)
    theta <- runif(1, 0.6, 3.0)
    delta <- runif(1, 1.1, 2.5)
    rot   <- ifelse(f == "bb1r180", 180, 0)
    list(name_base = "bb1r180",
         bic       = bicop_dist("bb1", rot, c(theta, delta)))
  } else if (f %in% c("bb7")) {
    # BB7: theta>=1, delta>0
    theta <- runif(1, 1.1, 3.0)
    delta <- runif(1, 0.6, 2.5)
    rot   <- 0
    list(name_base = "bb7",
         bic       = bicop_dist("bb7", rot, c(theta, delta)))
  } else if (f %in% c("bb7r180")) {
    # BB7: theta>=1, delta>0
    theta <- runif(1, 1.1, 3.0)
    delta <- runif(1, 0.6, 2.5)
    rot   <- ifelse(f == "bb7r180", 180, 0)
    list(name_base = "bb7r180",
         bic       = bicop_dist("bb7", rot, c(theta, delta)))
  } else stop("unknown family requested")
}


# build a 3D C-vine with random edge families/params
draw_vine_d3 <- function() {
  e12   <- draw_bicop()
  e13   <- draw_bicop()
  e23_1 <- draw_bicop()
  list(
    vc = vinecop_dist(
      pair_copulas = list(
        list(e12$bic, e13$bic),  # tree 1
        list(e23_1$bic)          # tree 2
      ),
      structure = cvine_structure(3:1)  # <-- full C-vine for d=3
    ),
    true_bases = c(e12$name_base, e13$name_base, e23_1$name_base),
    pair_list  = list(e12=e12, e13=e13, e23_1=e23_1)
  )
}

# clamp01_2 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)
# 
# # This is now your draw_bicop()
# draw_bicop <- function(
#     td_ranges = list(
#       # calm (symmetric, mild)
#       t   = c(0.02, 0.12),
# 
#       # contenders in the bulk: mild tails (both small)
#       bb1 = list(L = c(0.02, 0.15), U = c(0.02, 0.15)),
#       bb7 = list(L = c(0.02, 0.15), U = c(0.02, 0.15)),
# 
#       # stress: strong lower tail, weak upper tail (crash dependence)
#       bb1r180 = list(L = c(0.40, 0.80), U = c(0.01, 0.08)),
#       bb7r180 = list(L = c(0.40, 0.80), U = c(0.01, 0.08))
#     ),
#     nu_range = c(5, 30),   # calmer t when needed
#     eps = 1e-6
# ) {
#   # uses global fam_names like your original code
#   f <- sample(fam_names, 1)
# 
#   if (f == "t") {
#     lam <- clamp01_2(runif(1, td_ranges$t[1], td_ranges$t[2]), eps)
#     nu  <- runif(1, nu_range[1], nu_range[2])
#     rho <- t_tail2rho(lam, nu)
#     bic <- bicop_dist("t", 0, c(rho, nu))
#     return(list(name_base = "t", bic = bic, lambdaL = lam, lambdaU = lam, par = c(rho, nu)))
#   }
# 
#   draw_lams <- function(key) {
#     lamL <- clamp01_2(runif(1, td_ranges[[key]]$L[1], td_ranges[[key]]$L[2]), eps)
#     lamU <- clamp01_2(runif(1, td_ranges[[key]]$U[1], td_ranges[[key]]$U[2]), eps)
#     list(lamL = lamL, lamU = lamU)
#   }
# 
#   if (f == "bb1") {
#     l   <- draw_lams("bb1")
#     par <- bb1_tail2par(l$lamL, l$lamU)
#     bic <- bicop_dist("bb1", 0, par)
#     return(list(name_base = "bb1", bic = bic, lambdaL = l$lamL, lambdaU = l$lamU, par = par))
#   }
# 
#   if (f == "bb1r180") {
#     l   <- draw_lams("bb1r180")
#     par <- bb1_tail2par(l$lamL, l$lamU)
#     bic <- bicop_dist("bb1", 180, par)
#     return(list(name_base = "bb1r180", bic = bic, lambdaL = l$lamL, lambdaU = l$lamU, par = par))
#   }
# 
#   if (f == "bb7") {
#     l   <- draw_lams("bb7")
#     par <- bb7_tail2par(l$lamL, l$lamU)
#     bic <- bicop_dist("bb7", 0, par)
#     return(list(name_base = "bb7", bic = bic, lambdaL = l$lamL, lambdaU = l$lamU, par = par))
#   }
# 
#   if (f == "bb7r180") {
#     l   <- draw_lams("bb7r180")
#     par <- bb7_tail2par(l$lamL, l$lamU)
#     bic <- bicop_dist("bb7", 180, par)
#     return(list(name_base = "bb7r180", bic = bic, lambdaL = l$lamL, lambdaU = l$lamU, par = par))
#   }
# 
#   stop("unknown family requested")
# }


draw_vine_d2 <- function() {
  e12   <- draw_bicop()
  list(
    vc = vinecop_dist(
      pair_copulas = list(
        list(e12$bic)  # tree 1

      ),
      structure = cvine_structure(2:1)  # <-- full C-vine for d=3
    ),
    true_bases = c(e12$name_base),
    pair_list  = list(e12=e12)
  )
}





test_loglik <- function(model, Utest) sum(na.omit(log(dvinecop(Utest, model))))

q_sstd <- function(u) qsstd(u, mean = 0, sd = 1, nu = 3, xi = 0.9)

# solve u* from C(beta, u*) = alpha*beta using the TRUE bivariate copula
solve_u_star <- function(bic, alpha, beta) {
  target <- alpha * beta
  f <- function(u) pbicop(cbind(beta, u), family = bic$pair_copulas[[1]][[1]]$family, 
                          rotation = bic$pair_copulas[[1]][[1]]$rotation, 
                          parameters = bic$pair_copulas[[1]][[1]]$parameters )  - target
  uniroot(f, interval = c(1e-12, 1 - 1e-12))$root
}

############
rmse <- function(x, y, na.rm = TRUE) {
  if (na.rm) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
  }
  sqrt(mean((x - y)^2))
}

mae <- function(x, y, na.rm = TRUE) {
  if (na.rm) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
  }
  mean(abs(x - y))
}

#------------------------------------------------------------
# Compare out$CoVaR_tail (n_test x 2 x 4) against truth series
# - Uses only scenarios: a0.05b0.05 -> covar_true_5
#                         a0.1b0.1  -> covar_true_10
# - Aligns by taking last n_test points from truth (since truth includes train)
#------------------------------------------------------------
covar_rmse_mae <- function(out, dgp,
                           n_test = dim(out$CoVaR_tail)[1],
                           scen_names = dimnames(out$CoVaR_tail)[[3]],
                           map = list(
                             a0.05b0.05 = "covar_true_5",
                             a0.1b0.1   = "covar_true_10"
                           ),
                           na.rm = TRUE) {
  
  # --- sanity checks
  if (is.null(out$CoVaR_tail)) stop("out$CoVaR_tail not found.")
  if (length(dim(out$CoVaR_tail)) != 3) stop("out$CoVaR_tail must be a 3D array: (time x cond_asset x scenario).")
  if (dim(out$CoVaR_tail)[2] != 2) stop("Expected 2 conditioning assets in out$CoVaR_tail (dim[,2] = 2).")
  if (is.null(scen_names)) stop("Scenario names missing: add dimnames(out$CoVaR_tail)[[3]].")
  
  # build output holder
  res <- data.frame(
    scenario = character(0),
    truth_obj = character(0),
    cond_asset = integer(0),
    RMSE = numeric(0),
    MAE  = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # loop over selected scenarios
  for (sc in names(map)) {
    truth_name <- map[[sc]]
    
    if (!(sc %in% scen_names)) {
      stop(sprintf("Scenario '%s' not found in dimnames(out$CoVaR_tail)[[3]].", sc))
    }
    if (is.null(dgp$truth[[truth_name]])) {
      stop(sprintf("dgp$truth$%s not found.", truth_name))
    }
    
    truth_mat <- dgp$truth[[truth_name]]
    if (!is.matrix(truth_mat) || ncol(truth_mat) != 2) {
      stop(sprintf("dgp$truth$%s must be a matrix with 2 columns (assets).", truth_name))
    }
    
    # align: use last n_test rows from truth
    if (nrow(truth_mat) < n_test) {
      stop(sprintf("Truth series '%s' has only %d rows but need at least n_test=%d.",
                   truth_name, nrow(truth_mat), n_test))
    }
    truth_last <- truth_mat[(nrow(truth_mat) - n_test + 1):nrow(truth_mat), , drop = FALSE]
    
    # predicted slice for this scenario: (n_test x 2 cond_assets)
    pred_slice <- out$CoVaR_tail[, , sc, drop = FALSE][, , 1]
    
    # compute errors for each conditioning asset (dimension 2)
    for (j in 1:2) {
      pred_j  <- pred_slice[, j]
      truth_j <- truth_last[, j]
      
      res <- rbind(res, data.frame(
        scenario = sc,
        truth_obj = truth_name,
        cond_asset = j,
        RMSE = rmse(pred_j, truth_j, na.rm = na.rm),
        MAE  = mae(pred_j, truth_j, na.rm = na.rm),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  res
}


covar_rmse_mae_all <- function(out, dgp,
                           n_test = dim(out$CoVaR_tail)[1],
                           scen_names = dimnames(out$CoVaR_tail)[[3]],
                           map = NULL,
                           na.rm = TRUE) {
  
  # --- sanity checks
  if (is.null(out$CoVaR_tail)) stop("out$CoVaR_tail not found.")
  if (length(dim(out$CoVaR_tail)) != 3) stop("out$CoVaR_tail must be a 3D array: (time x cond_asset x scenario).")
  if (dim(out$CoVaR_tail)[2] != 2) stop("Expected 2 conditioning assets in out$CoVaR_tail (dim[,2] = 2).")
  if (is.null(scen_names)) stop("Scenario names missing: add dimnames(out$CoVaR_tail)[[3]].")
  
  # --- default mapping: ALL FOUR scenarios
  # assumes you stored truth as:
  # dgp$truth$covar_true_a05b05, covar_true_a05b10, covar_true_a10b10, covar_true_a10b05
  if (is.null(map)) {
    map <- list(
      a0.05b0.05 = "covar_true_a05b05",
      a0.05b0.1  = "covar_true_a05b10",
      a0.1b0.1   = "covar_true_a10b10",
      a0.1b0.05  = "covar_true_a10b05"
    )
  }
  
  # build output holder
  res <- data.frame(
    scenario   = character(0),
    truth_obj  = character(0),
    cond_asset = integer(0),
    RMSE       = numeric(0),
    MAE        = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # loop over selected scenarios
  for (sc in names(map)) {
    truth_name <- map[[sc]]
    
    if (!(sc %in% scen_names)) {
      stop(sprintf("Scenario '%s' not found in dimnames(out$CoVaR_tail)[[3]]. Available: %s",
                   sc, paste(scen_names, collapse = ", ")))
    }
    if (is.null(dgp$truth[[truth_name]])) {
      stop(sprintf("dgp$truth$%s not found.", truth_name))
    }
    
    truth_mat <- dgp$truth[[truth_name]]
    if (!is.matrix(truth_mat) || ncol(truth_mat) != 2) {
      stop(sprintf("dgp$truth$%s must be a matrix with 2 columns (assets).", truth_name))
    }
    
    # align: use last n_test rows from truth
    if (nrow(truth_mat) < n_test) {
      stop(sprintf("Truth series '%s' has only %d rows but need at least n_test=%d.",
                   truth_name, nrow(truth_mat), n_test))
    }
    truth_last <- truth_mat[(nrow(truth_mat) - n_test + 1):nrow(truth_mat), , drop = FALSE]
    
    # predicted slice for this scenario: (n_test x 2 cond_assets)
    pred_slice <- out$CoVaR_tail[, , sc, drop = FALSE][, , 1]
    
    # compute errors for each conditioning asset (dimension 2)
    for (j in 1:2) {
      pred_j  <- pred_slice[, j]
      truth_j <- truth_last[, j]
      
      res <- rbind(res, data.frame(
        scenario   = sc,
        truth_obj  = truth_name,
        cond_asset = j,
        RMSE       = rmse(pred_j, truth_j, na.rm = na.rm),
        MAE        = mae(pred_j, truth_j, na.rm = na.rm),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  res
}



#############


var_list <- vector("list", n_sim)
covar_list <- vector("list", n_sim)
logpred_list <- vector("list", n_sim) 
QL_list <- vector("list", n_sim)
FZL_list <- vector("list", n_sim)
wCRPS_list <- vector("list", n_sim)
rmse_mae_from_covar_list <- vector("list", n_sim)


alphas_covar <- c(0.05, 0.10)
betas_covar  <- c(0.05, 0.10)


################################################################################
# Regimes

for (s in 1:(n_sim)) { #################
  set.seed(1111 + s)
  
  Ttot <- n_train + n_test
  L_switch <- 100L  # <- choose your "every x observations" here (or put inside cfg)
  
  piece <- make_piecewise_dgp_d2(T = Ttot, L_switch = L_switch, draw_fun = draw_vine_d2)
  
  # "dgp" now becomes a container for the full regime path + truth
  dgp <- list(
    regimes     = piece$regimes,
    regime_id   = piece$regime_id,
    true_base_t = piece$true_base_t
  )
  
  U <- piece$U
  
  U_transformed <- qsstd(U, mean = 0, sd = 1, nu = 3, xi = 0.9)
  mu_fc    <- matrix(0,  nrow = Ttot, ncol = d)
  sig_fc   <- matrix(NA, nrow = Ttot, ncol = d)
  df_fc    <- matrix(rep(3,   d * Ttot), ncol = d)
  shape_fc <- matrix(rep(0.9, d * Ttot), ncol = d)
  
  omega = 1.785714e-05; alpha = 0.05; beta = 0.90
  
  simulate_garch <- function(n, initial_vol, mu_fc, U_transformed,
                             omega = 1.785714e-05, alpha = 0.05, beta = 0.9) {
    N <- n
    h <- matrix(NA, nrow = n, ncol = d)
    r <- matrix(NA, nrow = n, ncol = d)
    h[1, ] <- omega / (1 - alpha - beta)
    r[1, ] <- sqrt(h[1, ]) * U_transformed[1, ]
    for (t in 2:N) {
      h[t, ] <- omega + alpha * r[t - 1, ]^2 + beta * h[t - 1, ]
      r[t, ] <- sqrt(h[t, ]) * U_transformed[t, ]
    }
    return(list(returns = r, sig = h))
  }
  
  garch_result <- simulate_garch(Ttot, initial_vol, mu_fc, U_transformed)
  simulated_returns <- garch_result$returns
  h <- garch_result$sig
  sig_fc <- matrix(h, ncol = d)
  
  data <- list(
    U        = U,
    mu_fc    = mu_fc,
    sig_fc   = sig_fc,
    df_fc    = df_fc,
    shape_fc = shape_fc,
    y_real   = simulated_returns
  )
  
  # --- Build time-aligned truth (piecewise by regime) ---
  # u_star changes with the true copula, so we store u_star_t over time
  u_star_t <- matrix(NA_real_, nrow = Ttot, ncol = 4,
                     dimnames = list(NULL, c("a0.05b0.05","a0.05b0.10","a0.10b0.10","a0.10b0.05")))
  
  
  for (g in dgp$regimes) {
    u_star_t[g$idx, "a0.05b0.05"] <- g$u_star["a0.05b0.05"]
    u_star_t[g$idx, "a0.05b0.10"] <- g$u_star["a0.05b0.10"]
    u_star_t[g$idx, "a0.10b0.10"] <- g$u_star["a0.10b0.10"]
    u_star_t[g$idx, "a0.10b0.05"] <- g$u_star["a0.10b0.05"]
  }
  
  
  q_beta <- setNames(q_sstd(betas_covar), paste0("b", betas_covar))
  q_u_t <- apply(u_star_t, 2, q_sstd)
  
  covar_true_a05b05 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.05b0.05"])
  covar_true_a05b10 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.05b0.10"])
  covar_true_a10b10 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.10b0.10"])
  covar_true_a10b05 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.10b0.05"])
  
  
  dgp$truth <- list(
    regimes = lapply(dgp$regimes, function(g) {
      list(
        r = g$r, t_start = g$t_start, t_end = g$t_end,
        true_bases = g$true_bases,
        u_star = g$u_star
      )
    }),
    u_star_t     = u_star_t,
    alphas_covar = alphas_covar,
    betas_covar  = betas_covar,
    var_true_10  = sqrt(data$sig_fc) * q_beta[2],
    var_true_5   = sqrt(data$sig_fc) * q_beta[1],
    
    covar_true_a05b05 = covar_true_a05b05,
    covar_true_a05b10 = covar_true_a05b10,
    covar_true_a10b10 = covar_true_a10b10,
    covar_true_a10b05 = covar_true_a10b05
  )
  
  
  cfg <- modifyList(build_cfg(d = 2), list(M = 500, label = "M500"))
  # cfg$L_switch <- L_switch  # optional, if you want it accessible inside smc_simul()
  
  # --- SMC (unchanged call signature) ---
  # out <- smc_simul(data, cfg, dgp)
  # 
  # n_oos <- nrow(data$U) - cfg$W_predict
  # y_real_oos  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  # rp_real_oos <- rowMeans(y_real_oos)
  # 
  # eval_var <- port_var_backtest_df(out, rp_real_oos, cfg$alphas)
  # eval_covar <- covar_backtest_grid(
  #   out         = out,
  #   y_real_oos  = y_real_oos,
  #   rp_real_oos = rp_real_oos,
  #   cfg_alphas  = cfg$alphas,
  #   alphas_eval = c(0.10, 0.05),
  #   d           = d
  # )
  # rmse_mae_from_covar <- covar_rmse_mae(out, dgp)
  
  out <- naive_simul_d2_regimes(data, cfg, dgp)
  #out2 <- naive_simul_d2_regimes(data, cfg, dgp)

  n_oos <- nrow(data$U) - cfg$W_predict
  y_real_oos  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  rp_real_oos <- rowMeans(y_real_oos)

  eval_var <- port_var_backtest_df(out, rp_real_oos, cfg$alphas)
  eval_covar <- covar_backtest_grid(
    out         = out,
    y_real_oos  = y_real_oos,
    rp_real_oos = rp_real_oos,
    cfg_alphas  = cfg$alphas,
    alphas_eval = c(0.10, 0.05),
    d           = d,
  )

  rmse_mae_from_covar <- covar_rmse_mae_all(out, dgp)
  
  
  
  sim_id <- s  # or whatever your simulation index is
  eval_var$sim <- sim_id
  eval_covar$sim <- sim_id
  
  var_list[[sim_id]]  <- eval_var
  covar_list[[sim_id]] <- eval_covar
  logpred_list[[s]] <- out$log_pred
  QL_list[[s]] <- apply(out$port$QL, 2, sum)
  FZL_list[[s]] <- apply(out$port$FZL, 2, sum)
  wCRPS_list[[s]] <- sum(out$port$wCRPS)
  rmse_mae_from_covar_list[[s]] <- rmse_mae_from_covar
  
  out <- list(
    var_list     = var_list,
    covar_list   = covar_list,
    logpred_list = logpred_list,
    QL_list      = QL_list,
    FZL_list     = FZL_list,
    wCRPS_list   = wCRPS_list,
    rmse_mae_from_covar_list = rmse_mae_from_covar_list
  )
  out_dir <- "simul_results/2d_2"
  out_file <- file.path(out_dir, sprintf("results_naive_regimes_tail1_refit252_s%03d.rds", s))
  
  saveRDS(out, out_file)
}






################################################################################

for (s in 1:(n_sim)) {
  set.seed(1111 + s)
  dgp <- draw_vine_d2()
  true_bases <- dgp$true_bases
  U   <- rvinecop(n_train + n_test, dgp$vc)
          
  U_transformed <- qsstd(U, mean = 0, sd = 1, nu = 3, xi = 0.9)
  mu_fc   <- matrix(0, nrow = n_train + n_test, ncol = d)
  sig_fc  <- matrix(NA, nrow = n_train + n_test, ncol = d)
  df_fc   <- matrix(rep(3, d*(n_train + n_test)), ncol = d)  # Degrees of freedom (fixed at 3)
  shape_fc<- matrix(rep(0.9, d*(n_train + n_test)), ncol = d)  # Shap
  
  omega = 1.785714e-05; alpha = 0.05; beta = 0.90
  simulate_garch <- function(n, initial_vol, mu_fc, U_transformed, omega = 1.785714e-05, alpha = 0.05, beta = 0.9) {
    N <- n 
    h <- matrix(NA, nrow = n, ncol = d)
    r <- matrix(NA, nrow = n, ncol = d)
    h[1, ] <- omega / (1 - alpha - beta)  # Initial volatility (unconditional volatility)
    r[1, ] <- sqrt(h[1, ]) * U_transformed[1, ]  # First 
    for (t in 2:N) {
      h[t, ] <- omega + alpha * r[t - 1, ]^2 + beta * h[t - 1, ]
      r[t, ] <- sqrt(h[t, ]) * U_transformed[t, ]
    }
    return(list(returns = r, sig = h))
  }
  # Simulate returns and volatilities
  garch_result <- simulate_garch(n_train + n_test, initial_vol, mu_fc, U_transformed)  
  # Extract simulated returns and volatilities
  simulated_returns <- garch_result$returns
  h <- garch_result$sig
  sig_fc <- matrix(h, ncol = d)  # Assign var to sig_fc for each dimension
  data <- list()
  # Store transformed U and other variables
  data$U <- U
  data$mu_fc <- mu_fc
  data$sig_fc <- sig_fc
  data$df_fc <- df_fc
  data$shape_fc <- shape_fc
  data$y_real <- simulated_returns
  
  #NEW
   u_star <- c(
     a0.05b0.05 = solve_u_star(dgp$vc, alpha = 0.05, beta = 0.05),
     a0.10b0.10 = solve_u_star(dgp$vc, alpha = 0.10, beta = 0.10)
   )

   q_beta = setNames(q_sstd(betas_covar), paste0("b", betas_covar))
   q_u <- setNames(q_sstd(u_star), names(u_star))
   dgp$truth <- list(
     u_star = u_star,
     alphas_covar = alphas_covar,
     betas_covar  = betas_covar,
     var_true_10 = sqrt(data$sig_fc) * q_beta[2],
     var_true_5 = sqrt(data$sig_fc) * q_beta[1],
     covar_true_10 = sqrt(data$sig_fc) * q_u[2],
     covar_true_5 = sqrt(data$sig_fc) * q_u[1]
   )
  
  
  cfg <- modifyList(build_cfg(d = 2), list(M = 500, label = "M500"))
  
  # # # NAIVE
  # out <- naive_simul_d2(data, cfg, dgp)
  # 
  # n_oos <- nrow(data$U) - cfg$W_predict
  # y_real_oos  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  # rp_real_oos <- rowMeans(y_real_oos)
  # 
  # eval_var <- port_var_backtest_df(out, rp_real_oos, cfg$alphas)
  # eval_covar <- covar_backtest_grid(
  #   out         = out,
  #   y_real_oos  = y_real_oos,
  #   rp_real_oos = rp_real_oos,
  #   cfg_alphas  = cfg$alphas,
  #   alphas_eval = c(0.10, 0.05),
  #   d           = d,
  # )
  # 
  # rmse_mae_from_covar <- covar_rmse_mae(out, dgp)

  #SMC
  out <- smc_simul(data, cfg, dgp)
  n_oos <- nrow(data$U) - cfg$W_predict
  y_real_oos  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  rp_real_oos <- rowMeans(y_real_oos)

  eval_var <- port_var_backtest_df(out, rp_real_oos, cfg$alphas)
  eval_covar <- covar_backtest_grid(
    out         = out,
    y_real_oos  = y_real_oos,
    rp_real_oos = rp_real_oos,
    cfg_alphas  = cfg$alphas,
    alphas_eval = c(0.10, 0.05),
    d           = d,
  )
  rmse_mae_from_covar <- covar_rmse_mae(out, dgp)
  
  # TRUE MODEL
  
  # out_true <- tryCatch(
  #   true_model_simul(data, cfg, dgp),
  #   error = function(e) NULL
  # )
  # if (is.null(out_true)) next
  # n_oos <- nrow(data$U) - cfg$W_predict
  # y_real_oos  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  # rp_real_oos <- rowMeans(y_real_oos)
  # 
  # eval_var_true <- port_var_backtest_df(out_true, rp_real_oos, cfg$alphas)
  # eval_covar_true <- covar_backtest_grid(
  #   out         = out_true,
  #   y_real_oos  = y_real_oos,
  #   rp_real_oos = rp_real_oos,
  #   cfg_alphas  = cfg$alphas,
  #   alphas_eval = c(0.10, 0.05),
  #   d           = d,
  # )

  
  
  sim_id <- s  # or whatever your simulation index is
  eval_var$sim <- sim_id
  eval_covar$sim <- sim_id

  var_list[[sim_id]]  <- eval_var
  covar_list[[sim_id]] <- eval_covar
  logpred_list[[s]] <- out$log_pred
  QL_list[[s]] <- apply(out$port$QL, 2, sum)
  FZL_list[[s]] <- apply(out$port$FZL, 2, sum)
  wCRPS_list[[s]] <- sum(out$port$wCRPS)
  rmse_mae_from_covar_list[[s]] <- rmse_mae_from_covar
  
  out <- list(
    var_list     = var_list,
    covar_list   = covar_list,
    logpred_list = logpred_list,
    QL_list      = QL_list,
    FZL_list     = FZL_list,
    wCRPS_list   = wCRPS_list,
    rmse_mae_from_covar_list = rmse_mae_from_covar_list
  )
  out_dir <- "simul_results/2d_2"
  out_file <- file.path(out_dir, sprintf("results_smc_s%03d.rds", s))
  
  saveRDS(out, out_file)
}









## SAVE
out <- list(
  var_list     = var_list,
  covar_list   = covar_list,
  logpred_list = logpred_list,
  QL_list      = QL_list,
  FZL_list     = FZL_list,
  wCRPS_list   = wCRPS_list
)

#saveRDS(out, file = file.path(, "results.rds"))
saveRDS(var_list, "simul_results/2d/var_list_smc.rds")
saveRDS(covar_list, "simul_results/2d/covar_list_smc.rds")
  
  

var_all  <- rbindlist(var_list,  fill = TRUE)
covar_all <- rbindlist(covar_list, fill = TRUE)

scores <- do.call(rbind, records)
print(scores)
# mean(scores$ll_mle); if ("ll_smc" %in% names(scores)) mean(scores$ll_smc)




res <- list()
s=50

for (i in 1:s) {
  df <- xxx$covar_list[[i]]
  if (is.null(df)) next
  
  sub <- df[df$asset %in% c(1,2) &
              abs(df$alpha_j - 0.1) < 1e-12 &
              abs(df$alpha_port - 0.1) < 1e-12, , drop = FALSE]
  
  if (nrow(sub) == 0) next
  sub$sim <- i
  res[[length(res) + 1]] <- sub
}

res_01_01 <- do.call(rbind, res)

mean(res_01_01[res_01_01$asset==1,,drop=FALSE]$rate)
mean(res_01_01[res_01_01$asset==2,,drop=FALSE]$rate)


res <- list()

for (i in 1:s) {
  df <- xxx$covar_list[[i]]
  if (is.null(df)) next
  
  sub <- df[df$asset %in% c(1,2) &
              abs(df$alpha_j - 0.05) < 1e-12 &
              abs(df$alpha_port - 0.05) < 1e-12, , drop = FALSE]
  
  if (nrow(sub) == 0) next
  sub$sim <- i
  res[[length(res) + 1]] <- sub
}

res_01_01 <- do.call(rbind, res)

mean(res_01_01[res_01_01$asset==1,,drop=FALSE]$rate)
mean(res_01_01[res_01_01$asset==2,,drop=FALSE]$rate)





res_01_01 <- do.call(rbind, lapply(seq_along(covar_list), function(i) {
  df <- covar_list[[i]]
  if (is.null(df)) return(NULL)
  
  
  
  if (nrow(sub) == 0) return(NULL)
  
  sub$sim <- i
  sub
}))

res_01_01







out_dir <- "simul_results/d2_2"
s <- 30  # set this

res_list <- vector("list", s)
x <- rep(0, 30)

for (s in seq_len(s)) {
  out_file <- file.path(out_dir, sprintf("results_smc_s%03d.rds", s))
  res <- readRDS(out_file)
  x[s] <- res$QL_list

}



out_file <- file.path(out_dir, sprintf("results_naive_s%03d.rds", s))
res <- readRDS(out_file)data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC




xx <- matrix(0,s,8)
for (s in seq_len(s)) {
  xx[s,] <- xxx$rmse_mae_from_covar_list[[s]]$RMSE
}
colMeans(xx)


xx <- matrix(0,s,8)
for (s in seq_len(s)) {
  xx[s,] <- ccc$rmse_mae_from_covar_list[[s]]$RMSE
}
colMeans(xx)




out_file <- file.path(out_dir, sprintf("results_smc_s%03d.rds", s))
res <- readRDS(out_file)

xx_smc <- rep(0, s)
for (s in seq_len(s)) {
  xx_smc[s] <- sum(na.omit(res$logpred_list[[s]]))
}

mean(xx_smc)

