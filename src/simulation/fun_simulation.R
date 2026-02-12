

library(here)
library(fGarch)
source(here("src/R", "config.R"))   


alphas_covar <- c(0.025, 0.05, 0.10)
betas_covar  <- c(0.025, 0.05, 0.10)


fam_names <- c("bb1", "bb1r180", "bb7", "bb7r180", "t")

SCEN_COVAR <- c(
  "a0.05b0.05", "a0.05b0.1", "a0.1b0.1", "a0.1b0.05",
  "a0.05b0.025", "a0.025b0.05"
)

cvine_struct <- function(d) cvine_structure(order = 1:d)  # C-vine 1-2-3

scale_by_time <- function(mat, v) {
  stopifnot(nrow(mat) == length(v))
  mat * matrix(v, nrow = nrow(mat), ncol = ncol(mat))
}

draw_regime_lengths <- function(T, mean_len = 100L, min_len = 10L, max_len = 500L) {
  p <- 1 / mean_len  # geometric => E[L] = 1/p
  lens <- integer(0)
  tot <- 0L
  while (tot < T) {
    L <- rgeom(1, prob = p) + 1L
    L <- max(min_len, min(L, max_len))
    lens <- c(lens, as.integer(L))
    tot <- tot + L
  }
  if (tot > T) lens[length(lens)] <- lens[length(lens)] - (tot - T)
  lens
}

assign_states_by_time <- function(lens, p_extreme = 0.20) {
  T <- sum(lens)
  n_reg <- length(lens)
  
  # allow single-regime schedule: EXTREME with prob p_extreme
  if (n_reg == 1L) {
    state <- if (runif(1) < p_extreme) "EXTREME" else "NORMAL"
    return(list(
      states = state,
      frac_extreme_time = if (state == "EXTREME") 1 else 0
    ))
  }
  
  target <- as.integer(round(p_extreme * T))
  
  states <- rep("NORMAL", n_reg)
  ord <- sample.int(n_reg, n_reg)
  
  ext_time <- 0L
  for (k in ord) {
    if (ext_time >= target) break
    states[k] <- "EXTREME"
    ext_time <- ext_time + lens[k]
  }
  
  list(
    states = states,
    frac_extreme_time = sum(lens[states == "EXTREME"]) / T
  )
}

make_regime_schedule <- function(T, mean_len = 100L, p_extreme = 0.20,
                                 min_len = 10L, max_len = 500L) {
  lens <- draw_regime_lengths(T, mean_len = mean_len, min_len = min_len, max_len = max_len)
  st   <- assign_states_by_time(lens, p_extreme = p_extreme)
  list(
    lens = lens,
    states = st$states,
    mean_len_real = mean(lens),
    frac_extreme_time_real = st$frac_extreme_time
  )
}
draw_bicop_regime <- function(state, fam_names) {
  
  f <- sample(fam_names, 1)
  
  is_t   <- (f == "t")
  is_bb1 <- (f %in% c("bb1", "bb1r180"))
  is_bb7 <- (f %in% c("bb7", "bb7r180"))
  
  if (!is_t && !is_bb1 && !is_bb7) stop("Unknown family in fam_names.")
  
  rot <- if (grepl("r180$", f)) 180L else 0L
  fam_base <- if (is_bb1) "bb1" else if (is_bb7) "bb7" else "t"
  
  if (state == "NORMAL") {
  lamL_target <- runif(1, 0.2, 0.85)
  if (is_t) {
    lamU_target <- lamL_target
  } else {
    lamU_target <- runif(1, 0.2, 0.85)
  }
  } else{stop("Check state")}
  
  if (rot == 0L) {
    lamL_base <- lamL_target
    lamU_base <- lamU_target
  } else {
    lamL_base <- lamU_target
    lamU_base <- lamL_target
  }
  
  # invert tails -> parameters
  if (is_t) {
    nu  <- runif(1, 2, 30)
    rho <- abs(t_tail2rho(lamL_target, nu))   # lamL=lamU for t
    par <- c(rho, nu)
    bic <- bicop_dist("t", 0L, par)      # rotation irrelevant for symmetric t
    return(list(
      name_base = "t",
      bic = bic,
      state = state,
      lambdaL = lamL_target,
      lambdaU = lamU_target,
      par = par,
      rot = 0L
    ))
  }
  
  if (is_bb1) {
    par <- bb1_tail2par(lamL_base, lamU_base)
    # BB1 requires delta>=1; ensure numerically safe
    bic <- bicop_dist("bb1", rot, par)
    return(list(
      name_base = f,  # keep bb1 vs bb1r180 in truth labels
      bic = bic,
      state = state,
      lambdaL = lamL_target,
      lambdaU = lamU_target,
      par = par,
      rot = rot
    ))
  }
  
  if (is_bb7) {
    par <- bb7_tail2par(lamL_base, lamU_base)
    bic <- bicop_dist("bb7", rot, par)
    return(list(
      name_base = f,  # keep bb7 vs bb7r180
      bic = bic,
      state = state,
      lambdaL = lamL_target,
      lambdaU = lamU_target,
      par = par,
      rot = rot
    ))
  }
  
  stop("Unreachable.")
}

draw_vine_d2_regime <- function(state) {
  e12 <- draw_bicop_regime(state, fam_names = fam_names)
  
  list(
    vc = vinecop_dist(
      pair_copulas = list(list(e12$bic)),
      structure    = cvine_structure(2:1)
    ),
    true_bases = c(e12$name_base),
    pair_list  = list(e12 = e12),
    state      = state
  )
}

make_piecewise_dgp_d2 <- function(T, L_switch, draw_fun) {
  
  # ---schedule parsing (backward compatible) ---
  if (is.list(L_switch) && !is.null(L_switch$lens) && !is.null(L_switch$states)) {
    lens   <- as.integer(L_switch$lens)
    states <- as.character(L_switch$states)
    if (sum(lens) != T) stop("Schedule lens must sum to T.")
    if (length(states) != length(lens)) stop("Schedule states must have same length as lens.")
  } else {
    L_fixed <- as.integer(L_switch)
    n_reg <- ceiling(T / L_fixed)
    lens <- rep(L_fixed, n_reg)
    lens[length(lens)] <- T - sum(lens[-length(lens)])
    states <- rep("NORMAL", length(lens))
  }
  
  n_reg <- length(lens)
  
  U <- matrix(NA_real_, nrow = T, ncol = 2)
  regime_id <- integer(T)
  regimes <- vector("list", n_reg)
  
  # helper: does draw_fun accept a state argument?
  draw_has_arg <- (length(formals(draw_fun)) >= 1L)
  
  t0 <- 1L
  for (r in seq_len(n_reg)) {
    t1 <- t0 + lens[r] - 1L
    st <- states[r]
    
    dgp_r <- if (draw_has_arg) draw_fun(st) else draw_fun()
    
    U[t0:t1, ] <- rvinecop(lens[r], dgp_r$vc)
    U[t0:t1, ] <- pmin(U[t0:t1, ], 0.9999)             # avoid possible 1
    
    
    u_star_r <- c(
      a0.05b0.05   = solve_u_star(dgp_r$vc, alpha = 0.05,  beta = 0.05),
      a0.05b0.1    = solve_u_star(dgp_r$vc, alpha = 0.10,  beta = 0.05),
      a0.1b0.1     = solve_u_star(dgp_r$vc, alpha = 0.10,  beta = 0.10),
      a0.1b0.05    = solve_u_star(dgp_r$vc, alpha = 0.05,  beta = 0.10),
      a0.05b0.025  = solve_u_star(dgp_r$vc, alpha = 0.025, beta = 0.05),
      a0.025b0.05  = solve_u_star(dgp_r$vc, alpha = 0.05,  beta = 0.025)
    )
    
    regimes[[r]] <- list(
      r = r,
      state   = st,
      t_start = t0,
      t_end   = t1,
      idx     = t0:t1,
      vc      = dgp_r$vc,
      true_bases = dgp_r$true_bases,
      u_star  = u_star_r,
      
      # optional diagnostics (won't break anything)
      lambdaL = if (!is.null(dgp_r$pair_list$e12$lambdaL)) dgp_r$pair_list$e12$lambdaL else NA_real_,
      lambdaU = if (!is.null(dgp_r$pair_list$e12$lambdaU)) dgp_r$pair_list$e12$lambdaU else NA_real_
    )
    
    regime_id[t0:t1] <- r
    t0 <- t1 + 1L
  }
  
  true_base_by_regime <- vapply(regimes, function(g) g$true_bases[1], character(1))
  true_base_t <- true_base_by_regime[regime_id]
  
  list(U = U, regimes = regimes, regime_id = regime_id, true_base_t = true_base_t)
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

# Compare out$CoVaR_tail (n_test x 2 x 4) against truth series
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
      a0.05b0.05  = "covar_true_a05b05",
      a0.05b0.1   = "covar_true_a05b10",
      a0.1b0.1    = "covar_true_a10b10",
      a0.1b0.05   = "covar_true_a10b05",
      
      # NEW
      a0.05b0.025 = "covar_true_a05b0025",
      a0.025b0.05 = "covar_true_a0025b05"
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


# REGIMES PARALLED
run_one_sim <- function(s,
                        n_train, n_test, d,
                        mean_len = 100L,
                        p_extreme = 0.0,
                        tip_k = NA_integer_,
                        cfg=NULL,
                        seed=NULL,
                        out_dir = "simul_results/2d_smc") {

  tryCatch({

  Ttot <- n_train + n_test

  sched <- make_regime_schedule(
    T = Ttot,
    mean_len  = mean_len,   
    p_extreme = p_extreme,   
    min_len = 20L,
    max_len = 2000000L
  )

  piece <- make_piecewise_dgp_d2(
    T = Ttot,
    L_switch = sched,
    draw_fun = draw_vine_d2_regime
  )
  
  dur <- vapply(piece$regimes, function(g) length(g$idx), integer(1))
  st  <- vapply(piece$regimes, function(g) g$state, character(1))

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

  omega <- 1.785714e-05
  alpha <- 0.05
  beta  <- 0.90

  simulate_garch <- function(n, mu_fc, U_transformed,
                             omega = 1.785714e-05, alpha = 0.05, beta = 0.9) {
    h <- matrix(NA, nrow = n, ncol = d)
    r <- matrix(NA, nrow = n, ncol = d)
    h[1, ] <- omega / (1 - alpha - beta)
    r[1, ] <- sqrt(h[1, ]) * U_transformed[1, ]
    for (t in 2:n) {
      h[t, ] <- omega + alpha * r[t - 1, ]^2 + beta * h[t - 1, ]
      r[t, ] <- sqrt(h[t, ]) * U_transformed[t, ]
    }
    list(returns = r, sig = h)
  }

  garch_result <- simulate_garch(Ttot, mu_fc, U_transformed)
  simulated_returns <- garch_result$returns
  sig_fc <- matrix(garch_result$sig, ncol = d)

  data <- list(
    U        = U,
    mu_fc    = mu_fc,
    sig_fc   = sig_fc,
    df_fc    = df_fc,
    shape_fc = shape_fc,
    y_real   = simulated_returns
  )

  # --- Truth u_star_t over time ---
  u_star_t <- matrix(
    NA_real_,
    nrow = Ttot,
    ncol = length(SCEN_COVAR),
    dimnames = list(NULL, SCEN_COVAR)
  )

  cols <- colnames(u_star_t)
  for (g in dgp$regimes) {
    u_star_t[g$idx, cols] <- matrix(g$u_star[cols], nrow  = length(g$idx), ncol  = length(cols), byrow = TRUE)
  }

  alphas_covar <- c(0.025, 0.05, 0.10)
  betas_covar  <- c(0.025, 0.05, 0.10)

  q_sstd <- function(u) qsstd(u, mean = 0, sd = 1, nu = 3, xi = 0.9)
  scale_by_time <- function(mat, v) mat * matrix(v, nrow = nrow(mat), ncol = ncol(mat))

  q_u_t <- apply(u_star_t, 2, q_sstd)

  covar_true_a05b05 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.05b0.05"])
  covar_true_a05b10 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.05b0.1"])
  covar_true_a10b10 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.1b0.1"])
  covar_true_a10b05 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.1b0.05"])
  covar_true_a05b0025 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.05b0.025"])
  covar_true_a0025b05 <- scale_by_time(sqrt(data$sig_fc), q_u_t[, "a0.025b0.05"])

  q_beta <- setNames(q_sstd(betas_covar), paste0("b", betas_covar))

  dgp$truth <- list(
    u_star_t     = u_star_t,
    alphas_covar = alphas_covar,
    betas_covar  = betas_covar,
    var_true_10  = sqrt(data$sig_fc) * q_beta[3],
    var_true_5   = sqrt(data$sig_fc) * q_beta[2],
    var_true_025   = sqrt(data$sig_fc) * q_beta[1],
    covar_true_a05b05 = covar_true_a05b05,
    covar_true_a05b10 = covar_true_a05b10,
    covar_true_a10b10 = covar_true_a10b10,
    covar_true_a10b05 = covar_true_a10b05,
    covar_true_a05b0025 = covar_true_a05b0025,
    covar_true_a0025b05 = covar_true_a0025b05
  )

  #cfg <- modifyList(build_cfg(d = 2), list(M = 500, label = "M500", use_tail_informed_prior = TRUE, tip_k=25))

  # --- Run your method ---
  out <- naive_simul_d2_regimes(data, cfg, dgp)
  #out <- smc_simul_serial(data, cfg, dgp)

  n_oos <- nrow(data$U) - cfg$W_predict
  y_real_oos  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  rp_real_oos <- rowMeans(y_real_oos)

  eval_var <- port_var_backtest_df(out, rp_real_oos, cfg$alphas)
  grid_ab_use <- rbind(
    data.frame(alpha_j=0.10,  alpha_port=0.10),
    data.frame(alpha_j=0.10,  alpha_port=0.05),
    data.frame(alpha_j=0.05,  alpha_port=0.10),
    data.frame(alpha_j=0.05,  alpha_port=0.05),
    data.frame(alpha_j=0.05,  alpha_port=0.025),  # NEW
    data.frame(alpha_j=0.025, alpha_port=0.05)    # NEW
  )
  
  eval_covar <- covar_backtest_grid(
    out         = out,
    y_real_oos  = y_real_oos,
    rp_real_oos = rp_real_oos,
    cfg_alphas  = cfg$alphas,
    grid_ab     = grid_ab_use,
    d           = d
  )

  rmse_mae_from_covar <- covar_rmse_mae_all(out, dgp)
  
  eval_covar_asset <- covar_backtest_grid_asset(
    out         = out,
    y_real_oos  = y_real_oos,
    rp_real_oos = rp_real_oos,
    cfg_alphas  = cfg$alphas,
    grid_ab     = grid_ab_use,
    d           = d
  )

  eval_var$sim   <- s
  eval_covar$sim <- s

  
  result <- list(
    eval_var  = eval_var,
    eval_covar = eval_covar,
    eval_covar_asset = eval_covar_asset,
    log_pred  = out$log_pred,
    QL        = apply(out$port$QL, 2, sum),
    FZL       = apply(out$port$FZL, 2, sum),
    wCRPS     = sum(out$port$wCRPS),
    rmse_mae_from_covar = rmse_mae_from_covar,
    dgp       =dgp,
    
    # record grid identifiers in each file (handy later)
    mean_len  = mean_len,
    p_extreme = p_extreme,
    tip_k     = tip_k,
    sim       = s,
    sched = sched,
    piece = piece,
    data = data,
    dgp = dgp$truth,
    cfg = cfg
    
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  f_final <- file.path(out_dir, sprintf("sim_%03d.rds", s))
  f_tmp   <- paste0(f_final, ".tmp")
  
  saveRDS(result, f_tmp)
  file.rename(f_tmp, f_final)
  
  return(f_final)  # return only the filename (lightweight)
  
  }, error = function(e) {
    # save error so you can inspect failures later
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    f_final <- file.path(out_dir, sprintf("sim_%03d_ERROR.rds", s))
    saveRDS(list(
      sim = s,
      mean_len = mean_len,
      p_extreme = p_extreme,
      tip_k = tip_k,
      error = conditionMessage(e)
    ), f_final)
    return(f_final)
  })
}












