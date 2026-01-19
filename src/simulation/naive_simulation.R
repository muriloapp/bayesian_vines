library(rvinecopulib)
library(VineCopula)
library(here)
library(data.table)

source(here("src/R", "utils.R"))
source(here("src/R", "metrics.R"))


# build_cfg <- function(d,
#                       K = d * (d - 1) / 2) {
#   list(
#     d       = d,
#     K       = K,
#     W_predict    = 756L,                                 
#     alphas     = c(.1, .05, .025, .01),
#     nc       = max(parallel::detectCores()-1, 1)
#   )
# }

# map_families <- function(fams) {
#   mapping <- c(indep = 0, gaussian = 1, t=2, bb1 = 3, bb8 = 6)
#   return(unname(mapping[fams]))
# }

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

# extend_with_gaussian <- function(vinecop_t1, structure_full) {
#   # get dimension from structure (or from T1 length + 1)
#   d <- dim(structure_full)[1]
#   
#   obj <- vinecop_t1                       # keep class 'vinecop'
#   obj$structure <- structure_full         # ensure full structure set
#   
#   # ensure 'pair_copulas' has d-1 trees; keep Tree 1 as-is
#   pcs <- obj$pair_copulas
#   pcs <- c(pcs, vector("list", (d - 1) - length(pcs)))  # pad if needed
#   
#   # fill Trees 2..(d-1) with Gaussian placeholders
#   for (t in 2:(d - 1)) {
#     pcs[[t]] <- replicate(d - t, bicop_dist("gaussian", parameters = 0), simplify = FALSE)
#   }
#   obj$pair_copulas <- pcs
#   obj
# }
# 
# extend_with_t <- function(vinecop_t1, structure_full) {
#   # get dimension from structure (or from T1 length + 1)
#   d <- dim(structure_full)[1]
#   
#   obj <- vinecop_t1                       # keep class 'vinecop'
#   obj$structure <- structure_full         # ensure full structure set
#   
#   # ensure 'pair_copulas' has d-1 trees; keep Tree 1 as-is
#   pcs <- obj$pair_copulas
#   pcs <- c(pcs, vector("list", (d - 1) - length(pcs)))  # pad if needed
#   
#   # fill Trees 2..(d-1) with Gaussian placeholders
#   for (t in 2:(d - 1)) {
#     pcs[[t]] <- replicate(d - t, bicop_dist("t", parameters = c(0,5)), simplify = FALSE)
#   }
#   obj$pair_copulas <- pcs
#   obj
# }



# your name->code map (extend if needed)
name_to_code <- c(
  indep    = 0,
  gaussian = 1,
  t        = 2,
  bb1      = 7,
  bb7      = 9,
  bb1r180  = 17,
  bb7r180  = 19
)

# convert VineCopula BiCopSelect -> rvinecopulib bicop_dist + pieces for hbicop
vc_to_rvl <- function(fit_vc) {
  code <- fit_vc$family
  
  # family + rotation (VineCopula rotation is encoded in the integer code)
  fam_rot <- switch(as.character(code),
                    "0"  = list(fam="indep",    rot=0L),
                    "1"  = list(fam="gaussian", rot=0L),
                    "2"  = list(fam="t",        rot=0L),
                    "7"  = list(fam="bb1",      rot=0L),
                    "17" = list(fam="bb1",      rot=180L),
                    "9"  = list(fam="bb7",      rot=0L),
                    "19" = list(fam="bb7",      rot=180L),
                    stop("Unsupported VineCopula family code: ", code)
  )
  
  # parameters
  pars <- if (fam_rot$fam == "indep") numeric(0)
  else if (fam_rot$fam %in% c("gaussian")) c(fit_vc$par)
  else c(fit_vc$par, fit_vc$par2)  # t, bb1, bb7
  
  list(
    family     = fam_rot$fam,
    rotation   = fam_rot$rot,
    parameters = pars,
    dist       = bicop_dist(fam_rot$fam, fam_rot$rot, pars)
  )
}



#####################################################################################

naive_simul <- function(data, cfg, dgp) {

  true_bases <- dgp$true_bases
  U      <- data$U
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers    <- colnames(U); A <- length(cfg$alphas); t_train <- cfg$W_predict
  
  mu_fc   <- data$mu_fc[(nrow(data$mu_fc)   - n_oos + 1):nrow(data$mu_fc), , drop = FALSE]
  sig_fc  <- data$sig_fc[(nrow(data$sig_fc) - n_oos + 1):nrow(data$sig_fc), , drop = FALSE]
  df_fc   <- data$df_fc[(nrow(data$df_fc)   - n_oos + 1):nrow(data$df_fc), , drop = FALSE]
  shape_fc<- data$shape_fc[(nrow(data$shape_fc) - n_oos + 1):nrow(data$shape_fc), , drop = FALSE]
  y_real  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  
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
    
    if (t == 1)  skel <- make_skeleton_CVM(U[1:t_train, ], trunc_tree = cfg$trunc_tree, structure = dgp$vc$structure)
  
    if (t == 1 || (t %% 50) == 0) {
      allowed_names <- cfg$families_first
      
      # ---- edge (1,2)
      allowed_names_masked <- allowed_names[ (allowed_names) != (true_bases[1]) ]
      allowed_codes <- unname(name_to_code[allowed_names_masked])
      
      fit12_vc  <- BiCopSelect(u1=u_train[,1], u2=u_train[,2],
                               familyset=allowed_codes, selectioncrit="AIC", rotations=FALSE)
      fit12 <- vc_to_rvl(fit12_vc)
      
      # ---- edge (1,3)
      allowed_names_masked <- allowed_names[ (allowed_names) != (true_bases[2]) ]
      allowed_codes <- unname(name_to_code[allowed_names_masked])
      
      fit13_vc  <- BiCopSelect(u1=u_train[,1], u2=u_train[,3],
                               familyset=allowed_codes, selectioncrit="AIC", rotations=FALSE)
      fit13 <- vc_to_rvl(fit13_vc)
      
      # ---- pseudo obs for Tree 2 (C-vine with edges (1,2) and (1,3) => conditional edge uses |1)
      # IMPORTANT: use u_train (not u)
      u2_given_1 <- hbicop(cbind(u_train[,1], u_train[,2]),
                           cond_var=1, family=fit12$family, rotation=fit12$rotation, parameters=fit12$parameters)
      
      u3_given_1 <- hbicop(cbind(u_train[,1], u_train[,3]),
                           cond_var=1, family=fit13$family, rotation=fit13$rotation, parameters=fit13$parameters)
      
      
      # ---- Tree 2 edge selection on (U2|1, U3|1)
      allowed_names_masked <- allowed_names[ (allowed_names) != (true_bases[3]) ]
      allowed_codes <- unname(name_to_code[allowed_names_masked])
      
      fit23_1_vc <- BiCopSelect(u1=u2_given_1, u2=u3_given_1,
                                familyset=allowed_codes, selectioncrit="AIC", rotations=FALSE)
      fit23_1 <- vc_to_rvl(fit23_1_vc)
      
      # ---- assemble final rvinecopulib object (structure unchanged)
      pair_copulas <- list(
        list(fit12$dist, fit13$dist),   # tree 1
        list(fit23_1$dist)              # tree 2
      )
      
      model <- vinecop_dist(pair_copulas = pair_copulas, structure = dgp$vc$structure)
    } 
    
    
    out$fam_hist[t, ] <- name_to_code[unlist(get_all_families(model, trees = NA))]
    out$par1_hist[t, ] <- extract_params(model)$par1
    out$par2_hist[t, ] <- extract_params(model)$par2
    out$rotation_hist[t, ] <- extract_rotations(model)
    
  
    ## Predictive density
    u_test            <- U[test_idx, , drop = FALSE]
    out$log_pred[t]   <- dvinecop(u_test, model, cores = cfg$nc)
    
    ## Draws & risk metrics 
    draws <- rvinecop(n=5000, model, qrng = FALSE, cores = cfg$nc)
    Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])  
    R_t <- sweep(Z_pred, 2, as.numeric(sqrt(sig_fc[t, ])), `*`) + as.numeric(mu_fc[t, ])
    
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
  
    
    out$QL[t, , ]  <- pinball_matrix(t(as.matrix(y_real_t, drop=FALSE)), rs$VaR, cfg$alphas) # Quantile loss per asset & alpha using the VaR you already computed
    out$FZL[t, , ] <- fzl_pzc_matrix(y_real_t, rs$VaR, rs$ES, cfg$alphas) # FZL joint loss for (VaR, ES)
    out$wCRPS[t, ] <- wcrps_gr_matrix(R_t, t(as.matrix(y_real_t, drop=FALSE))) # Weighted CRPS from predictive draws 'R_t' and realization
    
    
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

out$cfg <- cfg
out
}



naive_simul_d2 <- function(data, cfg, dgp) {
  
  true_bases <- dgp$true_bases
  U      <- data$U
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers    <- colnames(U); A <- length(cfg$alphas); t_train <- cfg$W_predict
  
  mu_fc   <- data$mu_fc[(nrow(data$mu_fc)   - n_oos + 1):nrow(data$mu_fc), , drop = FALSE]
  sig_fc  <- data$sig_fc[(nrow(data$sig_fc) - n_oos + 1):nrow(data$sig_fc), , drop = FALSE]
  df_fc   <- data$df_fc[(nrow(data$df_fc)   - n_oos + 1):nrow(data$df_fc), , drop = FALSE]
  shape_fc<- data$shape_fc[(nrow(data$shape_fc) - n_oos + 1):nrow(data$shape_fc), , drop = FALSE]
  y_real  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  
  out <- list(
    log_pred    = numeric(n_oos),
    fam_hist      <- matrix(NA_integer_, n_oos, 1),
    par1_hist     <- matrix(NA_real_,    n_oos, 1),
    par2_hist     <- matrix(NA_real_,    n_oos, 1),
    rotation_hist <- matrix(NA_integer_, n_oos, 1),
    
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
    
    if (t == 1 || (t %% 50) == 0) {
      allowed_names <- cfg$families_first
      
      # ---- edge (1,2)
      allowed_names_masked <- allowed_names[ (allowed_names) != (true_bases[1]) ]
      allowed_codes <- unname(name_to_code[allowed_names_masked])
      
      fit12_vc  <- BiCopSelect(u1=u_train[,1], u2=u_train[,2],
                               familyset=allowed_codes, selectioncrit="AIC", rotations=FALSE)
      fit12 <- vc_to_rvl(fit12_vc)
      

      # ---- assemble final rvinecopulib object (structure unchanged)
      pair_copulas <- list(
        list(fit12$dist)   # tree 1
           # tree 2
      )
      
      model <- vinecop_dist(pair_copulas = pair_copulas, structure = dgp$vc$structure)
      model <- fit12$dist
    } 
    
    out$fam_hist[t]      <- name_to_code[model$family][1]
    out$par1_hist[t]     <- model$parameters[1]
    out$par2_hist[t]     <- ifelse(length(model$parameters) >= 2, model$parameters[2], NA_real_)
    out$rotation_hist[t] <- model$rotation
    
    u_test <- U[test_idx, 1:2, drop = FALSE]
    
    # density on copula scale
    out$log_pred[t] <- log(dbicop(u_test, model))
    
    # predictive draws on copula scale
    draws <- rbicop(5000, model)
    ## Predictive density
    # u_test            <- U[test_idx, , drop = FALSE]
    # out$log_pred[t]   <- dvinecop(u_test, model, cores = cfg$nc)
    # 
    # ## Draws & risk metrics 
    # draws <- rvinecop(n=5000, model, qrng = FALSE, cores = cfg$nc)
    Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])  
    R_t <- sweep(Z_pred, 2, as.numeric(sqrt(sig_fc[t, ])), `*`) + as.numeric(mu_fc[t, ])
    
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
    
    
    out$QL[t, , ]  <- pinball_matrix(t(as.matrix(y_real_t, drop=FALSE)), rs$VaR, cfg$alphas) # Quantile loss per asset & alpha using the VaR you already computed
    out$FZL[t, , ] <- fzl_pzc_matrix(y_real_t, rs$VaR, rs$ES, cfg$alphas) # FZL joint loss for (VaR, ES)
    out$wCRPS[t, ] <- wcrps_gr_matrix(R_t, t(as.matrix(y_real_t, drop=FALSE))) # Weighted CRPS from predictive draws 'R_t' and realization
    
    
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
  
  out
}







###########################################################

naive_simul_d2_regimes <- function(data, cfg, dgp) {
  
  U <- data$U
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d
  n_oos <- N - cfg$W_predict
  t_train <- cfg$W_predict
  A <- length(cfg$alphas)
  
  tickers <- colnames(U)
  if (is.null(tickers)) tickers <- paste0("V", seq_len(d))
  
  # OOS slices (as in your original)
  mu_fc    <- data$mu_fc[(nrow(data$mu_fc)       - n_oos + 1):nrow(data$mu_fc), , drop = FALSE]
  sig_fc   <- data$sig_fc[(nrow(data$sig_fc)     - n_oos + 1):nrow(data$sig_fc), , drop = FALSE]
  df_fc    <- data$df_fc[(nrow(data$df_fc)       - n_oos + 1):nrow(data$df_fc), , drop = FALSE]
  shape_fc <- data$shape_fc[(nrow(data$shape_fc) - n_oos + 1):nrow(data$shape_fc), , drop = FALSE]
  y_real   <- data$y_real[(nrow(data$y_real)     - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  
  # --- require time-varying truth from the piecewise DGP ---
  # dgp$true_base_t must be length N (same index as rows of U)
  if (is.null(dgp$true_base_t) || length(dgp$true_base_t) != N) {
    stop("naive_simul_d2: dgp$true_base_t must exist and have length nrow(data$U).")
  }
  has_regime <- !is.null(dgp$regime_id) && length(dgp$regime_id) == N
  
  # Exclusion rule: "window" (default) excludes all true bases present in the training window;
  #                "current" excludes only the current true base at test_idx.
  #exclude_rule <- if (is.null(cfg$exclude_rule)) "window" else cfg$exclude_rule
  exclude_rule <- "current"
  refit_every  <- if (is.null(cfg$refit_every)) 252L else as.integer(cfg$refit_every)
  
  out <- list(
    log_pred       = numeric(n_oos),
    fam_hist       = matrix(NA_integer_, n_oos, 1),
    par1_hist      = matrix(NA_real_,    n_oos, 1),
    par2_hist      = matrix(NA_real_,    n_oos, 1),
    rotation_hist  = matrix(NA_integer_, n_oos, 1),
    
    # helpful bookkeeping for later analysis
    true_base_hist    = character(n_oos),
    regime_hist       = if (has_regime) integer(n_oos) else NULL,
    blocked_bases_hist= vector("list", n_oos),
    
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
    
    CoVaR_tail = array(
      NA_real_, c(n_oos, d, 4),
      dimnames = list(NULL, tickers, c("a0.05b0.05","a0.05b0.1","a0.1b0.1","a0.1b0.05"))
    )
  )
  
  model <- NULL
  
  for (t in seq_len(n_oos)) {
    test_idx <- t_train + t
    idx_train <- (test_idx - t_train):(test_idx - 1)
    
    u_train  <- U[idx_train, , drop = FALSE]
    y_real_t <- y_real[t, ]
    
    # bookkeeping truth at this time
    out$true_base_hist[t] <- dgp$true_base_t[test_idx]
    if (has_regime) out$regime_hist[t] <- dgp$regime_id[test_idx]
    
    # (re)fit every refit_every steps (same idea as your original 50)
    if (t == 1 || (t %% refit_every) == 0) {
      
      allowed_names <- cfg$families_first
      
      # blocked bases according to rule
      blocked_bases <- if (exclude_rule == "current") {
        dgp$true_base_t[test_idx]
      } else {
        unique(dgp$true_base_t[idx_train])
      }
      blocked_bases <- blocked_bases[!is.na(blocked_bases)]
      out$blocked_bases_hist[[t]] <- blocked_bases
      
      allowed_names_masked <- setdiff(allowed_names, blocked_bases)
      
      # safety: if you accidentally block everything, fall back to unmasked set
      if (length(allowed_names_masked) == 0L) {
        allowed_names_masked <- allowed_names
      }
      
      allowed_codes <- unname(name_to_code[allowed_names_masked])
      
      fit12_vc <- BiCopSelect(
        u1 = u_train[, 1], u2 = u_train[, 2],
        familyset     = allowed_codes,
        selectioncrit = "AIC",
        rotations     = FALSE
      )
      
      fit12 <- vc_to_rvl(fit12_vc)
      model <- fit12$dist
    }
    
    out$fam_hist[t]      <- name_to_code[model$family][1]
    out$par1_hist[t]     <- model$parameters[1]
    out$par2_hist[t]     <- ifelse(length(model$parameters) >= 2, model$parameters[2], NA_real_)
    out$rotation_hist[t] <- model$rotation
    
    u_test <- U[test_idx, 1:2, drop = FALSE]
    
    # density on copula scale
    out$log_pred[t] <- log(dbicop(u_test, model))
    
    # predictive draws on copula scale
    draws <- rbicop(2000, model)
    
    # transform + build returns
    Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])
    R_t <- sweep(Z_pred, 2, as.numeric(sqrt(sig_fc[t, ])), `*`) + as.numeric(mu_fc[t, ])
    
    # Risk metrics
    rs <- risk_stats_full(R_t, cfg$alphas)
    
    out$risk$dates[t]     <- t
    out$risk$VaR [t, , ]  <- rs$VaR
    out$risk$ES  [t, , ]  <- rs$ES
    
    # EW-portfolio metrics
    r_p <- rowMeans(R_t)
    ps  <- port_stats(r_p, cfg$alphas)
    
    out$port$dates[t]    <- t
    out$port$VaR [t, ]   <- ps$VaR
    out$port$ES  [t, ]   <- ps$ES
    
    r_p_real <- mean(as.numeric(y_real_t))
    out$port$QL[t, ]     <- vapply(seq_along(cfg$alphas), function(k) pinball_loss(r_p_real, ps$VaR[k], cfg$alphas[k]), numeric(1))
    out$port$FZL[t, ]    <- vapply(seq_along(cfg$alphas), function(k) fzl_pzc_scalar(r_p_real, ps$VaR[k], ps$ES[k], cfg$alphas[k]), numeric(1))
    out$port$wCRPS[t]    <- wcrps_gr_scalar(r_p, r_p_real)
    
    out$QL[t, , ]        <- pinball_matrix(t(as.matrix(y_real_t, drop = FALSE)), rs$VaR, cfg$alphas)
    out$FZL[t, , ]       <- fzl_pzc_matrix(y_real_t, rs$VaR, rs$ES, cfg$alphas)
    out$wCRPS[t, ]       <- wcrps_gr_matrix(R_t, t(as.matrix(y_real_t, drop = FALSE)))
    
    # CoVaR
    k5  <- which.min(abs(cfg$alphas - 0.05))
    k10 <- which.min(abs(cfg$alphas - 0.10))
    
    VaRj_5  <- rs$VaR[, k5]
    VaRj_10 <- rs$VaR[, k10]
    
    covar5     <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.05, minN = 50)
    covar5b10  <- covar_tail_vec(R_t, r_p, VaRj_5,  port_alpha = 0.10, minN = 50)
    covar10    <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.10, minN = 50)
    covar10b5  <- covar_tail_vec(R_t, r_p, VaRj_10, port_alpha = 0.05, minN = 50)
    
    out$CoVaR_tail[t, , "a0.05b0.05"] <- covar5
    out$CoVaR_tail[t, , "a0.05b0.1"]  <- covar5b10
    out$CoVaR_tail[t, , "a0.1b0.1"]   <- covar10
    out$CoVaR_tail[t, , "a0.1b0.05"]  <- covar10b5
    
    print(t)
  }
  
  out
}





##########################################################

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



port_var_backtest_df <- function(out, rp_real_oos, alphas) {
  # out$port$VaR must be n_oos x length(alphas)
  
  eval_port_var <- lapply(seq_along(alphas), function(k) {
    a <- alphas[k]
    q <- out$port$VaR[, k]
    h <- var_hits(rp_real_oos, q)
    
    list(
      alpha   = a,
      n       = length(h),
      hits    = sum(h),
      rate    = mean(h),
      kupiec  = kupiec_test(h, a),
      chr_ind = christoffersen_ind_test(h),
      chr_cc  = christoffersen_cc_test(h, a)
    )
  })
  
  # flatten to a tidy data.frame (similar style to eval_covar)
  eval_port_df <- do.call(rbind, lapply(eval_port_var, function(x) {
    data.frame(
      alpha     = x$alpha,
      n         = x$n,
      hits      = x$hits,
      rate      = x$rate,
      kupiec_p  = x$kupiec$pval,
      ind_p     = x$chr_ind$pval,
      cc_p      = x$chr_cc$pval,
      row.names = NULL
    )
  }))
  
  eval_port_df
}



# ==== B) Asset VaR backtests (per asset, chosen alphas e.g. 5% & 10%) ====
# alphas_eval <- c(0.10, 0.05, 0.025, 0.01)
# eval_asset_var <- do.call(rbind, lapply(seq_len(d), function(j) {
#   do.call(rbind, lapply(alphas_eval, function(a) {
#     k <- which.min(abs(cfg$alphas - a))
#     qj <- out$risk$VaR[, j, k]
#     hj <- var_hits(y_real_oos[, j], qj)
#     data.frame(
#       asset = tickers[j], alpha = a,
#       n = length(hj), hits = sum(hj), rate = mean(hj),
#       kupiec_p = kupiec_test(hj, a)$pval,
#       ind_p    = christoffersen_ind_test(hj)$pval,
#       cc_p     = christoffersen_cc_test(hj, a)$pval,
#       row.names = NULL
#     )
#   }))
# }))

# 
# 
# alphas_eval <- c(0.10, 0.05)
# fmt_ab <- function(x) ifelse(abs(x - 0.05) < 1e-12, "0.05",
#                              ifelse(abs(x - 0.10) < 1e-12, "0.1", as.character(x)))
# grid_ab <- expand.grid(alpha_j = alphas_eval,
#                        alpha_port = alphas_eval,
#                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
# eval_covar <- do.call(rbind, lapply(seq_len(nrow(grid_ab)), function(i) {
#   a <- grid_ab$alpha_j[i]
#   b <- grid_ab$alpha_port[i]
#   
#   # FIX 1: index in cfg$alphas (not alphas_eval)
#   k <- which.min(abs(cfg$alphas - a))
#   
#   lab <- paste0("a", fmt_ab(a), "b", fmt_ab(b))
#   
#   # FIX 2: VaRj_oos is n_oos x d (no extra [,,1])
#   VaRj_oos  <- out$risk$VaR[, , k, drop = FALSE][, , 1]
#   CoVaR_oos <- out$CoVaR_tail[, , lab]
#   
#   cond_hits <- covar_hits_by_j(rp_real_oos, as.matrix(y_real_oos),
#                                VaRj_oos, CoVaR_oos, alpha = a)
#   
#   do.call(rbind, lapply(seq_len(d), function(j) {
#     hj <- cond_hits$hits[[j]]
#     if (length(hj) == 0) {
#       data.frame(asset=tickers[j], alpha_j=a, alpha_port=b, T_event=0,
#                  rate=NA_real_, kupiec_p=NA_real_, ind_p=NA_real_, cc_p=NA_real_,
#                  row.names=NULL)
#     } else {
#       data.frame(asset=tickers[j], alpha_j=a, alpha_port=b, T_event=length(hj),
#                  rate=mean(hj),
#                  kupiec_p=kupiec_test(hj, b)$pval,
#                  ind_p=christoffersen_ind_test(hj)$pval,
#                  cc_p=christoffersen_cc_test(hj, b)$pval,
#                  row.names=NULL)
#     }
#   }))
# }))





covar_backtest_grid <- function(out,
                                y_real_oos,
                                rp_real_oos,
                                cfg_alphas,
                                alphas_eval = c(0.10, 0.05),
                                d = ncol(y_real_oos),
                                tickers = NULL) {
  
  # label formatter to match dimnames like "a0.05b0.1"
  fmt_ab <- function(x) ifelse(abs(x - 0.05) < 1e-12, "0.05",
                               ifelse(abs(x - 0.10) < 1e-12, "0.1", as.character(x)))
  
  grid_ab <- expand.grid(alpha_j = alphas_eval,
                         alpha_port = alphas_eval,
                         KEEP.OUT.ATTRS = FALSE,
                         stringsAsFactors = FALSE)
  
  eval_covar <- do.call(rbind, lapply(seq_len(nrow(grid_ab)), function(i) {
    a <- grid_ab$alpha_j[i]      # VaR level for conditioning asset j
    b <- grid_ab$alpha_port[i]   # portfolio CoVaR tail level
    
    # IMPORTANT: VaR slice index lives in cfg_alphas (not alphas_eval)
    k <- which.min(abs(cfg_alphas - a))
    
    lab <- paste0("a", fmt_ab(a), "b", fmt_ab(b))
    
    VaRj_oos  <- out$risk$VaR[, , k, drop = FALSE][, , 1]   # n_oos x d
    CoVaR_oos <- out$CoVaR_tail[, , lab]                    # n_oos x d
    
    cond_hits <- covar_hits_by_j(rp_real_oos,
                                 as.matrix(y_real_oos),
                                 VaRj_oos,
                                 CoVaR_oos,
                                 alpha = a)
    
    do.call(rbind, lapply(seq_len(d), function(j) {
      hj <- cond_hits$hits[[j]]
      
      asset_name <- if (!is.null(tickers)) tickers[j] else j
      
      if (length(hj) == 0) {
        data.frame(
          asset      = asset_name,
          alpha_j    = a,
          alpha_port = b,
          T_event    = 0L,
          rate       = NA_real_,
          kupiec_p   = NA_real_,
          ind_p      = NA_real_,
          cc_p       = NA_real_,
          row.names  = NULL
        )
      } else {
        data.frame(
          asset      = asset_name,
          alpha_j    = a,
          alpha_port = b,
          T_event    = length(hj),
          rate       = mean(hj),
          kupiec_p   = kupiec_test(hj, b)$pval,
          ind_p      = christoffersen_ind_test(hj)$pval,
          cc_p       = christoffersen_cc_test(hj, b)$pval,
          row.names  = NULL
        )
      }
    }))
  }))
  
  eval_covar
}

# Tidy summaries ready to print/plot:
#print(eval_port_var)     # list (one per alpha) with stats
#print(eval_asset_var)    # data.frame per asset per alpha
#print(eval_covar)        # data.frame per conditioning asset per alpha







# apply(out$port$FZL, 2, sum)
# apply(out$port$QL, 2, sum)
# sum(out$port$wCRPS)


  


############################################



true_model_simul <- function(data, cfg, dgp) {
  
  true_bases <- dgp$true_bases
  U      <- data$U
  M <- cfg$M; K <- cfg$K; N <- nrow(U); d <- cfg$d; n_oos <- N - cfg$W_predict
  tickers    <- colnames(U); A <- length(cfg$alphas); t_train <- cfg$W_predict
  
  mu_fc   <- data$mu_fc[(nrow(data$mu_fc)   - n_oos + 1):nrow(data$mu_fc), , drop = FALSE]
  sig_fc  <- data$sig_fc[(nrow(data$sig_fc) - n_oos + 1):nrow(data$sig_fc), , drop = FALSE]
  df_fc   <- data$df_fc[(nrow(data$df_fc)   - n_oos + 1):nrow(data$df_fc), , drop = FALSE]
  shape_fc<- data$shape_fc[(nrow(data$shape_fc) - n_oos + 1):nrow(data$shape_fc), , drop = FALSE]
  y_real  <- data$y_real[(nrow(data$y_real) - n_oos + 1):nrow(data$y_real), , drop = FALSE]
  
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
    
    if (t == 1)  skel <- make_skeleton_CVM(U[1:t_train, ], trunc_tree = cfg$trunc_tree, structure = dgp$vc$structure)
    
    model <- dgp$vc 
    
    out$fam_hist[t, ] <- name_to_code[unlist(get_all_families(model, trees = NA))]
    out$par1_hist[t, ] <- extract_params(model)$par1
    out$par2_hist[t, ] <- extract_params(model)$par2
    out$rotation_hist[t, ] <- extract_rotations(model)
    
    
    ## Predictive density
    u_test            <- U[test_idx, , drop = FALSE]
    out$log_pred[t]   <- dvinecop(u_test, model, cores = cfg$nc)
    
    ## Draws & risk metrics 
    draws <- rvinecop(n=5000, model, qrng = FALSE, cores = cfg$nc)
    Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])  
    R_t <- sweep(Z_pred, 2, as.numeric(sqrt(sig_fc[t, ])), `*`) + as.numeric(mu_fc[t, ])
    
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
    
    
    out$QL[t, , ]  <- pinball_matrix(t(as.matrix(y_real_t, drop=FALSE)), rs$VaR, cfg$alphas) # Quantile loss per asset & alpha using the VaR you already computed
    out$FZL[t, , ] <- fzl_pzc_matrix(y_real_t, rs$VaR, rs$ES, cfg$alphas) # FZL joint loss for (VaR, ES)
    out$wCRPS[t, ] <- wcrps_gr_matrix(R_t, t(as.matrix(y_real_t, drop=FALSE))) # Weighted CRPS from predictive draws 'R_t' and realization
    
    
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
  
  out$cfg <- cfg
  out
}
  






















  


