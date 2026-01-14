
library(here)
library(fGarch)
source(here("src/R", "config.R"))   

n_sim   <- 30
d       <- 3
n_train <- 756
n_test  <- 3000
n       <- n_train+n_test

fam_names <- c("bb1", "bb1r180", "bb7", "bb7r180", "t")

cvine_struct <- function(d) cvine_structure(order = 1:d)  # C-vine 1-2-3

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



test_loglik <- function(model, Utest) sum(na.omit(log(dvinecop(Utest, model))))


var_list <- vector("list", n_sim)
covar_list <- vector("list", n_sim)
logpred_list <- vector("list", n_sim) 
QL_list <- vector("list", n_sim)
FZL_list <- vector("list", n_sim)
wCRPS_list <- vector("list", n_sim)


for (s in seq_len(n_sim)) {
  set.seed(1111 + s)
  dgp <- draw_vine_d3()
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
  
  
  cfg <- modifyList(build_cfg(d = 3), list(M = 1000, label = "M1000"))
  
  # # # NAIVE
  out <- tryCatch(
    naive_simul(data, cfg, dgp),
    error = function(e) NULL
  )
  if (is.null(out)) next
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

  # SMC
  # out <- smc_simul(data, cfg, dgp)
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

  
  out <- list(
    var_list     = var_list,
    covar_list   = covar_list,
    logpred_list = logpred_list,
    QL_list      = QL_list,
    FZL_list     = FZL_list,
    wCRPS_list   = wCRPS_list
  )
  out_dir <- "simul_results/new"
  out_file <- file.path(out_dir, sprintf("results_naive_s%03d.rds", s))
  
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
saveRDS(var_list, "simul_results/new/var_list_smc.rds")
saveRDS(covar_list, "simul_results/new/covar_list_smc.rds")
  
  

var_all  <- rbindlist(var_list,  fill = TRUE)
covar_all <- rbindlist(covar_list, fill = TRUE)

scores <- do.call(rbind, records)
print(scores)
# mean(scores$ll_mle); if ("ll_smc" %in% names(scores)) mean(scores$ll_smc)




res <- list()

for (i in 1:30) {
  df <- covar_list[[i]]
  if (is.null(df)) next
  
  sub <- df[df$asset %in% c(1,2,3) &
              abs(df$alpha_j - 0.1) < 1e-12 &
              abs(df$alpha_port - 0.1) < 1e-12, , drop = FALSE]
  
  if (nrow(sub) == 0) next
  sub$sim <- i
  res[[length(res) + 1]] <- sub
}

res_01_01 <- do.call(rbind, res)
res_01_01

res_01_01[res_01_01$asset==1,,drop=FALSE]
mean(res_01_01[res_01_01$asset==3,,drop=FALSE]$rate)



res_01_01 <- do.call(rbind, lapply(seq_along(covar_list), function(i) {
  df <- covar_list[[i]]
  if (is.null(df)) return(NULL)
  
  
  
  if (nrow(sub) == 0) return(NULL)
  
  sub$sim <- i
  sub
}))

res_01_01







out_dir <- "simul_results/new"
s <- 22  # set this

res_list <- vector("list", s)
x <- rep(0, 22)

for (s in seq_len(s)) {
  out_file <- file.path(out_dir, sprintf("results_naive_s%03d.rds", s))
  res <- readRDS(out_file)
  x[s] <- res$QL_list

}



out_file <- file.path(out_dir, sprintf("results_naive_s%03d.rds", s))
res <- readRDS(out_file)

xx <- rep(0, s)
for (s in seq_len(s)) {
  xx[s] <- sum((res$logpred_list[[s]]))
}
mean(xx)


out_file <- file.path(out_dir, sprintf("results_smc_s%03d.rds", s))
res <- readRDS(out_file)

xx_smc <- rep(0, s)
for (s in seq_len(s)) {
  xx_smc[s] <- sum(na.omit(res$logpred_list[[s]]))
}

mean(xx_smc)

