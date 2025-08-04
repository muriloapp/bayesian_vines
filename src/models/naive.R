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


data <- list(
  U      = readRDS("data/PIT.rds")[,1:3],
  mu_fc  = readRDS("data/returns_mean_forecast.rds")[,2:4],# [,-1],  # drop date col
  sig_fc = readRDS("data/returns_vol_forecast.rds")[,2:4],  #[,-1],
  df_fc = readRDS("data/df_fc.rds")[,2:4],#[,-1]
  shape_fc = readRDS("data/shape_fc.rds")[,2:4]#[,-1]
)


U         <- data$U
mu_fc     <- data$mu_fc
sig_fc    <- data$sig_fc
df_fc     <- data$df_fc 
shape_fc     <- data$shape_fc 

cfg <- build_cfg(d = ncol(U))
N <- nrow(U); K <- cfg$K; tickers <- colnames(U); A <- length(cfg$alphas)
n_oos <- N - cfg$W_predict; d <- cfg$d; t_train <- cfg$W_predict


out <- list(
  log_pred    = numeric(n_oos),
  fam_hist  = matrix(NA_integer_, n_oos,K),
  par1_hist = matrix(NA_real_, n_oos,K),
  par2_hist = matrix(NA_real_, n_oos,K),
  risk    =  list(
    dates = integer(n_oos),                                        # 1‑D
    mean  = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    var   = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    ci_lo = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    ci_hi = matrix(NA_real_, n_oos, d, dimnames = list(NULL, tickers)),
    VaR   = array(NA_real_, dim = c(n_oos, d, A),
                  dimnames = list(NULL, tickers, paste0("a", cfg$alphas))),
    ES    = array(NA_real_, dim = c(n_oos, d, A),
                  dimnames = list(NULL, tickers, paste0("a", cfg$alphas)))
  ),                  # will rbind() rows
  port = list(
    dates  = integer(n_oos),
    mean   = numeric(n_oos),
    VaR    = matrix(NA_real_, n_oos, A,
                    dimnames = list(NULL, paste0("a", cfg$alphas))),
    ES     = matrix(NA_real_, n_oos, A,
                    dimnames = list(NULL, paste0("a", cfg$alphas)))
  ) 
)



for (t in seq_len(n_oos)) {
  
  test_idx <- t_train + t
  u_t <- U[(test_idx - t_train):(test_idx - 1), , drop = FALSE]
  if (t == 1)  skel <- make_skeleton_CVM(u_t)
  model <- vinecop(u_t,
                   family_set = c("gaussian", "bb1"),
                   structure   = skel$structure,
                   allow_rotations = FALSE)
  
  # out$fam_hist [t, ] <- t(vapply(particles, `[[`, integer(K),"fam"))
  # out$par1_hist[t, ]  <- t(vapply(particles, `[[`, numeric(K),"th1"))
  # out$par2_hist[t, ]  <- t(vapply(particles, `[[`, numeric(K),"th2"))
  # 
  
  ## --- 1-step-ahead predictive density ---------------------------------
  u_test            <- U[test_idx, , drop = FALSE]
  out$log_pred[t]   <- dvinecop(u_test, model, cores = 7)
  
  ## --- draws & risk metrics -------------------------------------------
  draws <- rvinecop(n=5000, model, qrng = FALSE, cores = 7)
  
  Z_pred <- st_inv_fast(draws, shape_fc[t, ], df_fc[t, ])  
  R_t <- sweep(Z_pred, 2,
               as.numeric(sig_fc[t, ]), `*`) +   # ← aligned
    as.numeric(mu_fc[t, ])
  
  
  rs <- risk_stats_full(            
      R_t,
      cfg$alphas)
    
    out$risk$dates[t]   <- t
    out$risk$mean [t, ] <- rs$mean
    out$risk$var  [t, ] <- rs$var
    out$risk$ci_lo[t, ] <- rs$ci["lo", ]
    out$risk$ci_hi[t, ] <- rs$ci["hi", ]
    out$risk$VaR [t, , ] <- rs$VaR         
    out$risk$ES  [t, , ] <- rs$ES
    
    # ---------- EW-portfolio metrics -------------------
    r_p  <- rowMeans(R_t)                                            
    ps   <- port_stats(r_p, cfg$alphas)                             
    
    out$port$dates[t]   <- t                                     
    out$port$mean [t]   <- ps$mu
    out$port$VaR [t, ]  <- ps$VaR
    out$port$ES  [t, ]  <- ps$ES  
    
    
    
    print(t)
  
}
  
  
saveRDS(out, file = file.path("empirical_results", "out_naive_3.rds"))

  
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










out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/out_naive_3.rds")




