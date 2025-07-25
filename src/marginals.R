
library(data.table)   
library(rugarch)      
library(readxl)     
library(xts)

win_len     <- 756      # rolling window length
refit_every <- 252      # re‑estimate GARCH every 252 obs

spec_norm <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),  # AR(1)
  distribution.model = "norm"                                       # Normal
)

################################################################################
## 1.  Load data & reshape                                                     ##
################################################################################

dt_long <- read_excel('C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/data/data.xlsx')          
dt_long <- as.data.table(dt_long)
dt_long[, Date := as.Date(Date, format = "%d/%m/%Y")]
dt_wide <- dcast(dt_long, Date ~ Ticker, value.var = "Ret")
setorder(dt_wide, Date)                           

ret_xts <- xts(dt_wide[, -1, with = FALSE], order.by = dt_wide$Date)

# Excluding MS due to data issues and excluding one particular row wihtou values for all assets
ret_xts <- ret_xts[, !colnames(ret_xts) %in% "MS"]
ret_xts <- ret_xts[-2728, ]

tickers <- colnames(ret_xts)
K       <- length(tickers)


###############################################################################
## 2.  Helpers: fit once, then forecast 1‑step w/out re‑estimation           ##
###############################################################################
fit_garch_window <- function(x, end_idx) {
  ugarchfit(spec_norm, data = x[(end_idx - win_len + 1):end_idx])
}

forecast_next <- function(par_vec, x_window) {
  spec_fixed <- ugarchspec(
    variance.model     = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model         = list(armaOrder = c(1,0), include.mean = TRUE),
    distribution.model = "norm",
    fixed.pars         = as.list(par_vec)
  )
  fc <- ugarchforecast(spec_fixed, data = x_window, n.ahead = 1)
  c(mu    = as.numeric(fitted(fc)),        # µ_{t+1|t}
    sigma = as.numeric(sigma(fc)))         # σ_{t+1|t}
}

# -- 2c  **NEW**: filtered standardised residual for the *last* obs ---
stdres_last <- function(par_vec, x_window_plus_new) {
  spec_fixed <- ugarchspec(
    variance.model     = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model         = list(armaOrder = c(1,0), include.mean = TRUE),
    distribution.model = "norm",
    fixed.pars         = as.list(par_vec)
  )
  filt <- ugarchfilter(spec_fixed, data = x_window_plus_new)
  tail(residuals(filt, standardize = TRUE), 1)           # numeric(1)
}


################################################################################
## 3.  Initial in‑sample fit  (days 1 … 756)                                   ##
################################################################################
ins_end    <- win_len
garch_fit  <- lapply(tickers, function(sym)
  fit_garch_window(ret_xts[, sym], end_idx = ins_end))
names(garch_fit) <- tickers

## PITs for the *training* window  (size 756 × K)
pit_mat <- do.call(cbind,
                   lapply(garch_fit,
                          function(f) pnorm(residuals(f, standardize = TRUE))))
colnames(pit_mat) <- tickers


################################################################################
## 4.  Pre‑allocate storage for OOS results                                    ##
################################################################################
n_oos       <- nrow(ret_xts) - win_len
date_oos    <- index(ret_xts)[(win_len + 1):(win_len + n_oos)]  # t+1 dates

mu_fc       <- sigma_fc <- actual_ret <- array(
  NA_real_, dim = c(n_oos, K),
  dimnames = list(NULL, tickers)
)

## grow‑as‑we‑go PIT matrix: start with the 756 in‑sample rows
pit_full    <- pit_mat

###############################################################################
## 5.  Rolling‑window forecast loop                                          ##
###############################################################################
for (h in seq_len(n_oos)) {
  t_idx <- ins_end + h - 1          # last obs available = time t
  refit_now <- ((h - 1) %% refit_every) == 0
  u_row <- numeric(K) 
  
  for (k in seq_len(K)) {
    ## 5.1  Re‑estimate (756‑window) only every 252 days
    if (refit_now) {
      cat("Refitting in h equal to", h, "for asset", k, "\n")
      garch_fit[[k]] <- fit_garch_window(ret_xts[, k], t_idx)
  
    }
    
    ## 5.2  Forecast r_{t+1}: mean & vol with *fixed* parameters
    par_k   <- coef(garch_fit[[k]])
    x_win   <- ret_xts[(t_idx - win_len + 1):t_idx, k]
    fc_vals <- forecast_next(par_k, x_win)
    
    mu_fc    [h, k] <- fc_vals["mu"]
    sigma_fc [h, k] <- fc_vals["sigma"]
    
    ## 5.3  Store the realised return r_{t+1}
    r_tp1           <- coredata(ret_xts[t_idx + 1, k])
    actual_ret[h, k] <- r_tp1
    
    # 5d  **NEW** PIT based on filtered residual (in‑sample)
    x_win_plus <- ret_xts[(t_idx - win_len + 2):(t_idx + 1), k]  # 756 obs incl. r_{t+1}
    z_new      <- stdres_last(par_k, x_win_plus)                 # ε/σ at t+1
    u_row[k]       <- pnorm(z_new)                                   # PIT under Normal
    
  }
  # 5e  append one PIT row (xts) with the correct date index
  pit_row_xts <- xts(t(u_row),                   # make it a 1×K matrix
                     order.by = index(ret_xts)[t_idx + 1])  # date = t+1
  colnames(pit_row_xts) <- tickers
  
  pit_full <- rbind(pit_full, pit_row_xts)
}


###############################################################################
## 6.  Save                                                                  ##
###############################################################################

date_col <- data.table(Date = date_oos)

# realised returns  r_t              →  ./data/returns_actual.rds
actual_dt <- cbind(date_col, as.data.table(actual_ret))
saveRDS(actual_dt, file = file.path("data", "returns_actual.rds"))

# forecast means   μ̂_{t|t-1}        →  ./data/returns_mean_forecast.rds
mean_dt   <- cbind(date_col, as.data.table(mu_fc))
saveRDS(mean_dt,   file = file.path("data", "returns_mean_forecast.rds"))

# forecast vols    σ̂_{t|t-1}        →  ./data/returns_vol_forecast.rds
sigma_dt  <- cbind(date_col, as.data.table(sigma_fc))
saveRDS(sigma_dt,  file = file.path("data", "returns_vol_forecast.rds"))


# optional sanity check
print(list(head(actual_dt, 2), head(mean_dt, 2), head(sigma_dt, 2)))

## ------------ persist to disk ----------------------------------------------
saveRDS(pit_full,  file = file.path("data", "PIT.rds"))












