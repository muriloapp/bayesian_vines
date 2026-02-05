
library(data.table)   
library(rugarch)      
library(readxl)     
library(xts)
library(fGarch)



################################################################################
## 1.  Load data & reshape                                                     ##
################################################################################

dt_long <- read_excel('C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/data/data.xlsx')          

library(dplyr)

dt_long <- dt_long %>%
  filter(Ticker %in% c("JPM","BAC", "C","MS","BK", "STT", "GS"))


dt_long <- as.data.table(dt_long)
dt_long[, Date := as.Date(Date, format = "%d/%m/%Y")]
dt_wide <- dcast(dt_long, Date ~ Ticker, value.var = "Ret")
setorder(dt_wide, Date)                           

ret_xts <- xts(dt_wide[, -1], order.by = dt_wide$Date)
ret_xts <- log1p(ret_xts)

# Excluding MS due to data issues and excluding one particular row wihtou values for all assets
#ret_xts <- ret_xts[, !colnames(ret_xts) %in% "MS"]

tickers <- colnames(ret_xts)
K       <- length(tickers)


# ---- SETTINGS ---------------------------------------------------------------

win_len     <- 504
refit_every <- 21


EXPECTED <- c("mu","ar1","omega","alpha1","beta1","skew","shape")

# Log failures here
fail_log <- data.table(asset = character(),
                       h = integer(),
                       date_t1 = as.Date(character()),
                       reason = character(),
                       conv_code = integer())

# last good parameters per asset (filled after initial fit)
last_good <- vector("list", K)
names(last_good) <- tickers


spec_norm <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),  # AR(1)
  distribution.model = "sstd"                                       
)

# ---- small helpers ----------------------------------------------------------
default_pars <- function(x) {
  # conservative defaults consistent with unconditional var
  x <- as.numeric(x); x <- x[is.finite(x)]
  v  <- var(x, na.rm = TRUE); if (!is.finite(v)) v <- 1e-6
  alpha1 <- 0.05; beta1 <- 0.90
  omega  <- max(1e-12, v * (1 - alpha1 - beta1))
  c(mu = mean(x, na.rm = TRUE), ar1 = 0,
    omega = omega, alpha1 = alpha1, beta1 = beta1,
    skew = 1.0, shape = 8.0)
}

safe_coef <- function(fit) {
  # return named vector in EXPECTED order, or NULL
  if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
  if (isTRUE(fit@fit$convergence != 0)) return(NULL)
  p <- (coef(fit))
  p <- p[EXPECTED]
  if (length(p) != length(EXPECTED) || any(!is.finite(p))) return(NULL)
  p
}

spec_with_start <- function(start_pars) {
  ugarchspec(
    variance.model   = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model       = list(armaOrder = c(1,0), include.mean = TRUE),
    distribution.model = "sstd",
    start.pars       = as.list(start_pars)   # warm start (free params)
  )
}

spec_with_fixed <- function(pars) {
  ugarchspec(
    variance.model   = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model       = list(armaOrder = c(1,0), include.mean = TRUE),
    distribution.model = "sstd",
    fixed.pars       = as.list(pars)         # use as-is (no optimization)
  )
}

fit_garch_window_warm <- function(x, end_idx, start_pars = NULL) {
  xw <- as.numeric(x[(end_idx - win_len + 1):end_idx])
  stopifnot(all(is.finite(xw)))
  sp <- if (is.null(start_pars)) spec_norm else spec_with_start(start_pars)
  (
    ugarchfit(
      sp, data = xw,
      solver = "hybrid",
      solver.control = list(trace = 0),
      fit.control = list(scale = TRUE, eval.se = FALSE)  # no Hessian
    )
  )
}

forecast_next_fixed <- function(pars, x_win) {
  pars <- pars[EXPECTED]
  sp   <- spec_with_fixed(pars)
  fc   <- ugarchforecast(sp, data = as.numeric(x_win), n.ahead = 1)
  c(mu    = as.numeric(fitted(fc)),
    sigma = as.numeric(sigma(fc)),
    df    = as.numeric(pars["shape"]),
    xi    = as.numeric(pars["skew"]))
}

stdres_last_fixed <- function(pars, x_win_plus) {
  pars <- pars[EXPECTED]
  sp   <- spec_with_fixed(pars)
  filt <- ugarchfilter(sp, data = as.numeric(x_win_plus))
  tail(residuals(filt, standardize = TRUE), 1)
}

# ---- 1) INITIAL IN-SAMPLE FITS: fill last_good ------------------------------
ins_end <- win_len
garch_fit <- lapply(seq_along(tickers), function(i) {
  fit_garch_window_warm(ret_xts[, i], end_idx = ins_end, start_pars = NULL)
})
names(garch_fit) <- tickers

# extract last_good; if missing, synthesize defaults from window
for (k in seq_len(K)) {
  p <- safe_coef(garch_fit[[k]])
  if (is.null(p)) {
    xw <- ret_xts[(ins_end - win_len + 1):ins_end, k]
    p  <- default_pars(xw)
    fail_log <- rbind(
      fail_log,
      data.table(asset = tickers[k], h = 0,
                 date_t1 = index(ret_xts)[ins_end],
                 reason = "initial-fit/fallback-default", conv_code = NA_integer_)
    )
  }
  last_good[[k]] <- p
}

# ---- 2) TRAINING PITs (optional) -------------------------------------------
pit_mat <- do.call(cbind, lapply(seq_len(K), function(k) {
  z  <- residuals(garch_fit[[k]], standardize = TRUE)
  nu <- last_good[[k]]["shape"]; xi <- last_good[[k]]["skew"]
  psstd(z, mean = 0, sd = 1, nu = as.numeric(nu), xi = as.numeric(xi))
}))
colnames(pit_mat) <- tickers
pit_full <- xts(pit_mat, order.by = index(ret_xts)[1:win_len])

# ---- 3) ROLLING LOOP with warm start + logging + fallback -------------------
n_oos    <- nrow(ret_xts) - win_len
date_oos <- index(ret_xts)[(win_len + 1):(win_len + n_oos)]

df_fc     <- array(NA_real_, dim = c(n_oos, K), dimnames = list(NULL, tickers))
shape_fc  <- array(NA_real_, dim = c(n_oos, K), dimnames = list(NULL, tickers))
mu_fc     <- sigma_fc <- actual_ret <- array(NA_real_, dim = c(n_oos, K),
                                             dimnames = list(NULL, tickers))

for (h in seq_len(n_oos)) {
  t_idx     <- ins_end + h - 1       # last in-sample index (time t)
  refit_now <- ((h - 1) %% refit_every) == 0
  u_row     <- numeric(K)
  
  for (k in seq_len(K)) {
    # ---- 3.1 Refit (warm start) occasionally
    if (refit_now) {
      cat("Refitting in h equal to", h, "for asset", k, "\n")
      start_pars <- last_good[[k]]
      fit_k <- try(
        fit_garch_window_warm(ret_xts[, k], end_idx = t_idx, start_pars = start_pars),
        
      )
      
      p_new <- safe_coef(fit_k)
      if (is.null(p_new)) {
        # log failure and keep previous params
        fail_log <- rbind(
          fail_log,
          data.table(asset = tickers[k],
                     h = h,
                     date_t1 = index(ret_xts)[t_idx + 1],
                     reason = if (inherits(fit_k, "try-error")) "try-error" else "non-converged/invalid",
                     conv_code = if (inherits(fit_k, "try-error")) NA_integer_ else fit_k@fit$convergence)
        )
      } else {
        last_good[[k]] <- p_new
      }
    }
    
    # ---- 3.2 Forecast at t+1 using last_good (fixed.pars)
    par_k  <- last_good[[k]]
    x_win  <- ret_xts[(t_idx - win_len + 1):t_idx, k]
    fc_val <- forecast_next_fixed(par_k, x_win)
    
    mu_fc   [h, k] <- fc_val["mu"] # this is the forecast for this date h
    sigma_fc[h, k] <- fc_val["sigma"]
    df_fc   [h, k] <- fc_val["df"]
    shape_fc[h, k] <- fc_val["xi"]
    
    # realized r_{t+1}
    actual_ret[h, k] <- coredata(ret_xts[t_idx + 1, k])
    
    # ---- 3.3 PIT at t+1 using filtered std. residual under fixed params
    x_win_plus <- ret_xts[(t_idx - win_len + 2):(t_idx + 1), k]
    z_new      <- stdres_last_fixed(par_k, x_win_plus)
    u_row[k]   <- psstd(z_new, mean = 0, sd = 1,
                        nu = as.numeric(df_fc[h, k]),
                        xi = as.numeric(shape_fc[h, k]))
  }
  
  pit_row_xts <- xts(t(u_row), order.by = index(ret_xts)[t_idx + 1])
  colnames(pit_row_xts) <- tickers
  pit_full <- rbind(pit_full, pit_row_xts)
  
  print(h)
}

# ---- 4) OPTIONAL: inspect/save failures ------------------------------------
if (nrow(fail_log)) {
  print(fail_log)
  # saveRDS(fail_log, file = file.path("data", "garch_refit_failures.rds"))
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

df_fc  <- cbind(date_col, as.data.table(df_fc))
saveRDS(df_fc,  file = file.path("data", "df_fc.rds"))

shape_fc  <- cbind(date_col, as.data.table(shape_fc))
saveRDS(shape_fc,  file = file.path("data", "shape_fc.rds"))


## ------------ persist to disk ----------------------------------------------
saveRDS(pit_full,  file = file.path("data", "PIT.rds"))


saveRDS(index(ret_xts), file = file.path("data", "dates.rds"))












































###########################################################################
win_len     <- 504      # rolling window length
refit_every <- 21      # re‑estimate GARCH every 252 obs





###############################################################################
## 2.  Helpers: fit once, then forecast 1‑step w/out re‑estimation           ##
###############################################################################
fit_garch_window <- function(x, end_idx) {
  ugarchfit(spec_norm, data = x[(end_idx - win_len + 1):end_idx])
}

forecast_next <- function(pars, x_win) {
  spec_fix <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1,1)),
    mean.model    =list(armaOrder=c(1,0), include.mean=TRUE),
    distribution.model="sstd",
    fixed.pars = as.list(pars))
  fc <- ugarchforecast(spec_fix, data = x_win, n.ahead = 1)
  c(mu    = as.numeric(fitted(fc)),
    sigma = as.numeric(sigma(fc)),
    df    = as.numeric(pars["shape"]),
    xi    = as.numeric(pars["skew"]))     # save ν  (df)
}

stdres_last <- function(pars, x_win_plus) {
  spec_fix <- ugarchspec(
    variance.model=list(model="sGARCH", garchOrder=c(1,1)),
    mean.model    =list(armaOrder=c(1,0), include.mean=TRUE),
    distribution.model="sstd",
    fixed.pars = as.list(pars))
  filt <- ugarchfilter(spec_fix, data = x_win_plus)
  tail(residuals(filt, standardize = TRUE), 1)
}




################################################################################
## 3.  Initial in‑sample fit  (days 1 … 756)                                   ##
################################################################################
ins_end    <- win_len
garch_fit  <- lapply(tickers, function(sym)
  fit_garch_window(ret_xts[, sym], end_idx = ins_end))
names(garch_fit) <- tickers

## --- Student-t PITs for the 756-day training window ------------
pit_mat <- do.call(cbind, lapply(garch_fit, function(f) {
  z  <- residuals(f, standardize = TRUE)          # ε_t / σ_t
  ν  <- coef(f)["shape"]                          # df parameter
  xi <- coef(f)["skew"]
  pst(z, xi, df=ν)                        # PIT under t_ν
}))
colnames(pit_mat) <- tickers


################################################################################
## 4.  Pre‑allocate storage for OOS results                                    ##
################################################################################
n_oos       <- nrow(ret_xts) - win_len
date_oos    <- index(ret_xts)[(win_len + 1):(win_len + n_oos)]  # t+1 dates

df_fc <- array(NA_real_, dim = c(n_oos, K), dimnames = list(NULL, tickers))
shape_fc <- array(NA_real_, dim = c(n_oos, K), dimnames = list(NULL, tickers))
mu_fc <- sigma_fc <- actual_ret <- array(NA_real_, dim = c(n_oos, K), dimnames = list(NULL, tickers) )

## grow‑as‑we‑go PIT matrix: start with the 756 in‑sample rows
pit_full    <- pit_mat
pit_full <- xts(pit_full, order.by = index(ret_xts)[1:win_len]) 

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
    df_fc   [h, k]    <- fc_vals["df"]
    shape_fc   [h, k]    <- fc_vals["xi"]
    
    ## 5.3  Store the realised return r_{t+1}
    r_tp1           <- coredata(ret_xts[t_idx + 1, k])
    actual_ret[h, k] <- r_tp1
    
    # 5d  **NEW** PIT based on filtered residual (in‑sample)
    x_win_plus <- ret_xts[(t_idx - win_len + 2):(t_idx + 1), k]  # 756 obs incl. r_{t+1}
    z_new      <- stdres_last(par_k, x_win_plus)                 # ε/σ at t+1
    ν      <- df_fc[h, k]                       # just saved above
    shape  <- shape_fc[h, k] 
    #u_row[k] <- pt(z_new, df = ν)                              
    u_row[k] <- pst(z_new, shape=shape, df = ν)     
    
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

df_fc  <- cbind(date_col, as.data.table(df_fc))
saveRDS(df_fc,  file = file.path("data", "df_fc.rds"))

shape_fc  <- cbind(date_col, as.data.table(shape_fc))
saveRDS(shape_fc,  file = file.path("data", "shape_fc.rds"))


## ------------ persist to disk ----------------------------------------------
saveRDS(pit_full,  file = file.path("data", "PIT.rds"))


saveRDS(index(ret_xts), file = file.path("data", "dates.rds"))









