
library(here)
source(here("src/R", "config.R"))  


#### VAR FORECASTING
n_assets <- 1:3

out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/standard_tip.rds")

data <- list(
  U      = readRDS("data/PIT.rds")[,n_assets],
  mu_fc  = readRDS("data/returns_mean_forecast.rds")[,n_assets+1],# [,-1],  # drop date col
  sig_fc = readRDS("data/returns_vol_forecast.rds")[,n_assets+1],  #[,-1],
  df_fc = readRDS("data/df_fc.rds")[,n_assets+1],#[,-1]
  shape_fc = readRDS("data/shape_fc.rds")[,n_assets+1]#[,-1]
)

for (i in 1:5000){
vals <- out$fam_hist[,i , 3]
prop <- prop.table(table(vals))
print(prop)
print(i)
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

out <- append(out, cfg)


y_real_oos = readRDS("data/returns_actual.rds")[,n_assets+1, with=FALSE]
#y_real_oos = readRDS("data/returns_actual.rds")[,1:4, with=FALSE]

y_real_oos <- y_real_oos[(.N - nrow(out$port$VaR) + 1):.N]
rp_real_oos <- rowMeans(y_real_oos)


alphas_eval <- c(0.10, 0.05, 0.025, 0.01)
# ==== A) Portfolio VaR backtests for all alphas ====
eval_port_var <- lapply(seq_along(alphas_eval), function(k) {
  a <- alphas_eval[k]
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
d <- length(n_assets)
tickers <- colnames(data$U)

eval_asset_var <- do.call(rbind, lapply(seq_len(d), function(j) {
  do.call(rbind, lapply(alphas_eval, function(a) {
    k <- which.min(abs(alphas_eval - a))
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
  k   <- which.min(abs(alphas_eval - a))
  
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
  VaRj_oos  <- out$risk$VaR[, , k, drop = FALSE][, , 1]
  CoVaR_oos <- out$CoVaR_tail[, , lab]  # n_oos × d
  
  # hits when r_p ≤ CoVaR(j; a,b) conditional on asset j being in VaR event at level a
  cond_hits <- covar_hits_by_j(rp_real_oos, y_real_oos, VaRj_oos, CoVaR_oos, alpha = a)
  
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





a = 0.1
j = 5

k <- which.min(abs(alphas_eval - a))
qj <- out$risk$VaR[, j, k]
hj <- var_hits(y_real_oos[, j, with=FALSE], qj)

VaRj_oos   <- out$risk$VaR[, , k, drop = FALSE][, , 1]   # n_oos×d
#lab <- if (a == 0.05) "a0.05" else "a0.10"
lab = "a0.1b0.1"
CoVaR_oos  <- out$CoVaR_tail[, , lab]                    # n_oos×d
cond_hits  <- covar_hits_by_j(rp_real_oos, as.matrix(y_real_oos), VaRj_oos, CoVaR_oos, alpha = a)
co_hj <- cond_hits$hits[[j]]



y_real = readRDS("data/returns_actual.rds")[,1:4, with=FALSE]
#y_real_oos = readRDS("data/returns_actual.rds")[,1:4, with=FALSE]

y_real <- y_real[(.N - nrow(out$port$VaR) + 1):.N]
rp_real_oos <- rowMeans(y_real_oos)

df <- cbind(y_real, hj)

length(co_hj) == sum(df$hj)
df$co_hj <- NA              # create empty column

df$co_hj[df$hj == 1] <- co_hj

df <- df[Date >= as.Date("2006-01-01")]
#df <- df[Date <= as.Date("2010-01-01")]


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









































































