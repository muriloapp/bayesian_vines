

port_stats <- function(r_vec, alphas = c(.05, .025)) {
  mu  <- mean(r_vec)
  
  # downside risk
  losses <- r_vec                  # “loss = −return”
  VaR <- vapply(alphas, function(a)  quantile(losses, a), numeric(1))
  ES  <- vapply(seq_along(alphas), function(k) {
    thr <- VaR[k]
    mean(losses[losses <= thr])
  }, numeric(1))
  
  list(mu = mu, VaR = VaR, ES = ES)
}


risk_stats_full <- function(R_t,           # L × d matrix  (future PIT‑draws ⇒ returns)
                            #mu_fc,            # length‑d numeric
                            #sig_fc,           # length‑d numeric
                            alphas = c(.05, .025))
{
  #mu_fc  <- as.numeric(mu_fc)
  #sig_fc <- as.numeric(sig_fc)
  
  ## 1. back‑transform to returns  ---------------------------------------------
  #R_t <- sweep(R_pred, 2, sig_fc, `*`) + rep(mu_fc, each = nrow(R_pred))
  
  mu_hat  <- colMeans(R_t)
  var_hat <- apply(R_t, 2, var)
  CI_lo   <- apply(R_t, 2, quantile, 0.025)
  CI_hi   <- apply(R_t, 2, quantile, 0.975)
  
  ## 2.  VaR & ES  -------------------------------------------------------------
  losses <- R_t
  VaR <- sapply(alphas, function(a)
    apply(losses, 2, quantile, probs = a))  
  ES  <- sapply(seq_along(alphas), function(k) {       # loop over α’s
    vapply(seq_len(ncol(losses)), function(j) {  # loop over assets
      thr <- VaR[j, k]                           # ← correct index
      mean(losses[losses[, j] <= thr, j])        # tail average
    }, numeric(1))
  })
  
  ## 3.  Return as a clean list  ----------------------------------------------
  list(
    mean = mu_hat,                       # length‑d
    var  = var_hat,                      # length‑d
    ci   = rbind(lo = CI_lo, hi = CI_hi),  # 2 × d
    VaR  = VaR,                          # d × length(alphas)
    ES   = ES                            # d × length(alphas)
  )
}