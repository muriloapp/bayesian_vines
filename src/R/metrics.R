

port_stats <- function(r_vec, alphas) {

  # downside risk
  losses <- r_vec                  # “loss = −return”
  VaR <- vapply(alphas, function(a)  quantile(losses, a), numeric(1))
  ES  <- vapply(seq_along(alphas), function(k) {
    thr <- VaR[k]
    mean(losses[losses <= thr])
  }, numeric(1))
  
  list(VaR = VaR, ES = ES)
}


risk_stats_full <- function(R_t,           # L × d matrix  (future PIT‑draws ⇒ returns)
                            #mu_fc,            # length‑d numeric
                            #sig_fc,           # length‑d numeric
                            alphas)
{

  
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
    VaR  = VaR,                          # d × length(alphas)
    ES   = ES                            # d × length(alphas)
  )
}