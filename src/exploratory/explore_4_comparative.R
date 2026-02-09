
# ---- setup -------------------------------------------------------------------
library(data.table)
library(here)

source(here("src/R", "config.R"))   # must provide import_data()

# ---- tests -------------------------------------------------------------------
var_hits <- function(y, q) as.numeric(y <= q)

kupiec_test <- function(h, alpha) {
  n  <- length(h)
  n1 <- sum(h)
  n0 <- n - n1
  # guard against 0/1 proportions
  pi_hat <- if (n > 0) n1 / n else NA_real_
  if (is.na(pi_hat) || pi_hat %in% c(0,1)) {
    return(list(stat = NA_real_, pval = NA_real_, rate = pi_hat))
  }
  lr_uc <- -2 * ( n0*log1p(-alpha) + n1*log(alpha) -
                    ( n0*log1p(-pi_hat) + n1*log(pi_hat) ) )
  pval  <- 1 - pchisq(lr_uc, df = 1)
  list(stat = lr_uc, pval = pval, rate = pi_hat)
}

christoffersen_ind_test <- function(h) {
  if (length(h) < 2) return(list(stat = NA_real_, pval = NA_real_))
  hlag  <- c(NA, h[-length(h)])
  valid <- !is.na(hlag)
  h     <- h[valid]; hlag <- hlag[valid]
  n00 <- sum(hlag == 0 & h == 0)
  n01 <- sum(hlag == 0 & h == 1)
  n10 <- sum(hlag == 1 & h == 0)
  n11 <- sum(hlag == 1 & h == 1)
  n0  <- n00 + n01; n1 <- n10 + n11
  
  # add tiny eps to avoid log(0)
  eps <- 1e-12
  pi   <- (n01 + n11) / (n0 + n1 + eps)
  pi0  <- if (n0 > 0) n01 / (n0 + eps) else 0
  pi1  <- if (n1 > 0) n11 / (n1 + eps) else 0
  
  ll_ind <- (n00*log1p(-pi0+eps) + n01*log(pi0+eps) +
               n10*log1p(-pi1+eps) + n11*log(pi1+eps))
  ll_iid <- ((n0 + n1 - (n01 + n11))*log1p(-pi+eps) + (n01 + n11)*log(pi+eps))
  lr_ind <- -2 * (ll_iid - ll_ind)
  pval   <- 1 - pchisq(lr_ind, df = 1)
  list(stat = lr_ind, pval = pval)
}

christoffersen_cc_test <- function(h, alpha) {
  uc  <- kupiec_test(h, alpha)
  ind <- christoffersen_ind_test(h)
  if (any(is.na(c(uc$stat, ind$stat)))) {
    return(list(stat = NA_real_, pval = NA_real_))
  }
  stat <- uc$stat + ind$stat
  list(stat = stat, pval = 1 - pchisq(stat, df = 2))
}

# ---- helpers -----------------------------------------------------------------
# find alpha index in a numeric vector alphas for a target a
.which_alpha <- function(alphas, a) {
  i <- which.min(abs(alphas - a))
  if (length(i) == 0L || !is.finite(i)) stop("Cannot match alpha.")
  i
}
# label maker for CoVaR dimension names like "a0.05b0.1"
#.fmt_ab <- function(x) sub("0\\.(\\d\\d?)$", "0.\\1", formatC(x, format = "fg", digits = 3, drop0trailing = FALSE))
#.covar_lab <- function(a, b) paste0("a", .fmt_ab(a), "b", .fmt_ab(b))

# Build conditional hit series for CoVaR^{p|j}
covar_hits_by_j <- function(r_p_real, y_real_oos, VaRj_oos, CoVaR_oos) {
  d <- ncol(y_real_oos)
  hits_list <- vector("list", d)
  n_list    <- integer(d)
  for (j in seq_len(d)) {
    mask <- y_real_oos[[j]] <= VaRj_oos[, j]
    n_list[j] <- sum(mask)
    hits_list[[j]] <- if (n_list[j] > 0) as.numeric(r_p_real[mask] <= CoVaR_oos[mask, j]) else numeric(0)
  }
  list(hits = hits_list, n = n_list)
}

# Align data$y_real to out-of-sample arrays in `out` and filter by Date range
# Align data$y_real to out-of-sample arrays in `out` and filter by Date range
align_and_filter <- function(out, data,
                             date_from = as.Date("1900-01-01"),
                             date_to   = as.Date("2100-01-01"),
                             keep_cols = NULL) {
  stopifnot("y_real" %in% names(data))
  y_dt <- as.data.table(data$y_real)
  if (!"Date" %in% names(y_dt)) stop("data$y_real must have a Date column.")
  setkey(y_dt, Date)
  
  n_oos <- nrow(out$port$VaR)
  if (!is.finite(n_oos) || n_oos <= 0) stop("out$port$VaR must be n_oos x K.")
  
  # --- exact OOS block used by the model (by *position*) ---
  y_oos <- y_dt[(.N - n_oos + 1L):.N]
  stopifnot(nrow(y_oos) == n_oos)
  
  # Optional asset subset (keep order)
  if (!is.null(keep_cols)) {
    cols <- c("Date", keep_cols)
    missing <- setdiff(cols, names(y_oos))
    if (length(missing)) stop("Missing columns in y_real: ", paste(missing, collapse = ", "))
    y_oos <- y_oos[, ..cols]
  }
  
  # --- build indices of the OOS block that pass the date filter ---
  keep_idx <- which(y_oos$Date >= date_from & y_oos$Date <= date_to)
  
  # Drop Date for matrices; keep tickers and the original OOS dates for indexing
  tickers <- setdiff(colnames(y_oos), "Date")
  Y_full  <- as.matrix(y_oos[, ..tickers])   # n_oos x d
  
  list(
    dates_oos  = y_oos$Date,        # length n_oos
    keep_idx   = keep_idx,          # positions within OOS block
    Y_keep     = Y_full[keep_idx, , drop = FALSE],  # filtered Y
    tickers    = tickers
  )
}


# ---- single-model evaluator ---------------------------------------------------
evaluate_model <- function(out, data,
                           date_from = as.Date("2003-01-01"),
                           date_to   = as.Date("2100-01-01"),
                           assets    = NULL,
                           alphas_var   = c(0.10, 0.05, 0.025, 0.01),
                           alphas_covar = c(0.10, 0.05)) {
  
  # 1) Align & filter by date
  aligned  <- align_and_filter(out, data, date_from, date_to, keep_cols = assets)
  Y        <- aligned$Y_keep                     # (n_keep x d) filtered by dates
  keep_idx <- aligned$keep_idx                   # positions within the OOS block
  tickers  <- aligned$tickers
  d        <- ncol(Y)
  rp       <- rowMeans(Y)
  
  # 2) Alpha grid for VaR arrays
  alphas_out <- if (!is.null(out$cfg$alphas)) as.numeric(out$cfg$alphas)
  else stop("out$cfg$alphas missing or unknown.")
  
  # 3) Subset model arrays by the *same indices* (no tail())
  # portVaR: (n_oos x K) -> (n_keep x K)
  portVaR <- out$port$VaR[keep_idx, , drop = FALSE]
  
  # riskVaR: (n_oos x d x K) -> (n_keep x d x K)
  riskVaR <- out$risk$VaR[keep_idx, , , drop = FALSE]
  
  # CoVaR: (n_oos x d x L) -> (n_keep x d x L)
  CoVaR   <- out$CoVaR_tail[keep_idx, , , drop = FALSE]
  
  # (optional) ensure asset order matches tickers if dimnames exist
  dn_risk <- dimnames(riskVaR)
  if (!is.null(dn_risk) && length(dn_risk) >= 2 && !is.null(dn_risk[[2]])) {
    pos <- match(tickers, dn_risk[[2]])
    if (anyNA(pos)) stop("Some tickers not found in out$risk$VaR dimnames.")
    riskVaR <- riskVaR[, pos, , drop = FALSE]
    # also reorder CoVaR's asset dimension
    dn_covar <- dimnames(CoVaR)
    if (!is.null(dn_covar) && length(dn_covar) >= 2 && !is.null(dn_covar[[2]])) {
      CoVaR <- CoVaR[, pos, , drop = FALSE]
    }
  }
  
  # 4) Portfolio VaR backtests
  dt_port <- rbindlist(lapply(alphas_var, function(a) {
    k <- .which_alpha(alphas_out, a)
    q <- as.numeric(portVaR[, k])
    h <- var_hits(rp, q)
    data.table(
      alpha   = a,
      n       = length(h),
      hits    = sum(h),
      rate    = mean(h),
      kupiec  = kupiec_test(h, a)$pval,
      ind     = christoffersen_ind_test(h)$pval,
      cc      = christoffersen_cc_test(h, a)$pval
    )
  }))
  
  # 5) Asset VaR backtests
  dt_asset <- rbindlist(lapply(seq_len(d), function(j) {
    rbindlist(lapply(alphas_var, function(a) {
      k  <- .which_alpha(alphas_out, a)
      qj <- as.numeric(riskVaR[, j, k])
      hj <- var_hits(Y[, j], qj)
      data.table(
        asset = tickers[j], alpha = a,
        n = length(hj), hits = sum(hj), rate = mean(hj),
        kupiec = kupiec_test(hj, a)$pval,
        ind    = christoffersen_ind_test(hj)$pval,
        cc     = christoffersen_cc_test(hj, a)$pval
      )
    }))
  }))
  
  # 6) CoVaR backtests (α_j x α_port)
  grid <- CJ(alpha_j = alphas_covar, alpha_port = alphas_covar, sorted = TRUE)
  # available CoVaR labels in the 3rd dim:
  labs_avail <- if (!is.null(dimnames(CoVaR)[[3]])) dimnames(CoVaR)[[3]] else NULL
  
  dt_covar <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    a <- grid$alpha_j[i]
    b <- grid$alpha_port[i]
    k <- .which_alpha(alphas_out, a)
    
    # simple label (and a robust fallback if names are like "a0.10b0.10")
    lab_simple <- paste0("a", a, "b", b)
    lab <- lab_simple
    if (!is.null(labs_avail) && !(lab_simple %in% labs_avail)) {
      lab_alt <- paste0("a", formatC(a, format = "f", digits = 2),
                        "b", formatC(b, format = "f", digits = 2))
      if (lab_alt %in% labs_avail) lab <- lab_alt
      # else: keep lab_simple; will error clearly if not found
    }
    
    VaRj  <- matrix(riskVaR[, , k], ncol = d); colnames(VaRj) <- tickers
    CoVab <- CoVaR[, , lab, drop = FALSE][, , 1]; colnames(CoVab) <- tickers
    
    cond <- covar_hits_by_j(rp, as.data.table(Y), VaRj, CoVab)
    
    rbindlist(lapply(seq_len(d), function(j) {
      hj <- cond$hits[[j]]
      if (length(hj) == 0) {
        data.table(asset = tickers[j], alpha_j = a, alpha_port = b,
                   T_event = 0L, rate = NA_real_, kupiec = NA_real_, ind = NA_real_, cc = NA_real_)
      } else {
        data.table(asset = tickers[j], alpha_j = a, alpha_port = b,
                   T_event = length(hj), rate = mean(hj),
                   kupiec = kupiec_test(hj, b)$pval,
                   ind    = christoffersen_ind_test(hj)$pval,
                   cc     = christoffersen_cc_test(hj, b)$pval)
      }
    }))
  }))
  
  list(port_var = dt_port[], asset_var = dt_asset[], covar = dt_covar[])
}


# ---- multi-model wrapper ------------------------------------------------------
compare_models <- function(out_list_named,
                           data,
                           date_from = as.Date("2003-01-01"),
                           date_to   = as.Date("2100-01-01"),
                           assets    = NULL,
                           alphas_var   = c(0.10, 0.05, 0.025, 0.01),
                           alphas_covar = c(0.10, 0.05)) {
  stopifnot(!is.null(names(out_list_named)), all(nzchar(names(out_list_named))))
  res <- lapply(names(out_list_named), function(m) {
    out <- out_list_named[[m]]
    ans <- evaluate_model(out, data, date_from, date_to, assets, alphas_var, alphas_covar)
    list(
      port_var  = cbind(model = m, ans$port_var),
      asset_var = cbind(model = m, ans$asset_var),
      covar     = cbind(model = m, ans$covar)
    )
  })
  list(
    port_var  = rbindlist(lapply(res, `[[`, "port_var"),  use.names = TRUE, fill = TRUE),
    asset_var = rbindlist(lapply(res, `[[`, "asset_var"), use.names = TRUE, fill = TRUE),
    covar     = rbindlist(lapply(res, `[[`, "covar"),     use.names = TRUE, fill = TRUE)
  )
}

# ---- example usage ------------------------------------------------------------
# Load data once
data <- import_data(drop_first_col = FALSE, n_assets = 5)

# Load as many 'out' objects as you want and name them for comparison
out_list <- list(
  #smc126_5   = readRDS("empirical_results/standard_tip_w126_M2000_tipk5.rds"),
  #smc126_10   = readRDS("empirical_results/standard_tip_w126_M2000_tipk10.rds"),
  #smc126_15   = readRDS("empirical_results/standard_tip_w126_M2000_tipk15.rds"),
   #smc252_5   = readRDS("empirical_results/standard_tip_w252_M2000_tipk5.rds"),
   smc252_10   = readRDS("empirical_results/standard_tip_w252_M2000_tipk10.rds"),
   #smc252_15   = readRDS("empirical_results/standard_tip_w252_M2000_tipk15.rds"),
   #smc504_5   = readRDS("empirical_results/standard_tip_w504_M2000_tipk5.rds"),
   #smc504_10   = readRDS("empirical_results/standard_tip_w504_M2000_tipk10.rds"),
   #smc504_15   = readRDS("empirical_results/standard_tip_w504_M2000_tipk15.rds")
   
   smc252_10_flat   = readRDS("empirical_results/standard_5d_flat.rds"),
   alt_gaussian   = readRDS("empirical_results/naive_5d_gaussian.rds"),
   alt   = readRDS("empirical_results/out_naive_3_5d_extend_t.rds")
)

# Choose period & assets by name (Date handling is centralized and consistent)
date_from <- as.Date("2004-01-01")
date_to   <- as.Date("2025-01-01")
assets_sel <- NULL            # or c("AIG","AXP","BAC") to force a set/order

#out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/test.rds")
#ans <- evaluate_model(out, data, date_from, date_to)

# Run comparison
cmp <- compare_models(
  out_list_named = out_list,
  data           = data,
  date_from      = date_from,
  date_to        = date_to,
  assets         = assets_sel,
  alphas_var     = c(0.10, 0.05, 0.025, 0.01),
  alphas_covar   = c(0.10, 0.05)
)

cmp$covar[(alpha_j == 0.05)&(alpha_port  == 0.05)]

# Tidy results ready for tables/plots:
# cmp$port_var  # columns: model, alpha, n, hits, rate, kupiec, ind, cc
# cmp$asset_var # columns: model, asset, alpha, n, hits, rate, kupiec, ind, cc
df <- (cmp$covar[(alpha_j == 0.05)&(alpha_port  == 0.05)&(asset  == "AIG")][,c("model", "rate","kupiec","cc")])    # columns: model, asset, alpha_j, alpha_port, T_event, rate, kupiec, ind, cc
df$rate <- round(df$rate, 3)
df$kupiec <- round(df$kupiec, 3)
df$cc <- round(df$cc, 3)
df

df2 <- (cmp$covar[(alpha_j == 0.1)&(alpha_port  == 0.1)&(asset  == "AIG")][,c("model", "rate","kupiec","cc")])     # columns: model, asset, alpha_j, alpha_port, T_event, rate, kupiec, ind, cc
df2$rate <- round(df2$rate, 3)
df2$kupiec <- round(df2$kupiec, 3)
df2$cc <- round(df2$cc, 3)
df2

merged_df <- cbind(df[,2:4],df2[,2:4])
# Create the new data.frame
new_df <- data.frame(
  W = c(rep(126, 3), rep(252, 3), rep(504, 3)),
  k = c(5, 10, 15, 5, 10, 15, 5, 10, 15)
)
# Combine (new_df first)
merged_df <- cbind(new_df, merged_df)
aig <- merged_df























library(dplyr)
library(ggplot2)
library(tidyr)

aig <- aig[, -c(3, 4, 5)]
axp <- axp[, -c(3, 4, 5)]
bac <- bac[, -c(3, 4, 5)]
c <- c[, -c(3, 4, 5)]
cof <- cof[, -c(3, 4, 5)]


aig <- aig[, 1:5]
axp <- axp[, 1:5]
bac <- bac[, 1:5]
c <- c[, 1:5]
cof <- cof[, 1:5]


# --- Combine all data frames and label each asset
df_all <- rbind(
  cbind(asset = "AIG", aig[, c("W", "k", "rate")]),
  cbind(asset = "AXP", axp[, c("W", "k", "rate")]),
  cbind(asset = "BAC", bac[, c("W", "k", "rate")]),
  cbind(asset = "C",   c[,   c("W", "k", "rate")]),
  cbind(asset = "COF", cof[, c("W", "k", "rate")])
)

# Ensure correct types
df_all$asset <- factor(df_all$asset, levels = c("AIG", "AXP", "BAC", "C", "COF"))
df_all$W <- factor(df_all$W)  # treat W as categorical for plotting

library(ggplot2)

# Heatmap
p_heatmap <- ggplot(df_all, aes(x = factor(k), y = factor(W), fill = rate)) +
  geom_tile(color = "white") +
  facet_wrap(~ asset) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = expression(paste("Sensitivity of Rate to W and k (", alpha, " = ", beta, " = 0.05)")),
    x = "k",
    y = "Look-back window (W)",
    fill = "Violation rate"
  ) +
  theme_minimal(base_size = 13)

# Line plot
p_line <- ggplot(df_all, aes(x = k, y = rate, color = factor(W), group = W)) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~ asset) +
  labs(
    title = expression(paste("Effect of W and k on Rate (", alpha, " = ", beta, " = 0.05)")),
    x = "k",
    y = "Violation rate",
    color = "W"
  ) +
  theme_minimal(base_size = 13)

# Save the heatmap
ggsave(
  filename = "rate_sensitivity_heatmap_10.png",
  plot = p_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)

# Save the line plot
ggsave(
  filename = "rate_sensitivity_lineplot_10.png",
  plot = p_line,
  width = 8,
  height = 6,
  dpi = 300
)

# --- Combine all data frames and label each asset
df_all <- rbind(
  cbind(asset = "AIG", aig),
  cbind(asset = "AXP", axp),
  cbind(asset = "BAC", bac),
  cbind(asset = "C",   c),
  cbind(asset = "COF", cof)
)


alpha <- 0.1
df_all <- df_all %>%
  mutate(reject_kupiec = kupiec < 0.1)
df_all %>%
  filter(reject_kupiec)
df_all %>%
  group_by(asset) %>%
  summarise(rejections = sum(reject_kupiec),
            total = n(),
            rejection_rate = mean(reject_kupiec))
ggplot(df_all, aes(x = factor(k), y = factor(W), fill = reject_kupiec)) +
  geom_tile(color = "white") +
  facet_wrap(~ asset) +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "white")) +
  labs(
    title = "Kupiec Test Rejections (α = 0.1)",
    x = "k", y = "W (look-back window)",
    fill = "Reject?"
  ) +
  theme_minimal(base_size = 13)


p_heatmap <- ggplot(df_all, aes(x = factor(k), y = factor(W), fill = rate)) +
  geom_tile(color = "white") +
  geom_text(
    data = subset(df_all, reject_kupiec),
    aes(label = "✗"),        # symbol shown when rejected
    color = "red", size = 6
  ) +
  facet_wrap(~ asset) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = expression(paste("Sensitivity of Rate to W and k (", alpha, " = ", beta, " = 0.05)")),
    x = "k",
    y = "Look-back window (W)",
    fill = "Violation rate"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = "rate_sensitivity_heatmap_10.png",
  plot = p_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)



# ===== Helper: build conditional CoVaR hit time series (base-R friendly) =====
covar_hits_timeseries <- function(out, data, j, a = 0.10, b = 0.10,
                                  date_from = as.Date("2004-01-01"),
                                  date_to   = as.Date("2100-01-01"),
                                  assets    = NULL) {
  al <- align_and_filter(out, data, date_from, date_to, keep_cols = assets)
  Y        <- al$Y_keep
  keep_idx <- al$keep_idx
  dates    <- al$dates_oos[keep_idx]
  d        <- ncol(Y); tickers <- colnames(Y)
  if (j < 1L || j > d) stop("j out of bounds")
  rp <- rowMeans(Y)
  
  alphas_out <- if (!is.null(out$cfg$alphas)) as.numeric(out$cfg$alphas) else stop("cfg$alphas missing")
  k   <- .which_alpha(alphas_out, a)
  lab <- paste0("a", a, "b", b)
  
  # subset model arrays by SAME date indices
  riskVaR <- out$risk$VaR[keep_idx, , , drop = FALSE]
  CoVaR   <- out$CoVaR_tail[keep_idx, , , drop = FALSE]
  
  # optional: align asset order to data
  if (!is.null(dimnames(riskVaR)[[2]])) {
    pos <- match(tickers, dimnames(riskVaR)[[2]])
    if (anyNA(pos)) stop("Tickers not found in riskVaR dimnames.")
    riskVaR <- riskVaR[, pos, , drop = FALSE]
    if (!is.null(dimnames(CoVaR)[[2]])) CoVaR <- CoVaR[, pos, , drop = FALSE]
  }
  
  VaRj  <- matrix(riskVaR[, , k], ncol = d); colnames(VaRj) <- tickers
  labs_avail <- dimnames(CoVaR)[[3]]
  if (!is.null(labs_avail) && !(lab %in% labs_avail)) {
    lab_alt <- paste0("a", formatC(a, format = "f", digits = 2),
                      "b", formatC(b, format = "f", digits = 2))
    if (lab_alt %in% labs_avail) lab <- lab_alt
  }
  CoVab <- CoVaR[, , lab, drop = FALSE][, , 1]; colnames(CoVab) <- tickers
  
  # conditioning mask & hits (NA on non-event days)
  mask_j <- Y[, j] <= VaRj[, j]
  hits_j <- rep(NA_real_, length(dates))
  if (any(mask_j)) hits_j[mask_j] <- as.numeric(rp[mask_j] <= CoVab[mask_j, j])
  
  # cumulative stats along observed events only
  idx_obs    <- which(mask_j)
  cum_hits   <- rep(NA_real_, length(dates))
  cum_trials <- rep(NA_real_, length(dates))
  if (length(idx_obs)) {
    cum_hits[idx_obs]   <- cumsum(hits_j[idx_obs])
    cum_trials[idx_obs] <- seq_along(idx_obs)
  }
  cum_prop <- ifelse(is.na(cum_hits) | is.na(cum_trials), NA_real_, cum_hits / cum_trials)
  
  list(
    Date      = dates,
    asset_j   = tickers[j],
    event     = mask_j,      # TRUE on conditioning-event days
    hit       = hits_j,      # 1/0 on event days; NA else
    cum_hits  = cum_hits,    # cumulative hits on observed timeline
    cum_trials= cum_trials,  # cumulative events count
    cum_prop  = cum_prop     # cumulative proportion among observed
  )
}

# ===== Example: compare two named models in out_list (base R plots) ===========
# data <- import_data(drop_first_col = FALSE, n_assets = 3)
# out_list <- list(
#   smc = readRDS(".../test.rds"),
#   alt1     = readRDS(".../out_naive_2_extend_t.rds")
# )

j <- 1
a <- 0.1; b <- 0.1
date_from <- as.Date("2004-01-01")
date_to   <- as.Date("2025-01-01")

ts_smc <- covar_hits_timeseries(out_list$smc, data, j, a, b, date_from, date_to)
ts_alt1     <- covar_hits_timeseries(out_list$alt1,     data, j, a, b, date_from, date_to)

# ---------- Plot 1: cumulative number of hits (step) ----------
op <- par(no.readonly = TRUE); on.exit(par(op))
par(mar = c(4,4,3,8), xaxs = "i")  # extra right margin for legend

x1 <- ts_smc$Date
y1 <- ts_smc$cum_hits
x2 <- ts_alt1$Date
y2 <- ts_alt1$cum_hits

ylim <- range(y1, y2, na.rm = TRUE); if (!is.finite(ylim[1])) ylim <- c(0,1)

plot(x1, y1, type = "s", lwd = 2, xlab = "Date", ylab = "Cumulative hits",
     main = sprintf("CoVaR cumulative hits | j=%s, a=%.2f, b=%.2f", ts_smc$asset_j, a, b),
     ylim = ylim)
lines(x2, y2, type = "s", lwd = 2, lty = 2)
legend("right", inset = c(-0.02, 0), xpd = NA, bty = "n", lwd = 2, lty = c(1,2),
       legend = c("smc", "alt1"))

# ---------- Plot 2: cumulative conditional proportion (step) ----------
y1p <- ts_smc$cum_prop
y2p <- ts_alt1$cum_prop
ylim2 <- c(0, 0.2)

plot(x1, y1p, type = "s", lwd = 2, xlab = "Date",
     ylab = "Cumulative hit proportion (conditional)",
     main = "CoVaR cumulative conditional hit proportion",
     ylim = ylim2)
lines(x2, y2p, type = "s", lwd = 2, lty = 2, col="red")
abline(h = b, col = "gray60", lty = 3)  # benchmark proportion



idx_only_alt <- which(ts_alt1$hit == 1L & ts_smc$hit == 0L)
points(x1[idx_only_alt], y2p[idx_only_alt], pch = 24, bg = "red", col = "red")
# vertical gap between smc and alt1 at those times
segments(x1[idx_only_alt], y1p[idx_only_alt],
         x1[idx_only_alt], y2p[idx_only_alt],
         col = "red", lty = 3)

# subtle rug for quick scan
rug(x1[idx_only_alt], side = 3, col = "red")

legend("bottomright", inset = c(-0.02, 0), xpd = NA, bty = "n",
       lwd = 2, lty = c(1,2), col = c("black","red"),
       pch = c(NA, 24), pt.bg = c(NA, "red"),
       legend = c("smc", "alt1 (▲ = alt1-only hits)"))


# ---------- Plot 3 (optional): events vs hits summary (barplot) ----------
totals <- rbind(
  data.frame(model = "smc",
             events = max(ts_smc$cum_trials, na.rm = TRUE),
             hits   = max(ts_smc$cum_hits,   na.rm = TRUE)),
  data.frame(model = "alt1",
             events = max(ts_alt1$cum_trials, na.rm = TRUE),
             hits   = max(ts_alt1$cum_hits,   na.rm = TRUE))
)

# avoid -Inf from all-NA (no events case)
totals$events[!is.finite(totals$events)] <- 0
totals$hits[!is.finite(totals$hits)]     <- 0
misses <- pmax(totals$events - totals$hits, 0)

mat <- rbind(Hits = totals$hits, Misses = misses)
barplot(mat, beside = FALSE, col = c("gray20", "gray70"),
        main = "CoVaR: Total conditioning events vs hits",
        ylab = "Count", ylim = c(0, max(colSums(mat))*1.15))
text(x = 1:ncol(mat), y = colSums(mat)*1.05, labels = paste0("Events: ", colSums(mat)))
axis(1, at = 1:ncol(mat), labels = totals$model, tick = FALSE, line = -0.5)



#########################

# Gruped plots

# --- choose the three institutions you want to show ---------------------------
j_set <- c(1, 3, 5)   # <- pick the 3 j you want (indices of institutions)
a <- 0.1; b <- 0.1
date_from <- as.Date("2004-01-01")
date_to   <- as.Date("2025-01-01")
tickers <- setdiff((colnames(data$U)  ), "date")


png("covar_hits_j123.png", width = 1000, height = 800)

op <- par(
  mfrow = c(3, 1),
  mar   = c(3.2, 4.2, 2.2, 1.2),
  oma   = c(0, 0, 1.5, 0),
  cex.axis = 1.5   # <- enlarge x and y tick labels
)

for (i in seq_along(j_set)) {
  j <- j_set[i]
  
  ts_smc  <- covar_hits_timeseries(out_list$smc,  data, j, a, b, date_from, date_to)
  ts_alt1 <- covar_hits_timeseries(out_list$alt1, data, j, a, b, date_from, date_to)
  
  x1  <- ts_smc$Date
  x2  <- ts_alt1$Date
  y1p <- ts_smc$cum_prop
  y2p <- ts_alt1$cum_prop
  ylim2 <- c(0, 0.2)
  
  main_lbl <- if (!is.null(tickers)) {
    paste0("CoVaR cumulative conditional hit proportion — j = ", tickers[j])
  } else paste0("CoVaR cumulative conditional hit proportion — j = ", j)
  
  plot(x1, y1p, type = "s", lwd = 3,
       xlab = if (i == length(j_set)) "Date" else "",
       ylab = "",
       ylim = ylim2)
  
  lines(x2, y2p, type = "s", lwd = 3, lty = 2, col = "grey")
  abline(h = b, col = "gray60", lty = 3)
  
  idx_only_alt <- which(ts_alt1$hit == 1L & ts_smc$hit == 0L)
  if (length(idx_only_alt)) {
    rug(x1[idx_only_alt], side = 1, col = "red", lwd = 3, cex = 1.5)
  }
  
  if (i == 1) {
    legend("bottom", xpd = NA, bty = "n",
           lwd = 3, lty = c(1, 2), col = c("black", "grey"),
           legend = c("smc", "alt1"),
           cex = 1.8)  # make legend text larger
  }
}

if (!is.null(tickers)) {
  used_lbl <- paste(tickers[j_set], collapse = ", ")
} else {
  used_lbl <- paste("j =", paste(j_set, collapse = ", "))
}


mtext(paste("CoVaR cumulative conditional hit proportion - SMC vs ALT1. Tickers:", used_lbl), outer = TRUE, cex = 1.3)
par(op)

dev.off()




##########################

data$U$date
out = out_list$smc
# fam_hist: array [particles x time x edges]
fh <- out$fam_hist
dim(fh)
# [1] 1000 6289    3

n_edges <- dim(fh)[3]
n_time  <- dim(fh)[2]

# All labels present in the array
labels_all <- sort(unique(as.vector(fh)))

# ---- compute proportions -----------------------------------------------------
# result: list of length n_edges, each element is matrix [time x labels]
prop_list <- vector("list", n_edges)

for (e in 1:n_edges) {
  # slice array for edge e: matrix [particles x time]
  slice <- fh[, , e]
  
  # initialize matrix: n_time rows, one col per label
  mat <- matrix(NA_real_, nrow = n_time, ncol = length(labels_all))
  colnames(mat) <- labels_all
  
  # loop over time
  for (t in 1:n_time) {
    counts <- table(factor(slice[, t], levels = labels_all))
    mat[t, ] <- counts / sum(counts)
  }
  
  prop_list[[e]] <- mat
}

names(prop_list) <- paste0("edge", 1:n_edges)


op <- par(no.readonly = TRUE)
par(mfrow = c(3,1), mar = c(4,4,2,1))

time_index <- data$U$date
cols <- 1:length(labels_all)  # simple color set

for (e in 1:n_edges) {
  mat <- prop_list[[e]]
  plot(time_index, mat[,1], type = "l", ylim = c(0,1),
       xlab = "Time", ylab = "Proportion",
       main = paste("Edge", e),
       col = cols[1], lwd = 2)
  if (ncol(mat) > 1) {
    for (j in 2:ncol(mat)) {
      lines(time_index, mat[,j], col = cols[j], lwd = 2)
    }
  }
  legend("topright", legend = colnames(mat), col = cols, lwd = 2, bty = "n")
}

par(op)




#######################################


# --------- scoring primitives (vectorized) ----------
qs_quantile <- function(y, q, alpha) {
  # check loss / pinball loss
  (alpha - as.numeric(y <= q)) * (y - q)
}
qs_LM <- function(y, q) {
  # Lopez-Mancini
  ifelse(y < q, 1 + (y - q)^2, 0)
}
qs_LA <- function(y, q, alpha, ahat) {
  # Lopez-Aragon with inflation factor P based on empirical coverage ahat
  P <- if (is.finite(ahat) && ahat > alpha) exp((ahat - alpha) / alpha) else 1
  ifelse(y < q, P * (1 + abs(y - q)), abs(y - q))
}

# --------- build conditional CoVaR evaluation sample ----------
# Returns dates, mask of conditioning events, portfolio y, CoVaR threshold q
covar_eval_sample <- function(out, data, j, a = 0.10, b = 0.10,
                              date_from = as.Date("2004-01-01"),
                              date_to   = as.Date("2100-01-01"),
                              assets    = NULL) {
  al <- align_and_filter(out, data, date_from, date_to, keep_cols = assets)
  Y        <- al$Y_keep
  keep_idx <- al$keep_idx
  dates    <- al$dates_oos[keep_idx]
  d        <- ncol(Y)
  tickers  <- colnames(Y)
  if (j < 1L || j > d) stop("j out of bounds")
  
  rp <- rowMeans(Y)
  
  alphas_out <- if (!is.null(out$cfg$alphas)) as.numeric(out$cfg$alphas) else stop("cfg$alphas missing")
  k   <- .which_alpha(alphas_out, a)
  lab <- paste0("a", a, "b", b)
  
  riskVaR <- out$risk$VaR[keep_idx, , , drop = FALSE]
  CoVaR   <- out$CoVaR_tail[keep_idx, , , drop = FALSE]
  
  # align asset order if dimnames exist
  if (!is.null(dimnames(riskVaR)[[2]])) {
    pos <- match(tickers, dimnames(riskVaR)[[2]])
    if (anyNA(pos)) stop("Tickers not found in riskVaR dimnames.")
    riskVaR <- riskVaR[, pos, , drop = FALSE]
    if (!is.null(dimnames(CoVaR)[[2]])) CoVaR <- CoVaR[, pos, , drop = FALSE]
  }
  
  # fallback for labels like "a0.10b0.10"
  labs_avail <- dimnames(CoVaR)[[3]]
  if (!is.null(labs_avail) && !(lab %in% labs_avail)) {
    lab_alt <- paste0("a", formatC(a, format = "f", digits = 2),
                      "b", formatC(b, format = "f", digits = 2))
    if (lab_alt %in% labs_avail) lab <- lab_alt
  }
  
  VaRj  <- matrix(riskVaR[, , k], ncol = d); colnames(VaRj) <- tickers
  CoVab <- CoVaR[, , lab, drop = FALSE][, , 1]; colnames(CoVab) <- tickers
  
  mask <- Y[, j] <= VaRj[, j]       # conditioning events (asset j at level a)
  q    <- CoVab[, j]                # CoVaR threshold for portfolio
  y    <- rp                        # realized portfolio return
  
  list(dates = dates, mask = mask, y = y, q = q,
       asset_j = tickers[j], alpha_j = a, alpha_p = b)
}

# --------- compute CoVaR score time series + summary ----------
covar_scores <- function(out, data, j, a = 0.10, b = 0.10,
                         date_from = as.Date("2005-01-01"),
                         date_to   = as.Date("2100-01-01"),
                         assets    = NULL) {
  
  ev <- covar_eval_sample(out, data, j, a, b, date_from, date_to, assets)
  dates <- ev$dates; mask <- ev$mask; y <- ev$y; q <- ev$q
  
  # restrict to conditioning days to compute ahat (empirical conditional coverage)
  hits <- as.numeric(y[mask] <= q[mask])
  ahat <- if (any(mask)) mean(hits) else NA_real_
  
  # build time series (NA outside event days)
  qs  <- rep(NA_real_, length(dates))
  lm  <- rep(NA_real_, length(dates))
  la  <- rep(NA_real_, length(dates))
  if (any(mask)) {
    qs[mask] <- qs_quantile(y[mask], q[mask], b)   # quantile score at level b
    lm[mask] <- qs_LM(y[mask], q[mask])
    la[mask] <- qs_LA(y[mask], q[mask], b, ahat)
  }
  
  # aggregates over event days
  agg <- list(
    T_event = sum(mask),
    mean_qs = if (any(mask)) mean(qs[mask]) else NA_real_,
    mean_LM = if (any(mask)) mean(lm[mask]) else NA_real_,
    mean_LA = if (any(mask)) mean(la[mask]) else NA_real_,
    alpha_target = b,
    alpha_hat    = ahat
  )
  
  list(
    dates = dates, mask = mask, y = y, q = q,
    quantile_score = qs, LM = lm, LA = la,
    summary = agg,
    meta = list(asset_j = ev$asset_j, alpha_j = ev$alpha_j, alpha_p = ev$alpha_p)
  )
}


covar_scores_compare <- function(out_list_named, data, j, a = 0.10, b = 0.10,
                                 date_from = as.Date("2005-01-01"),
                                 date_to   = as.Date("2100-01-01"),
                                 assets    = NULL) {
  stopifnot(!is.null(names(out_list_named)), all(nzchar(names(out_list_named))))
  res <- lapply(names(out_list_named), function(m) {
    s <- covar_scores(out_list_named[[m]], data, j, a, b, date_from, date_to, assets)
    data.frame(
      model     = m,
      asset_j   = s$meta$asset_j,
      alpha_j   = s$meta$alpha_j,
      alpha_p   = s$meta$alpha_p,
      T_event   = s$summary$T_event,
      alpha_hat = s$summary$alpha_hat,
      mean_qs   = s$summary$mean_qs,
      mean_LM   = s$summary$mean_LM,
      mean_LA   = s$summary$mean_LA,
      row.names = NULL
    )
  })
  do.call(rbind, res)
}

# s1 <- covar_scores(out_list$smc, data, j=2, a=0.10, b=0.10,
#                    date_from=as.Date("2020-01-01"), date_to=as.Date("2022-01-01"))
# s2 <- covar_scores(out_list$alt1,     data, j=2, a=0.10, b=0.10,
#                    date_from=as.Date("2020-01-01"), date_to=as.Date("2022-01-01"))

plot_scores_series <- function(s, main_tag="model") {
  # step plots only on event days (NA elsewhere)
  op <- par(no.readonly=TRUE); on.exit(par(op))
  par(mfrow=c(3,1), mar=c(4,4,2,1))
  with(s, {
    plot(dates, quantile_score, type="s", lwd=2, xlab="Date",
         ylab="Quantile score", main=paste("Quantile score |", main_tag))
    abline(h=0, lty=3, col="gray60")
    plot(dates, LM, type="s", lwd=2, xlab="Date", ylab="LM",
         main=paste("LM |", main_tag))
    plot(dates, LA, type="s", lwd=2, xlab="Date", ylab="LA",
         main=paste("LA |", main_tag))
  })
}


# data <- import_data(drop_first_col = FALSE, n_assets = 3)
# out_list <- list(
#   smc = readRDS(".../test.rds"),
#   alt1     = readRDS(".../out_naive_2_extend_t.rds")
# )

# Single model, scores:
s1 <- covar_scores(out_list$smc, data, j=2, a=0.10, b=0.10,
                   date_from=as.Date("2003-01-01"), date_to=as.Date("2025-01-01"))
str(s1$summary)

# Compare models on mean scores:
tbl <- covar_scores_compare(out_list, data, j=2, a=0.10, b=0.10,
                            date_from=as.Date("2003-01-01"), date_to=as.Date("2025-01-01"))
print(tbl)

# Optional simple plots:
# plot_scores_series(s1, main_tag="smc")




##############################################
# ========= cumulative LM/LA plots for ONE asset j, comparing smc vs alt1 ======
# Settings ---------------------------------------------------------------------
# ========= cumulative LM & LA in ONE figure for ONE asset j ===================
# Settings ---------------------------------------------------------------------
j_sel     <- 2                  # <-- pick the conditioning asset index
a <- 0.10; b <- 0.10            # CoVaR^{b | j at level a}
date_from <- as.Date("2004-01-01")
date_to   <- as.Date("2025-01-01")
assets_sel <- NULL              # or c("AIG","AXP","BAC", ...)

# Use your existing score builder ------------------------------------------------
s_smc  <- covar_scores(out_list$smc,  data, j = j_sel, a = a, b = b,
                       date_from = date_from, date_to = date_to, assets = assets_sel)
s_alt1 <- covar_scores(out_list$alt1, data, j = j_sel, a = a, b = b,
                       date_from = date_from, date_to = date_to, assets = assets_sel)

# Align two date series and build cumulative sums (treat NA as 0 between events)
build_cum <- function(dates1, x1, dates2, x2) {
  DT <- data.table(Date = sort(unique(c(dates1, dates2))))
  setkey(DT, Date)
  DT[, smc  := NA_real_]
  DT[, alt1 := NA_real_]
  DT[.(dates1), smc  := x1]
  DT[.(dates2), alt1 := x2]
  DT[, cum_smc  := cumsum(replace(smc,  is.na(smc),  0))]
  DT[, cum_alt1 := cumsum(replace(alt1, is.na(alt1), 0))]
  DT[]
}

cum_LM <- build_cum(s_smc$dates, s_smc$LM, s_alt1$dates, s_alt1$LM)
cum_LA <- build_cum(s_smc$dates, s_smc$LA, s_alt1$dates, s_alt1$LA)

# Labels
asset_name <- s_smc$meta$asset_j
main_tag   <- sprintf("j = %s | a = %.2f, b = %.2f", asset_name, a, b)

# ---------- ONE FIGURE with two stacked panels ---------------------------------
png("covar_cum_LM_LA_smc_vs_alt1.png", width = 1200, height = 900)
op <- par(no.readonly = TRUE)
on.exit({par(op); dev.off()}, add = TRUE)

par(mfrow = c(2, 1),
    mar   = c(4.2, 5.2, 3.0, 10),  # extra right margin for legend on top panel
    oma   = c(0, 0, 2.0, 0),
    xaxs  = "i",
    cex.axis = 1.5, cex.lab = 1.6)

## Panel 1: cumulative LM
with(cum_LM, {
  ylim <- range(c(cum_smc, cum_alt1), na.rm = TRUE)
  if (!is.finite(ylim[1])) ylim <- c(0, 1)
  plot(Date, cum_smc, type = "s", lwd = 3,
       xlab = "", ylab = "Cumulative LM",
       main = paste("Cumulative LM —", main_tag),
       ylim = ylim)
  lines(Date, cum_alt1, type = "s", lwd = 3, lty = 2, col = "grey40")
})
legend("bottom", inset = c(-0.03, 0), xpd = NA, bty = "n",
       lwd = 3, lty = c(1, 2), col = c("black", "grey40"),
       legend = c("smc", "alt1"), cex = 1.5)

## Panel 2: cumulative LA
with(cum_LA, {
  ylim <- range(c(cum_smc, cum_alt1), na.rm = TRUE)
  if (!is.finite(ylim[1])) ylim <- c(0, 1)
  plot(Date, cum_smc, type = "s", lwd = 3,
       xlab = "", ylab = "Cumulative LA",
       main = paste("Cumulative LA —", main_tag),
       ylim = ylim)
  lines(Date, cum_alt1, type = "s", lwd = 3, lty = 2, col = "grey40")
})

mtext("CoVaR conditional scoring - AXP",
      outer = TRUE, cex = 1.8)
# ------------------------------------------------------------------------------


dev.off()







############################################################
# test rosenblatt transformations

library(dplyr)
library(tidyr)

# --- helper: run uniformity tests on one numeric vector
have_goftest <- requireNamespace("goftest", quietly = TRUE)

test_uniform <- function(x) {
  x <- x[is.finite(x)]
  p_ks <- tryCatch(stats::ks.test(x, "punif")$p.value, error = function(e) NA_real_)
  p_ad <- if (have_goftest) {
    tryCatch(goftest::ad.test(x, null = "punif")$p.value, error = function(e) NA_real_)
  } else NA_real_
  c(p_ks = p_ks, p_ad = p_ad)
}

# --- loop over (W, k) combos and build one tidy table
results <- list()
for (W in c(126, 252, 504)) {
  for (k in c(5, 10, 15)) {
    key <- paste0("smc", W, "_", k)
    mat <- out_list[[key]]$pit_rose               # matrix with 5 columns (assets)
    
    # set fallback asset names if none
    if (is.null(colnames(mat))) {
      colnames(mat) <- paste0("Asset", seq_len(ncol(mat)))
    }
    
    tmp <- apply(mat, 2, test_uniform)            # 2 = over columns
    df  <- as.data.frame(t(tmp))
    df$asset <- rownames(df)
    rownames(df) <- NULL
    df$W <- W; df$k <- k
    results[[key]] <- df
  }
}

pvals_long <- bind_rows(results) %>%
  relocate(W, k, asset, p_ks, p_ad) %>%
  mutate(reject_ks = p_ks < 0.05,
         reject_ad = ifelse(is.na(p_ad), NA, p_ad < 0.05))

# --- Wide tables (one row per (W,k), columns = assets)
pvals_ks <- pvals_long %>%
  select(W, k, asset, p_ks) %>%
  mutate(p_ks = round(p_ks, 3)) %>%
  tidyr::pivot_wider(names_from = asset, values_from = p_ks) %>%
  arrange(W, k)

pvals_ad <- pvals_long %>%
  select(W, k, asset, p_ad) %>%
  mutate(p_ad = round(p_ad, 3)) %>%
  tidyr::pivot_wider(names_from = asset, values_from = p_ad) %>%
  arrange(W, k)

# --- Print the main KS table (common choice for PIT uniformity)
print(pvals_ks)

# --- If you want LaTeX tables (optional):
# knitr::kable(pvals_ks, format = "latex", booktabs = TRUE,
#              caption = "PIT Uniformity p-values (KS) by W and k")
# knitr::kable(pvals_ad, format = "latex", booktabs = TRUE,
#              caption = "PIT Uniformity p-values (AD) by W and k")

# --- Quick rejection heatmap (optional)
# ggplot(pvals_long, aes(x = factor(k), y = factor(W), fill = reject_ks)) +
#   geom_tile(color = "white") +
#   facet_wrap(~ asset) +
#   scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey90", "NA"="grey70")) +
#   labs(title = "KS Rejections of PIT Uniformity (alpha = 0.05)",
#        x = "k", y = "W", fill = "Reject?") +
#   theme_minimal(base_size = 13)
