library(here)
source(here("src/R", "config.R"))         


n_assets <- 1:7
n_days <- 1:3000

run_empirical <- function() {
  
  # dat <- list(
  #   U      = readRDS("data/PIT.rds")[1:length(n_days),n_assets],
  #   mu_fc  = readRDS("data/returns_mean_forecast.rds")[n_days,(n_assets+1), with = FALSE],# [,-1],  # drop date col
  #   sig_fc = readRDS("data/returns_vol_forecast.rds")[n_days,(n_assets+1), with = FALSE],  #[,-1],
  #   df_fc = readRDS("data/df_fc.rds")[n_days,(n_assets+1), with = FALSE],#[,-1]
  #   shape_fc = readRDS("data/shape_fc.rds")[n_days,(n_assets+1), with = FALSE],
  #   y_real = readRDS("data/returns_actual.rds")[n_days,n_assets+1, with = FALSE]
  #   
  # )
  
  dat <- list(
    U      = readRDS("data/PIT.rds")[,n_assets],
    mu_fc  = readRDS("data/returns_mean_forecast.rds")[,(n_assets+1), with = FALSE],# [,-1],  # drop date col
    sig_fc = readRDS("data/returns_vol_forecast.rds")[,(n_assets+1), with = FALSE],  #[,-1],
    df_fc = readRDS("data/df_fc.rds")[,(n_assets+1), with = FALSE],#[,-1]
    shape_fc = readRDS("data/shape_fc.rds")[,(n_assets+1), with = FALSE],
    y_real = readRDS("data/returns_actual.rds")[,n_assets+1, with = FALSE]
    
  )
  
  #cfg_variants <- list(list(label = "test"))   
  cfg_variants <- list(
    list(label = "test"),
    list(label = "std",    use_weighted_ll = FALSE),
    list(label = "tailW_tauL0.2_taileps0.3",  use_weighted_ll = TRUE,
         tauL = 0.2, joint_k = 2L, tail_eps = 0.30),
    list(label = "tailW",  use_weighted_ll = TRUE,
         tauL = 0.05, joint_k = 2L, tail_eps = 0.10)
  )
  
  v=cfg_variants[1]
  for (v in cfg_variants) {
    cfg <- modifyList(build_cfg(ncol(dat$U)), v[ setdiff(names(v),"label") ])
    cfg$label <- v$label %||% "cfg"
    set.seed(cfg$seed)
    
    data = dat
    res <- smc_full(dat, cfg)
    res$alg <- if (isTRUE(cfg$use_weighted_ll)) "tail_weighted" else "standard"
    save_result(res, here("empirical_results"))
  }
}

run_empirical()











