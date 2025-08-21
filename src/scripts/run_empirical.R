library(here)
source(here("src/R", "config.R"))         


n_assets <- 1:7
n_days <- 1:3000

run_empirical <- function() {
  
  dat <- list(
    U      = readRDS("data/PIT.rds")[1:length(n_days),n_assets],
    mu_fc  = readRDS("data/returns_mean_forecast.rds")[n_days,(n_assets+1), with = FALSE],# [,-1],  # drop date col
    sig_fc = readRDS("data/returns_vol_forecast.rds")[n_days,(n_assets+1), with = FALSE],  #[,-1],
    df_fc = readRDS("data/df_fc.rds")[n_days,(n_assets+1), with = FALSE],#[,-1]
    shape_fc = readRDS("data/shape_fc.rds")[n_days,(n_assets+1), with = FALSE],
    y_real = readRDS("data/returns_actual.rds")[n_days,n_assets+1, with = FALSE]
    
  )
  
  cfg_variants <- list(list(label = "test"))   
  v=cfg_variants[1]
  for (v in cfg_variants) {
    cfg <- modifyList(build_cfg(ncol(dat$U)), v[ setdiff(names(v),"label") ])
    cfg$label <- v$label %||% "cfg"
    set.seed(cfg$seed)
    
    data = dat
    res <- smc_full(dat, cfg)
    res$alg <- "standard"
    save_result(res, here("empirical_results"))
  }
}

run_empirical()











