library(here)
source(here("src/R", "config.R"))         



run_empirical <- function() {
  
  dat <- list(
    U      = readRDS("data/PIT.rds")[,1:3],
    mu_fc  = readRDS("data/returns_mean_forecast.rds")[,2:4],# [,-1],  # drop date col
    sig_fc = readRDS("data/returns_vol_forecast.rds")[,2:4],  #[,-1],
    df_fc = readRDS("data/df_fc.rds")[,2:4],#[,-1]
    shape_fc = readRDS("data/shape_fc.rds")[,2:4]
  )
  
  cfg_variants <- list(list(label = "test"))   
  v=cfg_variants[1]
  for (v in cfg_variants) {
    cfg <- modifyList(build_cfg(ncol(dat$U)), v[ setdiff(names(v),"label") ])
    cfg$label <- v$label %||% "cfg"
    set.seed(cfg$seed)
    res <- smc_full(dat, cfg)
    res$alg <- "standard"
    save_result(res, here("empirical_results"))
  }
}

run_empirical()











