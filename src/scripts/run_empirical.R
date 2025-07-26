library(here)
source(here("src/R", "config.R"))         



run_empirical <- function() {
  
  dat <- list(
    U      = readRDS("data/PIT.rds"),
    mu_fc  = readRDS("data/returns_mean_forecast.rds")[,-1],  # drop date col
    sig_fc = readRDS("data/returns_vol_forecast.rds")[,-1]
  )
  
  cfg_variants <- list(list(label = "test"))          
  for (v in cfg_variants) {
    cfg <- modifyList(build_cfg(ncol(dat$U)), v[ setdiff(names(v),"label") ])
    cfg$label <- v$label %||% "cfg"
    res <- smc_full(dat, cfg)
    res$alg <- "standard"
    save_result(res, here("empirical_results"))
  }
}

run_empirical()











