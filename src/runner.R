# ──────────────────────────────────────────────────────────────────────────────
#  experiment_runner.R  ── batch-run C-vine SMC (standard & block)  ── 2025-06-18
# ──────────────────────────────────────────────────────────────────────────────

library(here)
source(here::here("src","config.R"))


dir.create(here::here("simul_results/static_dgp"), showWarnings = FALSE)

run_and_save <- function(U, cfg, alg = c("standard", "block"), tag = NULL) {
  alg <- match.arg(alg)
  
  res <- switch(
    alg,
    standard = run_standard_smc(U, cfg, type = "standard"),
    block    = run_block_smc(U, cfg, type = "block")
  )
  
  res$cfg <- cfg                                   
  
  tag   <- if (is.null(tag)) "" else paste0("_", tag)
  fname <- here::here("simul_results/static_dgp", paste0(alg, "_", ncol(U) , tag))

  saveRDS(res, fname)
  message(sprintf("✓ Saved %s", fname))
  
  cat("\n\n===== FINAL MODEL EVALUATION =====\n")
  cat(sprintf("Log Model Evidence: %.4f\n", res$log_model_evidence))
  
  rm(res); invisible(gc())                         
}


set.seed(126)
U  <- sim_static_cop_3(N = 200)              
d  <- ncol(U)

cfg_variants <- list(
  # list(
  #   label      = "nmh1_N200_fixed",
  #   n_mh       = 1,
  #   tau_prior  = "fixed"
  # ) ,
  list(
    label      = "nmh1_N200_ivgamma",
    n_mh       = 1,
    tau_prior  = "inv_gamma"
  )
)


system.time(
for (i in seq_along(cfg_variants)) {
  v    <- cfg_variants[[i]]          # pull the i-th inner list
  tag  <- v[["label"]]               # safe even if names = NULL
  
  tweaks <- v[ setdiff(names(v), "label") ]   # drop label for merging
  cfg    <- modifyList(build_cfg(d), tweaks)
  cfg$label <- tag
  
  for (alg in c("standard")) {
    set.seed(cfg$seed)
    run_and_save(U, cfg, alg, tag)
  }
})






