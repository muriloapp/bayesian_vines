# ──────────────────────────────────────────────────────────────────────────────
#  experiment_runner.R  ── batch-run C-vine SMC (standard & block)  ── 2025-06-18
# ──────────────────────────────────────────────────────────────────────────────

library(here)
source(here::here("src","config.R"))

library(profvis)



dir.create(here::here("simul_results/static_dgp"), showWarnings = FALSE)

run_and_save <- function(data, cfg, alg = c("standard", "block"), tag = NULL) {
  alg <- match.arg(alg)
  
  res <- switch(
    alg,
    standard = run_standard_smc(data, cfg, type = "standard"),
    block    = run_block_smc(data, cfg, type = "block")
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


set.seed(42)
data  <- sim_static_cop_3(N = 20)    
U <- data$U
d  <- ncol(U)

cfg_variants <- list(
  # list(
  #   label      = "nmh1_N200_fixed",
  #   n_mh       = 1,
  #   tau_prior  = "fixed"
  # ) ,
  list(
    label      = "N250_3F",
    tau_prior  = "inv_gamma",
    pi_prior   = "beta"
    #M           = 2000
  )
)


system.time(
for (i in seq_along(cfg_variants)) {
  v    <- cfg_variants[[i]]         
  tag  <- v[["label"]]               
  
  tweaks <- v[ setdiff(names(v), "label") ]   
  cfg    <- modifyList(build_cfg(d), tweaks)
  cfg$label <- tag
  
  for (alg in c("standard")) {
    set.seed(cfg$seed)
    #p <- profvis::profvis({
    run_and_save(data, cfg, alg, tag)
    #})
  }
})

# htmlwidgets::saveWidget(p, "profile.html")
# browseURL("profile.html")
