# ──────────────────────────────────────────────────────────────────────────────
#  mc_experiment.R  ── full Monte-Carlo for C-vine SMC  ── 2025-07-02
# ──────────────────────────────────────────────────────────────────────────────

library(here)
source(here("src", "config.R"))          # build_cfg(), sim_static_cop_6(), …
# source(here("src", "smc_kernels.R"))   # run_standard_smc(), run_block_smc()

## ---------- experiment specs -------------------------------------------------
n_sim     <- 10
N_obs     <- 1500
seed_base <- 42
dims      <- c(3, 6)                     # <-- asked for both

cfg_variants <- list(
  list(label = "M2000", M = 2000),
  list(label = "M1000", M = 1000)
)

## ---------- helpers ----------------------------------------------------------
out_root <- here("simul_results", "mc_static_dgp")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

save_fit <- function(res, label, alg, sim_tag, d) {
  path <- file.path(out_root, sprintf("d%d", d), label)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file.path(path, sprintf("%s_%s.rds", alg, sim_tag)))
  message(sprintf("saved   d=%d | %s | %s | %s",
                  d, label, alg, sim_tag))
}

## ---------- Monte-Carlo ------------------------------------------------------
for (sim in seq_len(n_sim)) {
  for (d in dims) {
    
    set.seed(seed_base + d + sim)          # unique per (d, sim)
    sim_tag <- sprintf("sim%03d", sim)
    
    ## 1️⃣  simulate data once per replication
    data <- sim_static_cop(d = d, N = N_obs)
    data_dir <- file.path(out_root, sprintf("d%d", d), "_data")
    dir.create(data_dir,  recursive = TRUE, showWarnings = FALSE)
    saveRDS(data, file.path(data_dir,
                            sprintf("data_%s.rds", sim_tag)))
    
    ## 2️⃣  run all cfg variants
    for (v in cfg_variants) {
      label  <- v$label
      cfg    <- modifyList(build_cfg(d), v)      # M overwritten
      cfg$label <- label
      
      ## === choose algorithms here =====================================
      for (alg in "standard") {
        t_sec <- system.time(
          res <- run_standard_smc(data, cfg, type = "standard")
        )[["elapsed"]]
        
        res$cfg      <- cfg
        res$sim_tag  <- sim_tag
        res$alg      <- alg
        res$elapsed  <- t_sec
        res$d        <- d
        
        save_fit(res, label, alg, sim_tag, d)
      }
    }
    
    rm(data); gc(verbose = FALSE)
  }
}



# 
# pkgs <- unique(renv::dependencies()$Package)
# 
# # create a vector of "package==version"
# lines <- sapply(pkgs, function(pkg) {
#   ver <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "not_installed")
#   paste0(pkg, "==", ver)
# })
# 
# writeLines(lines, "package_versions.txt")
# 


library(here)

# all block-SMC fits for the “M200_ivgamma_beta” variant
block_files <- list.files(
  here("simul_results/mc_static_dgp/d3/M2000/"),
  pattern = "^standard_.*\\.rds$", full.names = TRUE
)
res <- lapply(block_files, readRDS)

# corresponding data sets
data_files <- list.files(
  here("simul_results/mc_static_dgp/d3/_data"),
  pattern = "^data_.*\\.rds$", full.names = TRUE
)
sim_data <- lapply(data_files, readRDS)




# theta_mean <- apply(res[[1]]$theta_hist, MARGIN = c(2, 3), FUN = mean)
# theta_sd <- apply(res[[1]]$theta_hist, MARGIN = c(2, 3), FUN = sd)
# plot_theta_paths(tanh(theta_mean), theta_sd, k=12, theta_true = -0.37)




## ---- 1) make sure res and sim_data align by sim tag -------------------------
# If your list.files() already returns in "sim001..sim010" order, you can skip this.
tag_from_path <- function(p) {
  m <- regmatches(p, regexpr("sim\\d{3}", p))
  ifelse(length(m), m, NA_character_)
}

## ---- 2) mapper from sim codes -> fam_hist codes -----------------------------
map_sim_to_smc <- function(code) {
  ifelse(code == 0, 0L,
         ifelse(code == 1, 1L,
                ifelse(code == 7, 2L, NA_integer_)))
}

# ## ---- 3) core function: recovery for edge 1 ----------------------------------
# recover_true_family_edge1 <- function(res_list, sim_list, edge_idx = 1L, burn_in = 0L) {
#   out <- rep(NA_real_, 10)
#   
#   for (i in 1:10) {
#     fam_mat <- sim_list[[i]]$family
#     true_sim_code <- fam_mat[2, 1]                    # C-vine: edge (1,2)
#     true_code     <- map_sim_to_smc(true_sim_code)
#     
#     fh <- res_list[[i]]$fam_hist[, , edge_idx, drop = FALSE][, , 1]  # [M x T]
#     if (is.na(true_code)) { out[i] <- NA_real_; next }
#     
#     Tt <- ncol(fh)
#     t_idx <- seq.int(max(1L, burn_in + 1L), Tt)
#     # posterior mass on the true family at each t
#     mass_t <- apply(fh[, t_idx, drop = FALSE], 2, function(col) mean(col == true_code, na.rm = TRUE))
#     # average over time
#     out[i] <- mean(mass_t, na.rm = TRUE)
#   }
#   out
# }
# 
# ## ---- 4) run it and summarize ------------------------------------------------
# recov <- recover_true_family_edge1(res, sim_data, edge_idx = 1L, burn_in = 0)
# print(recov)                    # one value per simulation (probability in [0,1])
# cat(sprintf("mean=%.3f, median=%.3f, sd=%.3f\n",
#             mean(recov, na.rm = TRUE), median(recov, na.rm = TRUE), sd(recov, na.rm = TRUE)))


burn_in = 1260L
edge_idx = 3L
out <- rep(NA_real_, 10)

for (i in 1:10) {
  fam_mat <- sim_data[[i]]$family
  #true_sim_code <- fam_mat[2, 1]                    # C-vine: edge (1,2)
  #true_sim_code <- fam_mat[3, 1]                    # C-vine: edge (1,2)
  true_sim_code <- fam_mat[3, 2] 
  true_code     <- map_sim_to_smc(true_sim_code)
  
  fh <- res[[i]]$fam_hist[, , edge_idx, drop = FALSE][, , 1]  # [M x T]
  if (is.na(true_code)) { out[i] <- NA_real_; next }
  
  Tt <- ncol(fh)
  t_idx <- seq.int(max(1L, burn_in + 1L), Tt)
  # posterior mass on the true family at each t
  mass_t <- apply(fh[, t_idx, drop = FALSE], 2, function(col) mean(col == true_code, na.rm = TRUE))
  # average over time
  out[i] <- mean(mass_t, na.rm = TRUE)
}


mean(out)










burn_in  <- 1260L
edge_idx <- 2L
row_idx  <- 3L  # matches your fam_mat[3, 2]
col_idx  <- 1L
min_count <- 10L

# ---- output table ----
out_tab <- data.frame(sim = 1:10,
                      family_true = NA_integer_,
                      rms = NA_real_,
                      rms_theta = NA_real_,  # BB1 only
                      rms_delta = NA_real_,  # BB1 only
                      T_used = NA_integer_,
                      avg_npart = NA_real_)

for (i in 1:10) {
  # truth from simulator
  fam_true   <- sim_data[[i]]$family[row_idx, col_idx]
  theta_true <- sim_data[[i]]$theta [row_idx, col_idx]
  delta_true <- sim_data[[i]]$theta2[row_idx, col_idx]
  code_true  <- map_sim_to_smc(fam_true)
  out_tab$family_true[i] <- code_true
  if (is.na(code_true) || code_true == 0L) next  # skip indep/unknown
  
  # SMC outputs (particles x time) for the chosen edge
  fh <- res[[i]]$fam_hist[, , edge_idx, drop = FALSE][, , 1]
  p1 <- res[[i]]$par1_hist[, , edge_idx, drop = FALSE][, , 1]
  p2 <- if (!is.null(res[[i]]$par2_hist)) {
    res[[i]]$par2_hist[, , edge_idx, drop = FALSE][, , 1]
  } else if (!is.null(res[[i]]$par1_hist_2)) {
    res[[i]]$par1_hist_2[, , edge_idx, drop = FALSE][, , 1]
  } else NULL
  
  Tt <- ncol(fh)
  t_idx <- seq.int(max(1L, burn_in + 1L), Tt)
  
  mse_all <- numeric(0); mse_theta <- numeric(0); mse_delta <- numeric(0); np_t <- integer(0)
  
  for (t in t_idx) {
    sel <- which(fh[, t] == code_true)
    if (length(sel) < min_count) next
    
    if (code_true == 1L) {
      # Gaussian: par1 vs true correlation (stored in theta)
      dif1 <- p1[sel, t] - theta_true
      mse_all <- c(mse_all, mean(dif1^2, na.rm = TRUE))
      np_t    <- c(np_t, length(sel))
    } else if (code_true == 2L) {
      # BB1: par1=theta, par2=delta
      if (is.null(p2) || !is.finite(delta_true)) next
      dth <- p1[sel, t] - theta_true
      ddl <- p2[sel, t] - delta_true
      mse_theta <- c(mse_theta, mean(dth^2, na.rm = TRUE))
      mse_delta <- c(mse_delta, mean(ddl^2, na.rm = TRUE))
      mse_all   <- c(mse_all, 0.5 * (mean(dth^2, na.rm = TRUE) + mean(ddl^2, na.rm = TRUE)))
      np_t      <- c(np_t, length(sel))
    }
  }
  
  if (length(mse_all)) {
    out_tab$rms[i]       <- sqrt(mean(mse_all,   na.rm = TRUE))
    out_tab$T_used[i]    <- length(mse_all)
    out_tab$avg_npart[i] <- mean(np_t)
  }
  if (length(mse_theta)) out_tab$rms_theta[i] <- sqrt(mean(mse_theta, na.rm = TRUE))
  if (length(mse_delta)) out_tab$rms_delta[i] <- sqrt(mean(mse_delta, na.rm = TRUE))
}

out_tab
mean(out_tab$rms[4:10], na.rm = TRUE)


