
library(here)
library(future)
library(future.apply)

set.seed(11111)

plan(multisession, workers = max(1, parallel::detectCores() - 1))

source(here("src/R", "config.R"))
source(here("src/simulation/naive_simulation.R"))
source(here("src/simulation/main_simulation_nonparallel.R"))
source(here("src/simulation/fun_simulation.R"))




# ---- list ONLY the scenarios you want -----------------------------------------
scenarios <- list(
  #list(ml = 1e10, w = 252, wp = 1000),
  #list(ml = 500,  w = 252, wp = 1000),
  #list(ml = 500,  w = 126, wp = 1000),
  #list(ml = 1e10, w = 252, wp = 1000, re=1),
  #list(ml = 1e10, w = 252, wp = 1000, re=252, sim=100000),
  #list(ml = 1e10,  w = 252, wp = 1000,  re=252, sim=10000),
  list(ml = 1e10,  w = 252, wp = 1000,  re=252, sim=20000)
  #list(ml = 500,  w = 126, wp = 1000)
  #list(ml = 500,  w = 504, wp = 1000),
  # add more here, one line per combo
)

base_dir <- "simul_results/NAIVE"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)


make_tag <- function(s) {
  parts <- c()
  
  if (!is.null(s$ml)) parts <- c(parts, paste0("ml", s$ml))
  if (!is.null(s$w))  parts <- c(parts, sprintf("w%03d", s$w))
  if (!is.null(s$wp)) parts <- c(parts, sprintf("wp%04d", s$wp))
  if (!is.null(s$re)) parts <- c(parts, sprintf("re%03d", s$re))
  if (!is.null(s$sim)) parts <- c(parts, sprintf("sim%03d", s$sim))
  
  paste(parts, collapse = "_")
}

s <- scenarios[[1]]
for (s in scenarios) {
  
  # build cfg, using defaults when missing
  cfg_list <- list(M = 1000, label = "M1000", use_tail_informed_prior = TRUE, sim=s$sim)
  
  if (!is.null(s$wp)) cfg_list$W_predict <- s$wp
  if (!is.null(s$w))  cfg_list$W        <- s$w
  if (!is.null(s$re)) cfg_list$aic_refit_every <- s$re  # only if exists
  
  cfg <- modifyList(build_cfg(d = 2), cfg_list)
  
  n_train <- cfg$W_predict
  d       <- 2
  n_test  <- 20
  
  ml <- if (!is.null(s$ml)) s$ml else stop("scenario missing ml")  # or set a default

  tag <- make_tag(s)
  
  out_dir <- file.path(base_dir, tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n============================\n")
  cat("Running:", tag, "\n")
  cat("============================\n")
  
  sim_files <- future_lapply(
    X = 1:300,
    FUN = run_one_sim,
    n_train   = n_train,
    n_test    = n_test,
    d         = d,
    cfg       = cfg,
    mean_len  = ml,
    out_dir   = out_dir,
    future.seed = 11111
  )
  
  saveRDS(sim_files, file.path(out_dir, sprintf("sim_files_%s.rds", tag)))
}















# 
# 
# folder <- "simul_results/SMC/mlNA_w252_wp1000_re001"   # <- change this
# files  <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
# 
# # read all files (each file is assumed to be a list)
# obj_list <- lapply(files, readRDS)
# 
# # extract the element you want (file[[1]]$covar, file[[2]]$covar, ...)
# rmse_list <- lapply(obj_list, `[[`, "rmse_mae_from_covar")
# rmse_all <- do.call(rbind, rmse_list)
# 
# 
# with(rmse_all, mean(RMSE[cond_asset == 2 & scenario == "a0.1b0.1"], na.rm = TRUE))
# with(rmse_all, mean(MAE[cond_asset == 2 & scenario == "a0.1b0.1"], na.rm = TRUE))
# with(rmse_all, mean(RMSE[cond_asset == 2 & scenario == "a0.1b0.05"], na.rm = TRUE))
# with(rmse_all, mean(MAE[cond_asset == 2 & scenario == "a0.1b0.05"], na.rm = TRUE))
# with(rmse_all, mean(RMSE[cond_asset == 2 & scenario == "a0.05b0.05"], na.rm = TRUE))
# with(rmse_all, mean(MAE[cond_asset == 2 & scenario == "a0.05b0.05"], na.rm = TRUE))
# with(rmse_all, mean(RMSE[cond_asset == 2 & scenario == "a0.05b0.025"], na.rm = TRUE))
# with(rmse_all, mean(MAE[cond_asset == 2 & scenario == "a0.05b0.025"], na.rm = TRUE))
# 
# 
# 
# covar_list <- lapply(obj_list, `[[`, "eval_covar_asset")
# 
# covar_all <- do.call(rbind, covar_list)
# 
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.1 & alpha_port == 0.1], na.rm = TRUE))
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.1 & alpha_port == 0.05], na.rm = TRUE))
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.05 & alpha_port == 0.05], na.rm = TRUE))
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.05 & alpha_port == 0.025], na.rm = TRUE))
# 
# 
# 
# (covar_all[covar_all$asset==1 & covar_all$alpha_j==0.05 & covar_all$alpha_port==0.025,,drop=FALSE])
# 
# covar_all[covar_all$scenario =="a0.05b0.05",,drop=FALSE]
# 
# 
# 
# 
# 
# 
# 
# folder <- "simul_results/NAIVE/mlNA_w252_wp0252_re252"   # <- change this
# files  <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
# obj_list <- lapply(files, readRDS)
# covar_list <- lapply(obj_list, `[[`, "rmse_mae_from_covar")
# covar_list <- lapply(obj_list, `[[`, "eval_covar_asset")
# 
# covar_list <- Filter(Negate(is.null), covar_list)
# covar_all <- do.call(rbind, covar_list)
# 
# 
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.1 & alpha_port == 0.05], na.rm = TRUE))
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.05 & alpha_port == 0.05], na.rm = TRUE))
# with(covar_all, mean(rate[asset == 2 & alpha_j == 0.05 & alpha_port == 0.025], na.rm = TRUE))
# 
# 
# with(covar_all, mean(RMSE[cond_asset == 2 & scenario == "a0.1b0.05"], na.rm = TRUE))
# with(covar_all, mean(MAE[cond_asset == 2 & scenario == "a0.1b0.05"], na.rm = TRUE))
# with(covar_all, mean(RMSE[cond_asset == 2 & scenario == "a0.05b0.05"], na.rm = TRUE))
# with(covar_all, mean(MAE[cond_asset == 2 & scenario == "a0.05b0.05"], na.rm = TRUE))
# with(covar_all, mean(RMSE[cond_asset == 2 & scenario == "a0.05b0.025"], na.rm = TRUE))
# with(covar_all, mean(MAE[cond_asset == 2 & scenario == "a0.05b0.025"], na.rm = TRUE))
# 
# 
# 
# 
# 
# out2 <- covar_all[covar_all$cond_asset == 2 & covar_all$scenario =="a0.05b0.05",,drop=FALSE]$RMSE
# 
# 
# 



library(dplyr)
library(tidyr)
library(openxlsx)

folder_rel <- "NAIVE/ml1e+10_w252_wp1000_re252_sim100000"
#folder_rel <- "NAIVE_300/ml500_w252_wp0252_re063"

folder     <- file.path("simul_results", folder_rel)

scenarios_tbl <- tibble(
  alpha_j    = c(0.10, 0.10, 0.05, 0.05),
  alpha_port = c(0.10, 0.05, 0.05, 0.025)
) %>%
  mutate(
    scenario = sprintf("a%sb%s", alpha_j, alpha_port),
    Col      = sprintf("(%g,%g)", alpha_j, alpha_port)
  )

summarize_folder <- function(folder, cond_asset = 2, asset = 2, scenarios_tbl) {
  files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No .rds files found in: ", folder)
  
  obj <- lapply(files, readRDS)
  
  rmse_all  <- bind_rows(lapply(obj, `[[`, "rmse_mae_from_covar"))
  covar_all <- bind_rows(lapply(obj, `[[`, "eval_covar_asset"))
  
  rmse_sum <- rmse_all %>%
    filter(cond_asset == !!cond_asset) %>%
    group_by(scenario) %>%
    summarise(RMSE = mean(RMSE, na.rm = TRUE),
              MAE  = mean(MAE,  na.rm = TRUE), .groups = "drop")
  
  vr_sum <- covar_all %>%
    filter(asset == !!asset) %>%
    group_by(alpha_j, alpha_port) %>%
    summarise(ViolationRate = mean(rate, na.rm = TRUE), .groups = "drop") %>%
    mutate(scenario = sprintf("a%sb%s", alpha_j, alpha_port))
  
  scenarios_tbl %>%
    select(scenario, Col) %>%
    left_join(rmse_sum, by = "scenario") %>%
    left_join(vr_sum %>% select(scenario, ViolationRate), by = "scenario")
}

res <- summarize_folder(folder, scenarios_tbl = scenarios_tbl)

tab <- res %>%
  pivot_longer(c(RMSE, MAE, ViolationRate), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, ViolationRate = "Violation rate")) %>%
  select(Metric, Col, Value) %>%
  pivot_wider(names_from = Col, values_from = Value) %>%
  arrange(factor(Metric, levels = c("RMSE", "MAE", "Violation rate")))

out_file <- file.path("simul_results/", paste0("tables/", folder_rel, ".xlsx"))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.xlsx(tab, out_file, overwrite = TRUE)





#### SCATTER PLOT

folder_smc <- "NAIVE/ml1e+10_w252_wp1000_re252_sim100000"
folder_smc     <- file.path("simul_results", folder_smc)
files_smc <- list.files(folder_smc, pattern = "\\.rds$", full.names = TRUE)
obj_smc <- lapply(files_smc, readRDS)

folder_naive <- "NAIVE_300/mlNA_w252_wp0252_re252"
folder_naive     <- file.path("simul_results", folder_naive)
files_naive <- list.files(folder_naive, pattern = "\\.rds$", full.names = TRUE)
obj_naive <- lapply(files_naive, readRDS)



rmse_smc  <- bind_rows(lapply(obj_smc, `[[`, "rmse_mae_from_covar"))
covar_smc <- bind_rows(lapply(obj_smc, `[[`, "eval_covar_asset"))

rmse_naive  <- bind_rows(lapply(obj_naive, `[[`, "rmse_mae_from_covar"))
covar_naive <- bind_rows(lapply(obj_naive, `[[`, "eval_covar_asset"))


# subset to your condition
s1 <- rmse_smc  [rmse_smc$cond_asset == 2 & rmse_smc$scenario  == "a0.05b0.025", , drop=FALSE]
s2 <- rmse_naive[rmse_naive$cond_asset == 2 & rmse_naive$scenario == "a0.05b0.025", , drop=FALSE]

# vectors
rmse_A <- s1$RMSE; rmse_B <- s2$RMSE
mae_A  <- s1$MAE;  mae_B  <- s2$MAE

# zoom limits using 0.025 and 0.975
lims_rmse <- as.numeric(quantile(c(rmse_A, rmse_B), probs = c(0.02, 0.98), na.rm = TRUE))
lims_mae  <- as.numeric(quantile(c(mae_A,  mae_B ), probs = c(0.01, 0.99), na.rm = TRUE))

op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1) + 0.1)
on.exit(par(op), add = TRUE)

# --- RMSE panel ---
plot(rmse_A, rmse_B,
     xlim = lims_rmse, ylim = lims_rmse,
     asp  = 1,
     pch  = 19,
     xlab = "SMC",
     ylab = "Naive",
     main = "RMSE")
abline(0, 1, lwd = 1.5, col = gray(0.75), lty = 2)

# --- MAE panel ---
plot(mae_A, mae_B,
     xlim = lims_mae, ylim = lims_mae,
     asp  = 1,
     pch  = 19,
     xlab = "SMC",
     ylab = "Naive",
     main = "MAE")
abline(0, 1, lwd = 1.5, col = gray(0.75), lty = 2)
