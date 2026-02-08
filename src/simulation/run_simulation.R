
library(here)
library(future)
library(future.apply)

set.seed(11111)

plan(multisession, workers = max(1, parallel::detectCores() - 1))

source(here("src/R", "config.R"))
source(here("src/simulation/naive_simulation.R"))
source(here("src/simulation/main_simulation_nonparallel.R"))


mean_len_grid  <- c(1e10, 500)      
#p_extreme_grid <- c(0.00)   

#W_grid     <- c(252)  #AIC 
W_grid     <- c(252, 126, 504)  #SMC 

#aic_refit_every_grid <- c(252, 63) #AIC 
aic_refit_every_grid <- c(1) #SMC 

#W_predict_grid <- c(252) #AIC 
W_predict_grid <- c(1000) #SMC 

base_dir <- "simul_results/SMC"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)


ml <- mean_len <- 1e10
wp <- 756
re <- 1
w <- W <- 252

for (ml in mean_len_grid) {
  #for (pe in p_extreme_grid) {
  for (w in W_grid) {
    for (wp in W_predict_grid) {
      for (re in aic_refit_every_grid) {
        
        cfg <- modifyList(build_cfg(d = 2), list(M = 1000, label = "M1000", W_predict=wp, aic_refit_every = re, W=w,  use_tail_informed_prior = TRUE))
        
        
        n_train <- cfg$W_predict
        d       <- 2
        n_test  <- 2000
        n       <- n_train+n_test
        
        
        
        tag <- sprintf(
          "ml%s_w%03d_wp%04d_re%03d",
          formatC(ml, format = "d"),
          w,
          wp,
          re
        )
        
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
    }
  }
  #}
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

folder_rel <- "SMC/ml500_w126_wp1000_re001"
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



