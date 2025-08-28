
library(here)
source(here("src/R", "config.R"))  


#### VAR FORECASTING
n_assets <- 1:3

out <- readRDS("C:/Users/55419/Documents/Research/project_1/Code/Exploratory/smc_vines/empirical_results/standard_20250827_093916.rds")

data <- list(
  U      = readRDS("data/PIT.rds")[,n_assets],
  mu_fc  = readRDS("data/returns_mean_forecast.rds")[,n_assets+1],# [,-1],  # drop date col
  sig_fc = readRDS("data/returns_vol_forecast.rds")[,n_assets+1],  #[,-1],
  df_fc = readRDS("data/df_fc.rds")[,n_assets+1],#[,-1]
  shape_fc = readRDS("data/shape_fc.rds")[,n_assets+1]#[,-1]
)

U      <- data$U
mu_fc  <- data$mu_fc
sig_fc <- data$sig_fc
df_fc     <- data$df_fc 
shape_fc     <- data$shape_fc 

y_real <- data$y_real

v <- list(label = "tip",  use_tail_informed_prior = TRUE, tip_method = "EmpTC",          
          tip_k = NULL, tip_sd_logit = 0.025, q_flip = 0.2)
cfg <- modifyList(build_cfg(ncol(data$U)), v[ setdiff(names(v),"label") ])
t_train <- cfg$W_predict
skeleton <- make_skeleton_CVM(U[1:t_train, ], trunc_tree = cfg$trunc_tree)
cfg <- add_first_tree_map(cfg, skeleton)


K <- cfg$K
out <- vector("list", K)
if (is.null(cfg$edge_pair)) return(out)
for (e in seq_len(K)) {
  if (cfg$edge_tree[e] != 1L) {
    next
  } else {
    local_idx <- sum(cfg$edge_tree[seq_len(e)] == 1L)
    pair      <- cfg$edge_pair[local_idx, ]
    uv        <- data_up_to_t[, pair, drop = FALSE]
    emps      <- emp_tails_FRAPO(uv, method = cfg$tip_method, k = cfg$tip_k)
    out[[e]]  <- c(mL = as.numeric(emps["L"]), mU = as.numeric(emps["U"]))
  }
}





train = 126
data_up_to_t <- U[,pair]




out_U <- replicate(5033-train, NA)
out_L <- replicate(5033-train, NA)
for (i in 1:(5033-train)){
  uv        <- data_up_to_t[i:(i+train-1), pair, drop = FALSE]
  emps      <- emp_tails_FRAPO(uv, method = cfg$tip_method)
  out_L[i] <- emps["L"]
  out_U[i] <- emps["U"]
}


plot(out_L, type = "l", col="blue")
lines(out_U)


summary(out_U)
summary(out_L)

