library(rvinecopulib)
library(VineCopula)
library(data.table)
library(tictoc)
library(here)

#set.seed(42)


sim_static_cop_3 = function(N=400){
  d=3
  sim_matrix    <- matrix(c(1,2,3, 0,2,3, 0,0,3), 3, 3, byrow = FALSE)
  family_matrix <- matrix(c(0,0,7, 0,0,1, 0,0,0), 3, 3, byrow = FALSE)
  par_matrix  <- matrix(c(0,0.0, 1.47, 0,0,0.4, 0,0,0), 3, 3, byrow = FALSE)
  par2_matrix  <- matrix(c(0,0.0, 3.1, 0,0,0, 0,0,0), 3, 3, byrow = FALSE)
  #theta_matrix  <- matrix(c(0,0.1,-0.7, 0,0,0.7, 0,0,0), 3, 3, byrow = FALSE)
  RVM <- RVineMatrix(sim_matrix, family = family_matrix,  par    = par_matrix, par2   = par2_matrix)
  U   <- RVineSim(N, RVM);  colnames(U) <- paste0("U", 1:d)
  return(list(U      = U,           # simulated observations
              RVM    = RVM,         # full C-vine object
              family = family_matrix,
              par  = par_matrix,
              par2  = par2_matrix))
}



sim_static_cop_8 = function(N=200){
  d=8
  sim_matrix <- matrix(0, nrow = d, ncol = d)
  for(j in 1:d){
    sim_matrix[j:d, j] <- j:d
  }
  theta_matrix <- matrix(c(
    0, 0.66,  0.70,  0.53,  0.38,  0.70,  0.70,  0.58,
    0, 0.00, -0.34,  0.33, -0.10,  0.25,  0.30,  0.28,
    0, 0.00,  0.00,  0.03,  0.15,  0.20,  0.00,  0.00,
    0, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
    0, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
    0, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
    0, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
    0, 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00
  ), nrow = 8, ncol = 8, byrow = TRUE)
  
  family_matrix <- matrix(c(
    0, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, ncol = 8, byrow = TRUE)
  
  RVM <- RVineMatrix(Matrix = sim_matrix,
                     family = family_matrix,
                     par    = theta_matrix)
  
  U <- RVineSim(N, RVM)
  colnames(U) <- paste0("U", 1:d)
  return(U)
}




sim_ar1_copula_corr_3 <- function(N = 200) {
  library(VineCopula)
  
  ## Desired means (correlations)
  mu_corr <- c(0.0, 0.4, 0.6) # last edge equal zero
  names(mu_corr) <- c("x12", "x13", "x23")
  
  ## AR(1) parameters
  phi    <- c(0.99, 0.99, 0.99)                  # strong persistence
  sigma  <- c(0.00, 0.01, 0.01)                  # smaller variance for x23
  mu_lat <- atanh(mu_corr)                      # mean in latent space
  
  ## Simulate AR(1) latent paths
  x <- matrix(NA_real_, nrow = N, ncol = 3,
              dimnames = list(NULL, names(mu_corr)))
  x[1, ] <- mu_lat  # initialize at mean
  
  for (t in 2:N) {
    x[t, ] <- mu_lat + phi * (x[t - 1, ] - mu_lat) + rnorm(3, mean = 0, sd = sigma)
  }
  
  ## Transform to copula parameter space
  theta <- tanh(x)
  
  ## Vine structure (same as before)
  sim_matrix    <- matrix(c(1,2,3, 0,2,3, 0,0,3), 3, 3, byrow = FALSE)
  family_matrix <- matrix(c(0,1,1, 0,0,1, 0,0,0), 3, 3, byrow = FALSE)
  
  ## Simulate from time-varying copula
  U <- matrix(NA_real_, nrow = N, ncol = 3, dimnames = list(NULL, paste0("U", 1:3)))
  
  for (t in 1:N) {
    theta_mat <- matrix(c(
      0, theta[t, "x12"], theta[t, "x13"],
      0, 0,               theta[t, "x23"],
      0, 0,               0
    ), 3, 3, byrow = FALSE)
    
    RVM <- RVineMatrix(sim_matrix, family = family_matrix, par = theta_mat)
    U[t, ] <- RVineSim(1, RVM)
  }
  
  attr(U, "theta_path") <- theta
  return(U)
}


# U <- sim_ar1_copula_corr_3(N = 1000)
# θ <- attr(U, "theta_path")
# matplot(θ, type = "l", lty = 1, col = 1:3, ylab = "θ(t)", xlab = "time")
# legend("topright", legend = colnames(θ), col = 1:3, lty = 1)



sim_ar1_copula_corr_8 <- function(N = 200) {
  library(VineCopula)
  
  d <- 8
  sim_matrix <- matrix(0, d, d)
  for (j in 1:d) sim_matrix[j:d, j] <- j:d
  
  tree_level <- matrix(0, d, d)
  for (j in 1:(d - 1)) for (i in (j + 1):d) tree_level[i, j] <- j
  
  # ── target means & sds by tree ─────────────────────────────────────────────
  target_corr_means <- c(0.6, 0.3, 0.0, rep(0, d - 3))
  target_corr_sds   <- c(0.04, 0.04, 0.03, rep(0.0, d - 3))
  phi               <- 0.98                        # one scalar is enough
  
  theta_array  <- array(NA_real_, c(N, d, d))
  
  for (i in 2:d) {
    for (j in 1:(i - 1)) {
      tr        <- tree_level[i, j]
      mu_lat    <- atanh(target_corr_means[tr])
      sigma_lat <- target_corr_sds[tr]
      
      x <- numeric(N)
      x[1] <- mu_lat
      for (t in 2:N)
        x[t] <- mu_lat + phi * (x[t - 1] - mu_lat) + rnorm(1, 0, sigma_lat)
      
      theta_array[, i, j] <- tanh(x)
    }
  }
  
  # Set one edge from tree 2 to 0
  theta_array[,3,2] <- 0
  
  family_matrix <- matrix(0, d, d)
  family_matrix[lower.tri(family_matrix)] <- 1   # Gaussian everywhere
  
  U <- matrix(NA_real_, N, d)
  for (t in 1:N) {
    theta_t <- matrix(0, d, d)
    theta_t[lower.tri(theta_t)] <-
      theta_array[t, , ][lower.tri(theta_array[t, , ])]
    
    RVM <- RVineMatrix(sim_matrix, family_matrix, theta_t)
    U[t, ] <- RVineSim(1, RVM)
  }
  
  colnames(U) <- paste0("U", 1:d)
  attr(U, "theta_array") <- theta_array
  attr(U, "tree_level")  <- tree_level
  U
}

plot_tree_correlations <- function(U, tree = 1) {
  theta_array <- attr(U, "theta_array")
  tree_level  <- attr(U, "tree_level")
  
  if (is.null(theta_array) || is.null(tree_level))
    stop("The object you passed has no 'theta_array' or 'tree_level' attribute.")
  
  idx <- which(tree_level == tree, arr.ind = TRUE)
  if (nrow(idx) == 0)
    stop(paste("No edges found in tree", tree))
  
  corrs <- sapply(seq_len(nrow(idx)), function(k) {
    i <- idx[k, 1]; j <- idx[k, 2]
    theta_array[, i, j]
  })
  
  edge_names <- apply(idx, 1, \(p) paste0("U", p[2], "-U", p[1]))
  matplot(corrs, type = "l", lty = 1, col = seq_len(ncol(corrs)),
          ylab = "Correlation", xlab = "Time",
          main = paste("Tree", tree, "pairwise correlations"))
  legend("topright", legend = edge_names, lty = 1,
         col = seq_len(ncol(corrs)), bty = "n")
}

get_edge_path <- function(U, tree, edge) {
  theta_array <- attr(U, "theta_array")
  i <- tree + edge       # because i = j + k
  j <- tree
  theta_array[, i, j]    # time-series vector (length N)
}

edges_in_tree <- function(tree_level_mat, tree) {
  which(tree_level_mat == tree, arr.ind = TRUE)
}


#(example: tree 2, edge 4 → pair (i = 6, j = 2))
# U <- sim_ar1_copula_corr_8(500)   # your simulator
# plot_tree_correlations(U, tree = 2)
# path <- get_edge_path(U, tree = 1, edge = 2)   # (i=6, j=2)
# plot(path, type = "l")
# edges_t2 <- edges_in_tree(attr(U, "tree_level"), tree = 2)






sim_static_cop_6 <- function(N      = 200,
                             p_zero = 0.5,
                             rho_lo = -0.99,
                             rho_hi =  0.99) {
  d <- 6
  
  ## 1. C-vine structure matrix (same pattern as your 3-D prototype)
  sim_matrix <- matrix(0, d, d)
  for (j in 1:d)
    sim_matrix[j:d, j] <- j:d         # column j :  j, j+1, …, d
  
  ## 2. Family- and parameter matrices (lower-triangular part only)
  family_matrix <- matrix(0, d, d)
  theta_matrix  <- matrix(0, d, d)
  
  for (j in 1:(d - 1)) {              # columns
    for (i in (j + 1):d) {            # rows below the diagonal
      if (runif(1) > p_zero) {        # keep the edge?
        family_matrix[i, j] <- 1      # 1 = Gaussian copula
        theta_matrix[i,  j] <- runif(1, rho_lo, rho_hi)
      }
    }
  }
  
  ## 3. Build the R-vine object and simulate data
  RVM <- RVineMatrix(sim_matrix,
                     family = family_matrix,
                     par    = theta_matrix)
  
  U <- RVineSim(N, RVM)
  colnames(U) <- paste0("U", 1:d)
  
  ## 4. Return both the data and the true specification
  list(U      = U,           # simulated observations
       RVM    = RVM,         # full C-vine object
       family = family_matrix,
       theta  = theta_matrix)
}


# out <- sim_static_cop_6(N = 200, seed = 2)
# print(out$RVM)        # concise summary of structure, families, and θ
# out$family            # 6 × 6 family matrix (0 = independence, 1 = Gaussian)
# out$theta             # 6 × 6 correlation matrix (zeros on skipped edges)
# head(out$U)           # your data

