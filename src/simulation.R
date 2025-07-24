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




# ─────────────────────────────────────────────────────────────────────────────
#  sim_static_cop.R      generic static C-vine generator          2025-07-24
# ─────────────────────────────────────────────────────────────────────────────
#  First-tree edges ∈ { indep(0), gaussian(1), bb1(2) }
#  Deeper-tree edges ∈ { indep(0), gaussian(1) }
#
#  Arguments
#    d          : dimension (≥ 3)
#    N          : sample size
#    p_zero     : P(edge = independence) for *all* trees
#    p_bb1      : extra probability mass for BB1 in the first tree
#    beta_U/L   : Beta(α,β) hyper-pars for λU, λL  (upper / lower tail)
#    rho_lo/hi  : uniform range for Gaussian ρ
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
#  sim_static_cop()             rvinecopulib version                2025-07-24
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
#  sim_static_cop()              rvinecopulib-only                    2025-07-24
# ─────────────────────────────────────────────────────────────────────────────
library(rvinecopulib)

bb1_tail2par <- function(lambdaL, lambdaU) {
  theta <- log(2 - lambdaU) / (-log(lambdaL))
  delta <- log(2) / log(2 - lambdaU)
  c(theta, delta)
}

sim_static_cop <- function(d       = 6,
                           N       = 200,
                           p_zero  = 0.50,
                           p_bb1   = 0.50,              # only in first tree
                           beta_U  = c(2, 2),
                           beta_L  = c(2, 2),
                           rho_lo  = -0.99,
                           rho_hi  =  0.99) {
  
  stopifnot(d >= 3, p_zero >= 0, p_bb1 >= 0,
            p_zero + p_bb1 <= 1)
  
  ## 1️⃣  build a VALID C-vine structure matrix -------------------------------
  #   Row i (1-based) has (d-i+1) identical entries   d-i+1,
  #   followed by zeros ─ exactly the pattern rvinecopulib expects.
  struct <- matrix(0L, d, d)
  for (i in 1:d) {
    val <- d - i + 1
    struct[i, 1:val] <- val
  }
  # example d = 6 →
  # 6 6 6 6 6 6
  # 5 5 5 5 5 0
  # 4 4 4 4 0 0
  # 3 3 3 0 0 0
  # 2 2 0 0 0 0
  # 1 0 0 0 0 0
  
  ## 2️⃣  allocate family & parameter holders ---------------------------------
  fam     <- matrix(0L, d, d)          # 0 indep, 1 Gauss, 7 BB1
  theta1  <- theta2 <- matrix(0, d, d) # ρ  or  θ   |  δ for BB1
  
  ## 3️⃣  choose family / parameters for each lower-triangular position -------
  for (j in 1:(d - 1)) {               # column = tree (j = 1 ⇒ first tree)
    first_tree <- (j == 1)
    
    for (i in (j + 1):d) {
      if (runif(1) < p_zero)
        next                                 # keep independence (fam = 0)
      
      if (first_tree) {
        code <- sample(c(1L, 7L), 1, prob = c(1 - p_bb1, p_bb1))
      } else {
        code <- 1L                           # Gaussian only deeper down
      }
      fam[i, j] <- code
      
      if (code == 1L) {                      # Gaussian
        theta1[i, j] <- runif(1, rho_lo, rho_hi)
        
      } else {                               # BB1
        lambdaU <- rbeta(1, beta_U[1], beta_U[2])
        lambdaL <- rbeta(1, beta_L[1], beta_L[2])
        tp      <- bb1_tail2par(lambdaL, lambdaU)
        theta1[i, j] <- tp[1]                # θ
        theta2[i, j] <- tp[2]                # δ
      }
    }
  }
  
  ## 4️⃣  convert to rvinecopulib objects -------------------------------------
  pcs <- vector("list", d - 1)
  for (tr in 1:(d - 1)) {
    pcs[[tr]] <- vector("list", d - tr)
    for (ed in 1:(d - tr)) {
      
      i <- tr + ed     # row index in (fam, theta) matrices
      j <- tr          # column / tree level
      
      f <- fam[i, j]
      pcs[[tr]][[ed]] <-
        if (f == 0L) {
          bicop_dist("indep")
        } else if (f == 1L) {
          bicop_dist("gaussian", parameters = theta1[i, j])
        } else {        # BB1  (code 7)
          bicop_dist("bb1", parameters = c(theta1[i, j], theta2[i, j]))
        }
    }
  }
  
  vc <- vinecop_dist(pcs, structure = struct)
  
  ## 5️⃣  simulate data --------------------------------------------------------
  U <- rvinecop(N, vc)
  colnames(U) <- paste0("U", 1:d)
  
  ## 6️⃣  return --------------------------------------------------------------
  list(U      = U,
       RVM    = vc,          # vinecop_dist object
       family = fam,
       theta  = theta1,
       theta2 = theta2)
}


