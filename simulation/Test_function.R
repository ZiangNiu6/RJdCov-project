# first detect if Rtools is installed
library(pkgbuild)
library(mvtnorm)
find_rtools(T)

require(clue, quietly=T)
require(energy, quietly = T)
require(randtoolbox, quietly = T)
require(pracma, quietly = T)
#require(HHG, quietly = T)
require(dHSIC, quietly = T)
#require(copula, quietly = T)
require(IndepTest, quietly = T)
require(SpatialNP, quietly = T)
require(EnvStats, quietly = T)
require(LaplacesDemon, quietly = T)
require(MASS, quietly = T)
require(rngWELL, quietly = T)
require(jdcov, quietly = T)
require(steadyICA, quietly = T)



gensamdistrhodcov <- function(N, dim_list, d = length(dim_list), niter=15000,
                              fixgrid = halton(N, sum(dim_list))){
  fixgrid_list <- list()
  for (k in 1:d) {
    if(k==1){
      fixgrid_list[[k]] <- as.matrix(fixgrid[, (1:dim_list[k])])
    }else{
      fixgrid_list[[k]] <- as.matrix(fixgrid[, ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))])  
    }
  }
  tstat = numeric(niter)
  for(i in 1:niter)
  {
    per_list <- list()
    per_list[[1]] <- as.matrix(fixgrid_list[[1]])
    for (k in (1:(d-1))) {
      per_list[[k+1]] <- as.matrix(fixgrid_list[[k+1]][sample(N), ])
    }
    tstat[i] = hodcov(per_list, type = "V")
  }
  return(tstat)
}



gensamdistrjdcov <- function(N, dim_list, d = length(dim_list), niter=15000,
                             fixgrid = halton(N, sum(dim_list)), c = 1){
  fixgrid_list <- list()
  for (k in 1:d) {
    if(k==1){
      fixgrid_list[[k]] <- as.matrix(fixgrid[, (1:dim_list[k])])
    }else{
      fixgrid_list[[k]] <- as.matrix(fixgrid[, ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))])  
    }
  }
  tstat = numeric(niter)
  for(i in 1:niter)
  {
    per_list <- list()
    per_list[[1]] <- as.matrix(fixgrid_list[[1]])
    for (k in (1:(d-1))) {
      per_list[[k+1]] <- as.matrix(fixgrid_list[[k+1]][sample(N), ])
    }
    tstat[i] = jdcov.stat(per_list, stat.type = "V", cc = c)
  }
  return(tstat)
}


computestatisticrdcov <- function(X, dim_list, n = nrow(X), d = length(dim_list),
                                  gridch = halton(n, sum(dim_list)))
{
  data_list <- list()
  X <- as.matrix(X)
  gridch <- as.matrix(gridch)
  for (k in 1:d) {
    distmat <- matrix(0, nrow = n, ncol = n)
    if (k == 1){
      for(i in 1:n){
        distmat[i,] = apply((X[i, (1:dim_list[k])]-t(gridch[, (1:dim_list[k])])), 2, Norm)^2
      }
      assignmentFUN <- solve_LSAP(distmat)
      assignmentSOL <- cbind(seq_along(assignmentFUN), assignmentFUN)
      data_list[[k]] <- as.matrix(gridch[assignmentSOL[,2], (1:dim_list[k])])
    }else{
      for(i in 1:n){
        distmat[i,] = apply((X[i, ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))]
                             -t(gridch[, ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))])), 2, Norm)^2
      }
      assignmentFUN <- solve_LSAP(distmat)
      assignmentSOL <- cbind(seq_along(assignmentFUN), assignmentFUN)
      data_list[[k]] <- as.matrix(gridch[assignmentSOL[,2], ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))])
    }
  }
  randcovSTAT = hodcov(data_list, type = "V")
  return(randcovSTAT)
}

computestatisticjdcov <- function(X, dim_list, n = nrow(X), d = length(dim_list),
                                  gridch = halton(n, sum(dim_list)), c = 1)
{
  data_list <- list()
  X <- as.matrix(X)
  gridch <- as.matrix(gridch)
  for (k in 1:d) {
    distmat <- matrix(0, nrow = n, ncol = n)
    if (k == 1){
      for(i in 1:n){
        distmat[i,] = apply((X[i, (1:dim_list[k])]-t(gridch[, (1:dim_list[k])])), 2, Norm)^2
      }
      assignmentFUN <- solve_LSAP(distmat)
      assignmentSOL <- cbind(seq_along(assignmentFUN), assignmentFUN)
      data_list[[k]] <- as.matrix(gridch[assignmentSOL[,2], (1:dim_list[k])])
    }else{
      for(i in 1:n){
        distmat[i,] = apply((X[i, ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))]
                             -t(gridch[, ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))])), 2, Norm)^2
      }
      assignmentFUN <- solve_LSAP(distmat)
      assignmentSOL <- cbind(seq_along(assignmentFUN), assignmentFUN)
      data_list[[k]] <- as.matrix(gridch[assignmentSOL[,2], ((sum(dim_list[1:k-1])+1):(sum(dim_list[1:k])))])
    }
  }
  randcovSTAT = jdcov.stat(data_list, stat.type = "V", cc = c)
  return(randcovSTAT)
}


### compute the Matteson's estimator using Bootstrap

boot_matteson <- function(X, B, group){
  test_stat <- gmultidcov(X, group = group)
  boot <- numeric(B)
  for (b in 1:B) {
    boot[b] <- gmultidcov(boot_sample(X, n, d = max(group), group), group = group)
  }
  length(which(boot >= test_stat)) / (B + 1)
}


### compute the threshold in higher order distance covariance via bootstrap
boot_hodcov <- function(X, B){
  test_stat <- hodcov(X)
  boot <- numeric(B)
  n <- nrow(X[[1]])
  for (b in 1:B) {
    boot[b] <- hodcov(boot_sample_list(X, n, d = length(X)))
  }
  length(which(boot >= test_stat)) / (B + 1)
}


### Transport Dependency (Nies, Staudt, Munk 2021) — OT-based independence test
# cost = "euclidean": sum of per-group Euclidean distances (Otter.jl default)
#   c(z_i, z_j) = sum_k ||X^(k)_i - X^(k)_j||
# cost = "sq_euclidean": squared Euclidean on the joint space
#   c(z_i, z_j) = ||z_i - z_j||^2 = sum_k ||X^(k)_i - X^(k)_j||^2

# Helper: compute n x n squared Euclidean distance matrix using BLAS-accelerated tcrossprod
sq_euclidean_dist_matrix <- function(A, B) {
  normA <- rowSums(A^2)
  normB <- rowSums(B^2)
  cross <- tcrossprod(A, B)
  D2 <- outer(normA, normB, "+") - 2 * cross
  D2[D2 < 0] <- 0  # clamp floating-point artifacts
  D2
}

# Compute transport dependency statistic for a single dataset
# X_list: list of d matrices, each n x p_j (same interface as jdcov.test, dhsic.test)
# k: number of inner permutations to approximate the product measure
# cost: "euclidean" or "sq_euclidean"
compute_transport_dep <- function(X_list, k = 1, cost = "euclidean") {
  d <- length(X_list)
  n <- nrow(X_list[[1]])

  costs <- numeric(k)
  for (j in 1:k) {
    # Generate independent permutations for groups 2,...,d
    perms <- lapply(2:d, function(g) sample(n))

    # Build cost matrix group-by-group
    C <- matrix(0, n, n)
    for (g in 1:d) {
      A <- X_list[[g]]
      if (g == 1) {
        B <- A  # group 1 is unpermuted in both joint and product
      } else {
        B <- A[perms[[g - 1]], , drop = FALSE]
      }
      D2 <- sq_euclidean_dist_matrix(A, B)
      if (cost == "euclidean") {
        C <- C + sqrt(D2)
      } else {
        C <- C + D2
      }
    }
    assignment <- clue::solve_LSAP(C)
    costs[j] <- sum(C[cbind(1:n, as.integer(assignment))]) / n
  }
  mean(costs)
}

# Permutation test for transport dependency (self-contained, like boot_matteson)
# X_list: list of d matrices, each n x p_j
# B: number of outer permutations (null draws)
# k: number of inner permutations per statistic evaluation
# cost: "euclidean" or "sq_euclidean"
boot_transport_dep <- function(X_list, B = 500, k = 1, cost = "euclidean") {
  d <- length(X_list)
  n <- nrow(X_list[[1]])

  test_stat <- compute_transport_dep(X_list, k = k, cost = cost)

  boot <- numeric(B)
  for (b in 1:B) {
    # Outer permutation: destroy dependence between groups
    X_perm <- list(X_list[[1]])
    for (g in 2:d) {
      X_perm[[g]] <- X_list[[g]][sample(n), , drop = FALSE]
    }
    boot[b] <- compute_transport_dep(X_perm, k = k, cost = cost)
  }
  length(which(boot >= test_stat)) / (B + 1)
}

## bootstrap sample generator with list output
boot_sample_list <- function(X, n, d){
  boot_idx <- matrix(0, n, d)
  resample_X <- list()
  for (i in 1:d) {
    boot_idx[, i] <- sample(1:n, n, replace = TRUE)
    resample_X[[i]] <- X[[i]][boot_idx[, i], ]
  }
  resample_X
}




## bootstrap sample generator 
boot_sample <- function(X, n, d, group){
  boot_idx <- matrix(0, n, d)
  resample_X <- matrix(0, n, length(group))
  for (i in 1:d) {
    boot_idx[, i] <- sample(1:n, n, replace = TRUE)
    resample_X[, which(group == i)] <- X[boot_idx[, i], which(group == i)]
  }
  resample_X
}


# banded matrix
sig_B <- function(rho, d){
  sig_B <- matrix(0, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      if(abs(i-j)<=2){
        if( i != j){
          sig_B[i,j] <- rho
        }else{
          sig_B[i, j] <- 1
        }
      }
    }
  }
  sig_B
}

# Toeplitz matrix
sig_T <- function(rho, d){
  sig_T <- matrix(0, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in 1:d) {
      sig_T[i,j] <- rho**(abs(i-j))
    }
  }
  sig_T
}

# generate sphere distribution
generate_sphere_sample <- function(n, d, num_variable_per_group){
  X <- NULL
  for(group in 1:(d / num_variable_per_group)){
    X_unscaled <- matrix(rnorm(n*num_variable_per_group), nrow = n,
                         ncol = num_variable_per_group)
    row_norms <- sqrt(apply(X_unscaled^2, 1, sum))
    
    # cbind the matrix
    X <- cbind(X, sweep(X_unscaled, 1, row_norms, FUN = "/"))
  }
  return(X)
}

# generate low-rank covariance matrix
low_rank_covmat <- function(d, M){
  # compute the low-rank matrix
  low_rank_mat <- sapply(1:M, function(m){
    runif(d)
  })
  
  # return the low rank matrix
  return(low_rank_mat %*% diag(1, M) %*% t(low_rank_mat))
}

# generate Gaussian data with low-dim rank
generate_low_dim_Gaussian <- function(n, d, num_variable_per_group, cov_mat){
  X <- NULL
  for(group in 1:(d / num_variable_per_group)){
    X_low <- as.matrix(rmvnorm(n, mean = rep(0, d), 
                               sigma = cov_mat))
    
    # cbind the matrix
    X <- cbind(X, X_low)
  }
  return(X)
}

# generate multivariate t-distribution data with low-dim rank
generate_low_dim_mvt <- function(n, d, num_variable_per_group, cov_mat, df){
  X <- NULL
  for(group in 1:(d / num_variable_per_group)){
    X_low <- as.matrix(mvtnorm::rmvt(n = n,
                                     sigma = cov_mat, df = df))

    # cbind the matrix
    X <- cbind(X, X_low)
  }
  return(X)
}

# Generate center-outward (Hallin) reference grid on the unit ball B_p
# Following Hallin (2017) and Shi, Drton & Han (JASA 2022)
# Returns an n x p matrix of grid points inside B_p
generate_outward_grid <- function(n, p) {
  if (p == 1) {
    # Quantiles of Uniform(-1, 1)
    grid <- matrix(2 * (1:n) / (n + 1) - 1, ncol = 1)
    return(grid)
  }

  # p >= 2: concentric shells on B_p
  K <- round(n^(1/p))          # number of shells
  n_per_shell <- floor(n / K)  # base points per shell
  remainder <- n - n_per_shell * K

  grid <- matrix(0, nrow = 0, ncol = p)
  for (j in 1:K) {
    # number of points on this shell (distribute remainder to first shells)
    nj <- n_per_shell + ifelse(j <= remainder, 1, 0)

    # radius: quantile of radial distribution of Uniform(B_p)
    rj <- ((j - 0.5) / K)^(1/p)

    # random unit directions
    Z <- matrix(rnorm(nj * p), nrow = nj, ncol = p)
    row_norms <- sqrt(rowSums(Z^2))
    directions <- Z / row_norms

    # scale by radius
    shell_points <- rj * directions
    grid <- rbind(grid, shell_points)
  }

  return(grid)
}
