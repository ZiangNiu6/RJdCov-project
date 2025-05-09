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
