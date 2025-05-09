# first detect if Rtools is installed
library(pkgbuild)
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
# Mdcov_test <- function(X){
#   if (is.list(X)==FALSE){
#     stop("Input should be list")
#   }
#   if(length(X) < 2){
#     stop("Wrong Dimension")
#   }
#   n <- dim(X[[1]])[1]
#   print(n)
#   d <- length(X)
#   ## compute U function
#   U <- matrix(1, nrow = n, ncol = n)
#   total_expec <- rep(0, d)
#   for (k in 1:d) {
#     dist_mat <- matrix(0, nrow = n, ncol = n)
#     for (i in 1:n) {
#       for (j in 1:n) {
#         dist_mat[i,j] <- dist(rbind(X[[k]][i,], X[[k]][j,]))[1]
#       }
#     }
#     total_expec[k] <- sum(dist_mat)/(n**2)
#   }
#   for (i in 1:n) {
#     for (j in 1:n) {
#       for (k in d) {
#         a <- sum(apply((X[[k]][i,]-t(X[[k]])), 2, Norm))/n+sum(apply((X[[k]][j,]-t(X[[k]])), 2, Norm))/n
#         U[i,j] <- U[i,j]*(a - total_expec[k] - dist(rbind(X[[k]][i,], X[[k]][j,]))[1])
#       }
#     }
#   }
#   return(sum(U)/n)
# }



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
    tstat[i] = jdcov.stat(per_list, stat.type = "V")
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
  randcovSTAT = jdcov.stat(data_list, stat.type = "V")
  return(randcovSTAT)
}
