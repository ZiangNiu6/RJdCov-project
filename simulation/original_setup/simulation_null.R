########Type-I error comparison################
#################################################
source("simulation/Test_function.R")
########Gaussian Setting################

# compute threshold
n <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)

# specify the group of variables
group <- c(1,1,1,2,2,2,3,3,3)
type1_rdhdcov_G <- rep(0, 500)
type1_hdcov_G <- rep(0, 500)
type1_matteson <- rep(0, 500)
type1_dhsic <- rep(0, 500)

# 1000 replications
for (k in 1:500) {
  X <- matrix(0, nrow=n, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X[,a] <- rmvnorm(n, mean = rep(0,3), sigma = diag(3))
  X[,b] <- rmvnorm(n, mean = rep(0,6), sigma = diag(6))
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  type1_rdhdcov_G[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  # JdCov
  statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
  type1_hdcov_G[k] <- statis$p.value
  
  ## matteson's test
  type1_matteson[k] <- boot_matteson(X, B =500, group = group)
  
  ## dHSIC test
  type1_dhsic[k] <- dhsic.test(X_list, B = 500)$p.value
  
}

# store the result
Gaussian_type1 <- data.frame(rdhdcov = type1_rdhdcov_G, hdcov = type1_hdcov_G, matteson = type1_matteson, dHSIC = type1_dhsic)

# save the result
saveRDS(Gaussian_type1, "simulation/results/Gaussian_type1.rds")

############################### compute the rejection rate ####################
# compute the rejection rate
result <- readRDS("simulation/results/Gaussian_type1.rds")
alpha <- 0.05
print("MVN simulation done!")
print(apply(result, 2, function(x) mean(x <= alpha)))

################Copula Setting######################
n <- 500
set.seed(1)
group <- c(1,1,1,2,2,2,3,3,3)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
type1_rdhdcov_CG <- rep(0, 500)
type1_hdcov_CG <- rep(0, 500)
type1_matteson_CG <- rep(0, 500)
type1_dhsic_CG <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(0, nrow=500, ncol = 9)
  a <- c(1, 4, 7)
  b <- c(2,3,5,6,8,9)
  X[,a] <- (rmvnorm(500, mean = rep(0,3), sigma = diag(3)))**3
  X[,b] <- (rmvnorm(500, mean = rep(0,6), sigma = diag(6)))**3
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  type1_rdhdcov_CG[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)
  
  # JdCov
  statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
  type1_hdcov_CG[k] <- statis$p.value
  
  ## matteson's test
  type1_matteson_CG[k] <- boot_matteson(X, B =500, group = group)
  
  ## dHSIC test
  type1_dhsic_CG[k] <- dhsic.test(X_list, B = 500)$p.value
}

# sotre the result
CopulaGaussian_type1 <- data.frame(rdhdcov = type1_rdhdcov_CG, hdcov = type1_hdcov_CG, matteson = type1_matteson_CG, dHSIC = type1_dhsic_CG)

# save the result
saveRDS(CopulaGaussian_type1, "simulation/results/CopulaGaussian_type1.rds")

############################### compute the rejection rate ####################
# compute the rejection rate
result <- readRDS("simulation/results/CopulaGaussian_type1.rds")
alpha <- 0.05
print("Copula simulation done!")
print(apply(result, 2, function(x) mean(x <= alpha)))

########Cauchy Setting################
n <- 500
set.seed(1)
group <- c(1,1,1,2,2,2,3,3,3)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
type1_rdhdcov_C <- rep(0, 500)
type1_hdcov_C <- rep(0, 500)
type1_matteson_C <- rep(0, 500)
type1_dhsic_C <- rep(0, 500)
for (k in 1:500) {
  X <- matrix(rcauchy(9*500, 0, 1), nrow=n, ncol=9)
  X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
  
  # RJdCov
  statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
  type1_rdhdcov_C[k] <- length(which(emprical >= statis)) / (length(emprical) + 1)

  # JdCov
  statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
  type1_hdcov_C[k] <- statis$p.value
  
  ## matteson's test
  type1_matteson_C[k] <- boot_matteson(X, B =500, group = group)
  
  ## dHSIC test
  type1_dhsic_C[k] <- dhsic.test(X_list, B = 500)$p.value
}

# store the result
Cauchy_type1 <- data.frame(rdhdcov = type1_rdhdcov_C, hdcov = type1_hdcov_C, matteson = type1_matteson_C, dHSIC = type1_dhsic_C)

# save the result
saveRDS(Cauchy_type1, "simulation/results/Cauchy_type1.rds")

############################### compute the rejection rate ####################
# compute the rejection rate
result <- readRDS("simulation/results/Cauchy_type1.rds")
alpha <- 0.05
print("Cauchy simulation done!")
print(apply(result, 2, function(x) mean(x <= alpha)))
