########Power Analysis for low-dimension supported distribution################
source("simulation/Test_function.R")
#######banded covariance matrix setting########
n <- 500
B <- 500
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
group <- c(1,1,1,2,2,2,3,3,3)
method_list <- c("rdhdcov", "hdcov", "matteson", "dhsic")
signal_list <- seq(0, 0.25, length.out = 6)
power_sphere <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                      dimnames = list(
                        reps = 1:B,
                        method = method_list,
                        signal = signal_list
                        ))
for (signal in signal_list) {
  for (b in 1:B) {
    X <- generate_sphere_sample(n = n, d = length(group), num_variable_per_group = 3)
    X[, 4:6] <- X[, 4:6] + signal * X[,1:3]
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(3, 3))
    power_sphere[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1) 
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_sphere[b, "hdcov", as.character(signal)] <- statis$p.value
    
    ## matteson's test
    power_sphere[b, "matteson", as.character(signal)] <- boot_matteson(X, B = 500, group = group)
    
    ## dHSIC test
    power_sphere[b, "dhsic", as.character(signal)] <- dhsic.test(X_list, B = 500)$p.value
    
    print(power_sphere[b, , as.character(signal)])
  }
  print("finish one round")
}

# save the simulation results
saveRDS(power_sphere, "simulation/results/SP_power.rds")

#######low-dimensional Gaussian#########
n <- 500
B <- 500
d <- 10
M <- 3
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(d, 3), niter=1000)
group <- rep(1:3, each = d)
method_list <- c("rdhdcov", "hdcov", "matteson", "dhsic")
signal_list <- seq(0, 0.15, length.out = 6)
power_L <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                 dimnames = list(
                   reps = 1:B,
                   method = method_list,
                   signal = signal_list
                 ))

# generate covariance matrix
cov_mat <- low_rank_covmat(d = d, M = M)

for (signal in signal_list) {
  for (b in 1:B) {
    X <- generate_low_dim_Gaussian(n = n, 
                                   d = d, 
                                   num_variable_per_group = 3, 
                                   cov_mat = cov_mat)
    X[, (d+1):(2*d)] <- signal * X[,1:d] + X[, (d+1):(2*d)] - signal * X[, (2*d + 1):(3*d)]
    X_list <- list(X[,1:d], X[, (d+1):(2*d)], X[, (2*d + 1):(3*d)])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(d, 3))
    power_L[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1) 
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_L[b, "hdcov", as.character(signal)] <- statis$p.value
    
    ## matteson's test
    power_L[b, "matteson", as.character(signal)] <- boot_matteson(X, B = 500, group = group)
    
    ## dHSIC test
    power_L[b, "dhsic", as.character(signal)] <- dhsic.test(X_list, B = 500)$p.value
    
    print(power_L[b, , as.character(signal)])
  }
  print("finish one round")
}


saveRDS(power_L, "simulation/results/L_power.rds")


#######low-dimensional Gaussian project to Gaussian#########
n <- 500
B <- 500
d <- 10
M <- 3
set.seed(1)
reference_data <- matrix(rnorm(n * 3 * d), nrow = n, ncol = 3 * d)
emprical <- gensamdistrjdcov(n, fixgrid = reference_data, 
                             dim_list = rep(d, 3), niter=1000)
group <- rep(1:3, each = d)
method_list <- c("rdhdcov", "hdcov", "matteson", "dhsic")
signal_list <- seq(0, 0.15, length.out = 6)
power_L <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                 dimnames = list(
                   reps = 1:B,
                   method = method_list,
                   signal = signal_list
                 ))

# generate covariance matrix
cov_mat <- low_rank_covmat(d = d, M = M)

for (signal in signal_list) {
  for (b in 1:B) {
    X <- generate_low_dim_Gaussian(n = n, 
                                   d = d, 
                                   num_variable_per_group = 3, 
                                   cov_mat = cov_mat)
    X[, (d+1):(2*d)] <- signal * X[,1:d] + X[, (d+1):(2*d)] - signal * X[, (2*d + 1):(3*d)]
    X_list <- list(X[,1:d], X[, (d+1):(2*d)], X[, (2*d + 1):(3*d)])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(d, 3), 
                                    gridch = reference_data)
    power_L[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1) 
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_L[b, "hdcov", as.character(signal)] <- statis$p.value
    
    ## matteson's test
    power_L[b, "matteson", as.character(signal)] <- boot_matteson(X, B = 500, group = group)
    
    ## dHSIC test
    power_L[b, "dhsic", as.character(signal)] <- dhsic.test(X_list, B = 500)$p.value
    
    print(power_L[b, , as.character(signal)])
  }
  print("finish one round")
}


saveRDS(power_L, "simulation/results/L_power_normal_ref.rds")

#######low-dimensional mutlivariate t-distribution#########
n <- 500
B <- 500
d <- 10
M <- 3
df <- 1
set.seed(1)
emprical <- gensamdistrjdcov(n, dim_list = rep(d, 3), niter=1000)
group <- rep(1:3, each = d)
method_list <- c("rdhdcov", "hdcov", "matteson", "dhsic")
signal_list <- seq(0, 0.15, length.out = 6)
power_T <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                 dimnames = list(
                   reps = 1:B,
                   method = method_list,
                   signal = signal_list
                 ))

# generate covariance matrix
cov_mat <- low_rank_covmat(d = d, M = M)

for (signal in signal_list) {
  for (b in 1:B) {
    X <- generate_low_dim_mvt(n = n, 
                              d = d, 
                              num_variable_per_group = 3, 
                              cov_mat = cov_mat, df = df)
    X[, (d+1):(2*d)] <- signal * X[,1:d] + X[, (d+1):(2*d)] - signal * X[, (2*d + 1):(3*d)]
    X_list <- list(X[,1:d], X[, (d+1):(2*d)], X[, (2*d + 1):(3*d)])
    
    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(d, 3))
    power_T[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1) 
    
    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500, alpha = 0.05)
    power_T[b, "hdcov", as.character(signal)] <- statis$p.value
    
    ## matteson's test
    power_T[b, "matteson", as.character(signal)] <- boot_matteson(X, B = 500, group = group)
    
    ## dHSIC test
    power_T[b, "dhsic", as.character(signal)] <- dhsic.test(X_list, B = 500)$p.value
    
    print(power_T[b, , as.character(signal)])
  }
  print("finish one round")
}


saveRDS(power_T, "simulation/results/T_power.rds")

