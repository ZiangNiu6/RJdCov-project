########Power Analysis — center-outward reference measure####
######################################################################
source("simulation/Test_function.R")
n <- 500
B <- 500
group <- c(1,1,1,2,2,2,3,3,3)
method_list <- c("rdhdcov")

###Gaussian Toeplitz### vary rho within {0.05,0.1,0.15,0.2,0.25}
set.seed(1)
reference_data <- cbind(
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3)
)
emprical <- gensamdistrjdcov(n, fixgrid = reference_data,
                             dim_list = rep(3, 3), niter=1000)
signal_list <- seq(0.05, 0.25, length.out = 5)
power_GT <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                  dimnames = list(
                    reps = 1:B,
                    method = method_list,
                    signal = signal_list
                  ))

# loop over signal
for (signal in signal_list) {
  for (b in 1:B) {
    X <- as.matrix(rmvnorm(n, mean = rep(0,9), sigma = sig_T(signal, 9)))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    # RJdCov with center-outward reference
    statis <- computestatisticjdcov(X, gridch = reference_data,
                                    dim_list = rep(3, 3))
    power_GT[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1)

    print(power_GT[b, , as.character(signal)])
  }
  print("finish one round")
}

# save the simulation results
saveRDS(power_GT, "simulation/results/GT_power_outward_ref.rds")




#######banded covariance matrix setting########
set.seed(1)
reference_data <- cbind(
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3)
)
emprical <- gensamdistrjdcov(n, fixgrid = reference_data,
                             dim_list = rep(3, 3), niter=1000)
signal_list <- seq(0.05, 0.25, length.out = 5)
power_GB <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                  dimnames = list(
                    reps = 1:B,
                    method = method_list,
                    signal = signal_list
                  ))
for (signal in signal_list) {
  for (b in 1:B) {
    X <- as.matrix(rmvnorm(n, mean = rep(0,9), sigma = sig_B(signal, 9)))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])

    # RJdCov with center-outward reference
    statis <- computestatisticjdcov(X, gridch = reference_data,
                                    dim_list = rep(3, 3))
    power_GB[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1)

    print(power_GB[b, , as.character(signal)])
  }
  print("finish one round")
}

# save the simulation results
saveRDS(power_GB, "simulation/results/GB_power_outward_ref.rds")




#######heavy-tailed regression#########
set.seed(1)
reference_data <- cbind(
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3)
)
emprical <- gensamdistrjdcov(n, fixgrid = reference_data,
                             dim_list = rep(3, 3), niter=1000)
signal_list <- seq(0.02, 0.1, length.out = 5)
power_H <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                 dimnames = list(
                   reps = 1:B,
                   method = method_list,
                   signal = signal_list
                 ))
for (signal in signal_list) {
  for (b in 1:B) {
    X <- matrix(rcauchy(500*9, rep(0, 500*9), rep(1, 500*9)), nrow = 500, ncol = 9) + signal*rcauchy(500)
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])

    # RJdCov with center-outward reference
    statis <- computestatisticjdcov(X, gridch = reference_data,
                                    dim_list = rep(3, 3))
    power_H[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1)

    print(power_H[b, , as.character(signal)])
  }
  print("finish one round")
}


saveRDS(power_H, "simulation/results/H_power_outward_ref.rds")



### additive sin dependence#####
set.seed(1)
reference_data <- cbind(
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3),
  generate_hallin_grid(n, 3)
)
emprical <- gensamdistrjdcov(n, fixgrid = reference_data,
                             dim_list = rep(3, 3), niter=1000)
signal_list <- seq(0.1, 0.5, length.out = 5)
power_S <- array(NA, dim = c(B, length(method_list), length(signal_list)),
                 dimnames = list(
                   reps = 1:B,
                   method = method_list,
                   signal = signal_list
                 ))
for (signal in signal_list) {
  for (b in 1:B) {
    X <- matrix(rcauchy(500*9, rep(0, 500*9), rep(1, 500*9)), nrow = 500, ncol = 9) + sin(signal*rcauchy(500))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])

    # RJdCov with center-outward reference
    statis <- computestatisticjdcov(X, gridch = reference_data,
                                    dim_list = rep(3, 3))
    power_S[b, "rdhdcov", as.character(signal)] <- length(which(emprical >= statis)) / (length(emprical) + 1)

    print(power_S[b, , as.character(signal)])
  }
  print("finish one round")
}

saveRDS(power_S, "simulation/results/S_power_outward_ref.rds")
