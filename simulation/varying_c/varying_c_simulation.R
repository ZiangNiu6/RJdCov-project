########Power Analysis################
#######################################
source("simulation/Test_function.R")
n <- 500
B <- 500
group <- c(1,1,1,2,2,2,3,3,3)
c_list <- c(0.2, 1, 5)
names(c_list) <- c("low", "medium", "high")

###Gaussian### vary rho within \{0,0.05,0.1,0.15,0.2,0.25\}
set.seed(1)
reference_data <- matrix(rnorm(n*9), nrow = n, ncol = 9)
emprical <- data.frame(
  low = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, c = 0.2),
  medium = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, c = 1),
  high = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, c = 5)
)

signal_list <- seq(0, 0.25, length.out = 6)
power_GT <- array(NA, dim = c(B, length(c_list), length(signal_list)),
                  dimnames = list(
                    reps = 1:B,
                    c = c_list,
                    signal = signal_list
                  ))

# loop over signal 
for (signal in signal_list) {
  for (b in 1:B) {
    X <- as.matrix(rmvnorm(n, mean = rep(0,9), sigma = sig_T(signal, 9)))
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    
    for(size_c in names(c_list)){
      c = c_list[size_c]
      statis <- computestatisticjdcov(X, dim_list = rep(3, 3), c = c)
      reference_data <- emprical[, size_c]
      power_GT[b, as.character(c), as.character(signal)] <- length(which(reference_data >= statis)) / (length(reference_data) + 1) 
    }
    
    print(power_GT[b, , as.character(signal)])
  }
  print("finish one round")
}

# save the simulation results
saveRDS(power_GT, "simulation/results/GT_power_varying_c.rds")



#######heavy-tailed regression#########
set.seed(1)
reference_data <- matrix(rnorm(n*9), nrow = n, ncol = 9)
emprical <- data.frame(
  low = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, c = 0.2),
  medium = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, c = 1),
  high = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000, c = 5)
)
signal_list <- seq(0, 0.1, length.out = 6)
power_H <- array(NA, dim = c(B, length(c_list), length(signal_list)),
                 dimnames = list(
                   reps = 1:B,
                   c = c_list,
                   signal = signal_list
                 ))
for (signal in signal_list) {
  for (b in 1:B) {
    X <- matrix(rcauchy(500*9, rep(0, 500*9), rep(1, 500*9)), nrow = 500, ncol = 9) + signal*rcauchy(500)
    X_list <- list(X[,1:3], X[, 4:6], X[, 7:9])
    
    for(size_c in names(c_list)){
      c = c_list[size_c]
      statis <- computestatisticjdcov(X, dim_list = rep(3, 3), c = c)
      reference_data <- emprical[, size_c]
      power_H[b, as.character(c), as.character(signal)] <- length(which(reference_data >= statis)) / (length(reference_data) + 1) 
    }
    
    print(power_H[b, , as.character(signal)])
  }
  print("finish one round")
}


saveRDS(power_H, "simulation/results/H_power_varying_c.rds")

