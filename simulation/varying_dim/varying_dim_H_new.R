# This is a script for Cauchy setup (n=50, signal=0.9, d=20-100);
set.seed(1)
source("simulation/Test_function.R")

# specify the number of replications and method list
B <- 300
method_list <- c("rdhdcov", "hdcov", "matteson", "dhsic", "tdep")

# generate the data
generating_resample <- function(n, d, num_group = 3, B = 1000){
  resampled_data <- gensamdistrjdcov(n, dim_list = rep(d, num_group), niter = B)
  return(resampled_data)
}

generating_data <- function(n, d, num_group = 3, signal = 0.9){

  # compute how many total samples we need;
  total_sample <- n*d*num_group

  # obtain the combined matrix
  X <- matrix(rcauchy(total_sample, rep(0, total_sample), rep(1, total_sample)),
              nrow = n, ncol = d*num_group) + signal*rcauchy(n)

  # split the matrix to different elements in the list
  X_list <- lapply(seq(1, ncol(X), by = d), function(i) X[, i:(i + d - 1)])

  # return the result
  return(list(data_mat = X,
              data_list = X_list))
}

## Increasing the dimension
n <- 50
num_group <- 3
dim_list <-  seq(20, 100, length.out = 5)
power_H <- array(NA, dim = c(B, length(method_list), length(dim_list)),
                 dimnames = list(
                   reps = 1:B,
                   method = method_list,
                   dimension = dim_list
                 ))
for (d in dim_list) {
  # specify the group and get the resampled data
  group <- rep(1:num_group, each = d)
  resampled_data <- generating_resample(n, d, num_group)
  for (b in 1:B) {

    # generate data
    data <- generating_data(n, d, num_group, signal = 0.9)
    X <- data$data_mat
    X_list <- data$data_list

    # RJdCov
    statis <- computestatisticjdcov(X, dim_list = rep(d, num_group))
    power_H[b, "rdhdcov", as.character(d)] <- length(which(resampled_data >= statis)) / (length(resampled_data) + 1)

    # JdCov
    statis <- jdcov.test(X_list, stat.type = "V", B = 500)
    power_H[b, "hdcov", as.character(d)] <- statis$p.value

    ## matteson's test
    power_H[b, "matteson", as.character(d)] <- boot_matteson(X, B = 500, group = group)

    ## dHSIC test
    power_H[b, "dhsic", as.character(d)] <- dhsic.test(X_list, B = 500)$p.value

    ## Transport Dependency test
    power_H[b, "tdep", as.character(d)] <- boot_transport_dep(X_list, B = 500, k = 1, cost = "euclidean")

    print(power_H[b, , as.character(d)])
  }
  print(sprintf("finish one round with d = %d", d))
}


saveRDS(power_H, "simulation/results/H_power_varying_dim_new.rds")
