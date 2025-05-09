set.seed(1)
reference_dir <- "dependence_recovery/results/reference"
library(MASS)
correlation <- seq(0.05, 0.25, by = 0.05)
for (aa in 1:length(correlation)) {
  set.seed(1)
  a <- correlation[aa]

  group <- c(1,1,2,2,3,3,4,4)
  sig_T <- matrix(c(1, 0.2, 0, 0, 0, 0, a, a,
                    0.2, 1, 0, 0, 0, 0, a, a,
                    0.0, 0, 1, 0, 0, 0, a, a,
                    0.0, 0, 0, 1, 0, 0, a, a,
                    0, 0, 0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 1, 0, 0,
                    a, a, a, a, 0, 0, 1, 0,
                    a, a, a, a, 0, 0, 0, 1), nrow = 8, ncol = 8)
  result <- list()
  r <- 4
  N <- c(40, 80, 120, 160, 200)
  for (k in 1:length(N)) {
    n <- N[k]
    result_n <- list()
    for (b in 1:200) {
      supersetinclusion <- 0
      data <- as.matrix(mvrnorm(n, mu = rep(0, length(group)), Sigma = sig_T))
      # dependence structure detection algorithm
      S <- list()
      S[[1]] <- matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), nrow = 2, ncol = 6)
      S[[2]] <- matrix(c(1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4), nrow = 3, ncol = 4)
      S[[3]] <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
      TT <- list()
      for (i in 1:(r-1)) {
        S_i <- ncol(S[[i]])
        if(supersetinclusion == 0){
          TT[[i]] <- rep(0, S_i)
        }
        for (j in 1:S_i) {
          if(TT[[i]][j] == 1){
            next
          }
          # specify the group index
          group_index <- S[[i]][, j]
          index <- sort(c(2*group_index-1, 2*group_index))
          
          # extract the test data from the full data matrix
          data_test <- data[, index]
          
          # read the empirical distribution from RDS file
          empirical <- readRDS(sprintf("%s/empirical_%d_%d_%d.rds", 
                                       reference_dir,
                                       n, 2, nrow(S[[i]])))
          
          # compute the test statiusti and p-value the p-value
          test <- computestatisticrdcov(data_test, dim_list = rep(2, nrow(S[[i]])))
          p_value <- 1 - length(which(test >= empirical)) / (length(empirical) + 1)
          if(p_value <= 0.05){
            TT[[i]][j] <- 1
            if(i == r-1){
              break
            }else{
              for (l in (i+1):(r-1)) {
                S_l <- ncol(S[[l]])
                if(tryCatch({length(TT[[l]])}, error= function(e) 1e4) == 1e4)
                {
                  TT[[l]] <- rep(0, S_l)
                }
                for (ll in 1:S_l) {
                  if(all(S[[i]][,j] %in% S[[l]][,ll])){
                    TT[[l]][ll] <- 1
                    supersetinclusion <- 1
                  }
                }
              }
            }
          }
        }
      }
      result_n[[b]] <- TT
      print(b)
    }
    result[[k]] <- result_n
    print(n)
  }
  
  # creat result directory
  results_dir <- "dependence_recovery/results/Gaussian"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
    cat("Directory created:", results_dir, "\n")
  } else {
    cat("Directory already exists:", results_dir, "\n")
  }
  
  saveRDS(result, sprintf("%s/result-Gaussian-%.2f.rds", results_dir, a))
  print(aa)
}