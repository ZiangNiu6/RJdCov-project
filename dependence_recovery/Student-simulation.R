set.seed(1)
# set analysis parameters
df <- 2
group <- c(1,1,2,2,3,3,4,4)
r <- 4
N <- c(40, 80, 120, 160, 200)
noise_level <- seq(0, 1.5, length.out = 4)
reference_dir <- "dependence_recovery/results/reference"

# loop over noise 
for (sd_add in noise_level) {
  
  # loop over number of sample
  result <- list()
  for (k in 1:length(N)) {
    n <- N[k]
    result_n <- list()
    
    for (b in 1:200) {
      supersetinclusion <- 0
      X <- matrix(rt(n*8, df), nrow = n, ncol = 8)
      
      # create joint dependence within 123 (first coordinate)
      sign_list <- X[, 1]*X[, 3]*X[, 5]
      plus_index <- which(sign_list<=0)
      neg_index <- which(sign_list>0)
      X[, 5][plus_index] <- -X[, 5][plus_index]
      X[, 5][neg_index] <- X[, 5][neg_index]
      
      # create joint dependence within 34 
      X[, 7] <- abs(X[, 7]) * ifelse(X[, 6] > 0, 1, -1)
      
      # constrct the final data
      data <- X + matrix(rnorm(n * 8, sd = sd_add), nrow = n, ncol = 8)
      
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
      print(TT)
      print(b)
    }
    result[[k]] <- result_n
    print(n)
  }
  
  # create result directory
  results_dir <- "dependence_recovery/results/Student"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
    cat("Directory created:", results_dir, "\n")
  } else {
    cat("Directory already exists:", results_dir, "\n")
  }
  
  # save the result
  saveRDS(result, sprintf("%s/result-Student-noise%.1f.rds", results_dir, sd_add))
}
