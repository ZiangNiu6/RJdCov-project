
correlation <- seq(0.1, 0.25, by = 0.05)
num_diff <- matrix(0, nrow = length(correlation), ncol = 5)
colnames(num_diff) <- as.character(c(40, 80, 120, 160, 200))
rownames(num_diff) <- as.character(correlation)
consistency_indication <- matrix(0, nrow = length(correlation), ncol = 5)
colnames(consistency_indication) <- as.character(c(40, 80, 120, 160, 200))
rownames(consistency_indication) <- as.character(correlation)
fdr <- matrix(0, nrow = length(correlation), ncol = 5)
colnames(fdr) <- as.character(c(40, 80, 120, 160, 200))
rownames(fdr) <- as.character(correlation)
for (r in 1:length(correlation)) {
  # compute the power with result-Gaussian
  results_dir <- "dependence_recovery/results/Gaussian"
  result_Gaussian <- readRDS(sprintf("%s/result-Gaussian-%.2f.rds",
                                     results_dir,
                                     correlation[r]))
  
  # 1,4 2,4
  # power result: proportion rate
  diff_set <- matrix(0, nrow = 200, ncol = 5)
  for (i in 1:5) {
    result <- result_Gaussian[[i]]
    for (j in 1:200) {
      omiited <- 6 - (result[[j]][[1]][3] + 
                        result[[j]][[1]][5] + 
                        sum(result[[j]][[2]][2:4]) + 
                        result[[j]][[3]])
      redundant <- sum(result[[j]][[1]][1:2]) + result[[j]][[1]][4]
                    + result[[j]][[1]][6] + result[[j]][[2]][1]
      diff_set[j, i] <- omiited + redundant
    }
  }
  num_diff[r, ] <- apply(diff_set, 2, mean)
  
  # power result: consistency result
  consistency <- matrix(0, nrow = 200, ncol = 5)
  for (i in 1:5) {
    result <- result_Gaussian[[i]]
    for (j in 1:200) {
      if((result[[j]][[1]][3] + 
          result[[j]][[1]][5] + 
          sum(result[[j]][[2]][2:4]) + 
          result[[j]][[3]]) == 6){
        consistency[j, i] <- 1
      }
    }
  }
  
  consistency_indication[r, ] <- apply(consistency, 2, mean)
  
  # false dicovery event: FDE
  false_discovery <- matrix(0, nrow = 200, ncol = 5)
  for (i in 1:5) {
    result <- result_Gaussian[[i]]
    for (j in 1:200) {
      if(sum(unlist(result[[j]])) == 0){
        false_discovery[j, i] <- 0
      }else{
        false_discovery[j, i] <- (sum(result[[j]][[1]][1:2]) + result[[j]][[1]][4]
                                  + result[[j]][[1]][6] + result[[j]][[2]][1]) / sum(unlist(result[[j]]))
      }
    }
  }
  
  fdr[r, ] <- apply(false_discovery, 2, mean)
}

# create summary data folder
summary_dir <- "dependence_recovery/results/summary"
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir)
  cat("Directory created:", summary_dir, "\n")
} else {
  cat("Directory already exists:", summary_dir, "\n")
}

saveRDS(num_diff, sprintf("%s/Gaussian_diff_set.rds", summary_dir))
saveRDS(consistency_indication, sprintf("%s/Gaussian_consistency.rds", summary_dir))
saveRDS(fdr, sprintf("%s/Gaussian_fdr.rds", summary_dir))
