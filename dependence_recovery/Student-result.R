noise_level <- seq(0, 1.5, length.out = 4)
num_diff <- matrix(0, nrow = length(noise_level), ncol = 5)
colnames(num_diff) <- as.character(c(40, 80, 120, 160, 200))
rownames(num_diff) <- as.character(noise_level)
consistency_indication <- matrix(0, nrow = length(noise_level), ncol = 5)
colnames(consistency_indication) <- as.character(c(40, 80, 120, 160, 200))
rownames(consistency_indication) <- as.character(noise_level)
fdr <- matrix(0, nrow = length(noise_level), ncol = 5)
colnames(fdr) <- as.character(c(40, 80, 120, 160, 200))
rownames(fdr) <- as.character(noise_level)
for (r in 1:length(noise_level)) {
  # compute the power with result-Student
  results_dir <- "dependence_recovery/results/Student"
  result_Student <- readRDS(sprintf("%s/result-Student-noise%.1f.rds", results_dir, noise_level[r]))
  
  # power result: proportion rate
  diff_set <- matrix(0, nrow = 200, ncol = 5)
  consistency <- matrix(0, nrow = 200, ncol = 5)
  false_discovery <- matrix(0, nrow = 200, ncol = 5)
  for (i in 1:5) {
    result <- result_Student[[i]]
    for (j in 1:200) {
      ommited <- 5 - (result[[j]][[1]][6] + sum(result[[j]][[2]][c(1, 3, 4)]) + result[[j]][[3]][1])
      redundant <- (sum(result[[j]][[1]][1:5]) + result[[j]][[2]][2])
      diff_set[j, i] <- ommited + redundant
      
      # if ommited is zero, then there is consistency
      if(ommited == 0){
        consistency[j, i] <- 1
      }
      
      # compute FDR
      if(sum(unlist(result[[j]])) == 0){
        false_discovery[j, i] <- 0
      }else{
        false_discovery[j, i] <- redundant / (sum(unlist(result[[j]])))
      }
    }
  }
  num_diff[r, ] <- apply(diff_set, 2, mean)
  consistency_indication[r, ] <- apply(consistency, 2, mean)
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

saveRDS(num_diff, sprintf("%s/Student_diff_set_noise.rds", summary_dir))
saveRDS(consistency_indication, sprintf("%s/Student_consistency_noise.rds", summary_dir))
saveRDS(fdr, sprintf("%s/Student_fdr_noise.rds", summary_dir))
