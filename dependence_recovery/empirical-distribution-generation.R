# create results directory
dir.create("results")

# generate reference data
set.seed(1)
n_list <- seq(40, 200, by = 40)
source("Test_function.R")

# make a new directory called reference
results_dir <- "results/reference"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  cat("Directory created:", results_dir, "\n")
} else {
  cat("Directory already exists:", results_dir, "\n")
}

for(n in n_list){
  
  empirical_2_2 <- gensamdistrhodcov(n, dim_list = rep(2, 2), niter=1000)
  empirical_2_3 <- gensamdistrhodcov(n, dim_list = rep(2, 3), niter=1000)
  empirical_2_4 <- gensamdistrhodcov(n, dim_list = rep(2, 4), niter=1000)
  
  # save the reference data
  saveRDS(empirical_2_2, sprintf("%s/empirical_%d_2_2.rds", results_dir, n))
  saveRDS(empirical_2_3, sprintf("%s/empirical_%d_2_3.rds", results_dir, n))
  saveRDS(empirical_2_4, sprintf("%s/empirical_%d_2_4.rds", results_dir, n))
}