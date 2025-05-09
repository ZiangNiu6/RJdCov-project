# ICA simulation
# create figure folder
results_dir <- "ICA/results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  cat("Directory created:", results_dir, "\n")
} else {
  cat("Directory already exists:", results_dir, "\n")
}

# start the simulation
d <- 4
n <- 500
B <- 200
method_list <- c("proposed", "MT")
set.seed(2022)
emprical <- gensamdistrjdcov(n, dim_list = rep(1, d), niter=1000)

# generated random square mixing matrix
M <- mixmat(d)
distribution_parameter <- tibble(
  mixture = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  hyperparameter = list(
    list(mean = 0, var = 1, power = 3),
    list(mean = 0, var = 1, power = 5),
    list(shape = 5),
    list(shape = 10),
    list(mixprob = c(0.3, 0.7), rate = c(1, 5)),
    list(mixprob = c(0.3, 0.7), mean = c(-2, 2), var = c(1, 1), power = 1),
    list(lower = 0, upper = 1),
    list(mixprob = c(0.7, 0.3), mean = c(-2, 2), var = c(3, 1), power = 1),
    list(mixprob = c(0.5, 0.5), mean = c(-2, 2), var = c(2, 2), power = 1),
    list(mixprob = c(0.5, 0.5), mean = c(-2, 2), var = c(2, 2), power = 3),
    list(mixprob = c(0.5, 0.5), mean = c(-2, 2), var = c(2, 2), power = 5),
    list(mixprob = c(0.5, 0.5), mean = c(-2, 2), var = c(2, 2), power = 7)
  ),
  type = c("normal", "normal", "gamma", "gamma", "mix_exp", "mix_gauss", "unif", "mix_gauss", 
           "mix_gauss", "mix_gauss", "mix_gauss", "mix_gauss"),
  result = vector("list", 12)
)

# create a function generating the data accorinding to the information in distribution_parameter
generate_data <- function(type, hyperparameter, mixture, n, d){
  
  # generate data based on type and mixture
  # generate data based on type and mixture
  data <- if (mixture) {
    switch(type,
           mix_exp = {
             # extract information from hyperparameter
             mixprob <- hyperparameter$mixprob
             rate <- hyperparameter$rate
             mixtools::rexpmix(n * d, lambda = mixprob, rate = rate)
           },
           mix_gauss = {
             # extract information from hyperparameter
             mixprob <- hyperparameter$mixprob
             mean <- hyperparameter$mean
             variance <- hyperparameter$var  # renamed from var to avoid conflicts
             power <- hyperparameter$power
             (mixtools::rnormmix(n * d, lambda = mixprob, mu = mean, sigma = sqrt(variance)))^power
           }
    )
  } else {
    switch(type,
           normal = {
             # extract information from hyperparameter
             mean <- hyperparameter$mean
             variance <- hyperparameter$var  # renamed from var to variance
             power <- hyperparameter$power
             (rnorm(n * d, mean = mean, sd = sqrt(variance)))^power
           },
           gamma = {
             # extract information from hyperparameter
             shape <- hyperparameter$shape
             rgamma(n * d, shape = shape)
           },
           unif = {
             # extract information from hyperparameter
             lower <- hyperparameter$lower
             upper <- hyperparameter$upper
             runif(n * d, min = lower, max = upper)
           }
    )
  }
  
  # transform to matrix
  return(matrix(data, nrow = n, ncol = d))
}

# loop over distribution parameter
for (dgp_id in 1:nrow(distribution_parameter)) {
  
  # create the empty matrix
  result_mat <- matrix(NA, nrow = B, ncol = length(method_list),
                       dimnames = list(
                         reps = 1:B,
                         method = method_list
                       ))
  
  # extract necessary information from distribution_parameter
  hyperparameter <- distribution_parameter$hyperparameter[[dgp_id]]
  mixture <- distribution_parameter$mixture[dgp_id]
  type <- distribution_parameter$type[dgp_id]
  
  # loop over B
  for (b in 1:B) {
    set.seed(b)
    
    # generate data
    Z <- generate_data(type = type, hyperparameter = hyperparameter,
                       mixture = mixture, n = n, d = d)
    
    # transform the Z matrix
    a <- whitener(Z%*%M, n.comp = d)
    
    # Whiteout procedure
    W.true <- solve(M%*%a$whitener)
    X <- Z%*%M
    
    # compute ICA with our method
    suppressMessages(out <- JdCovICA(X = X, whiten=TRUE, n.comp = d, 
                                     verbose = TRUE, PIT = TRUE, maxit = 200, 
                                     eps = 1e-10, alpha.eps = 1e-10))
    print(sprintf("Ours: %.2f", amari.error(t(out$W), W.true)))
    result_mat[b, "proposed"] <- amari.error(t(out$W), W.true)
    
    # compute MT estimator
    suppressMessages(out_matteson <- steadyICA(X = X, whiten=TRUE, n.comp = d, 
                                               verbose = TRUE, 
                                               maxit = 200, 
                                               eps = 1e-10, alpha.eps = 1e-10, PIT = TRUE))
    print(sprintf("MT: %.2f", amari.error(t(out_matteson$W), W.true)))
    result_mat[b, "MT"] <- amari.error(t(out_matteson$W), W.true)
  }
  
  # save the results
  distribution_parameter$result[[dgp_id]] <- result_mat
  print(sprintf("Finish number %d DGP: Ours: %.2f, MT: %.2f", dgp_id,
                mean(result_mat[, "proposed"]), mean(result_mat[, "MT"])))
}

# save the result 
saveRDS(distribution_parameter, sprintf("%s/ICA_results.rds", results_dir))