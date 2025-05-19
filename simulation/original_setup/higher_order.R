########higher-order dependence########

###generate data######
generate_data_hod <- function(n, d){
  
  # check the type
  data <- list(
    gaussian = matrix(rnorm(n * d), n, d),
    student_3 = matrix(rt(n*d, df = 3), n, d),
    student_2 = matrix(rt(n*d, df = 2), n, d),
    cauchy = matrix(rcauchy(n*d), n, d)
    )
  # return the data
  return(data)
}

# specify the parameters
n <- 500
dist_list <- c("gaussian", "student_3", "student_2", "cauchy")
B <- 500
d <- 9
set.seed(1)
source("simulation/Test_function.R")

## test the pairwise independence
emprical_data <- list(
  emprical_pairwise = gensamdistrhodcov(n, dim_list = rep(3, 2), niter=1000),
  emprical_higherorder = gensamdistrhodcov(n, dim_list = rep(3, 3), niter=1000),
  emprical_full = gensamdistrjdcov(n, dim_list = rep(3, 3), niter=1000)
)

# create empty matrices
result <- list(
  proposed_pairwise = matrix(NA, nrow = B, ncol = length(dist_list),
                             dimnames = list(
                               rep = 1:B,
                               method = dist_list
                             )),
  alternative_pairwise = matrix(NA, nrow = B, ncol = length(dist_list),
                                dimnames = list(
                                  rep = 1:B,
                                  method = dist_list
                                )),
  proposed_higherorder = matrix(NA, nrow = B, ncol = length(dist_list),
                                dimnames = list(
                                  rep = 1:B,
                                  method = dist_list
                                )),
  alternative_higherorder = matrix(NA, nrow = B, ncol = length(dist_list),
                                   dimnames = list(
                                     rep = 1:B,
                                     method = dist_list
                                   )),
  proposed_full = matrix(NA, nrow = B, ncol = length(dist_list),
                         dimnames = list(
                           rep = 1:B,
                           method = dist_list
                         )),
  alternative_full = matrix(NA, nrow = B, ncol = length(dist_list),
                            dimnames = list(
                              rep = 1:B,
                              method = dist_list
                            ))
)

# loop over B
type_list <- c("pairwise", "higherorder", "full")
for (k in 1:B) {
  
  # generate data
  data <- generate_data_hod(n = n, d = d)
  
  # loop over distribution
  for (dist in dist_list) {
    X <- data[[dist]]
    sign_list <- X[, 1]*X[, 4]*X[, 7]
    plus_index <- which(sign_list<=0)
    neg_index <- which(sign_list>0)
    X[, 7][plus_index] <- -X[, 7][plus_index]
    X[, 7][neg_index] <- X[, 7][neg_index]
    
    # loop over type list
    for (type in type_list) {
      
      # extract the data as list
      if(type == "pairwise"){
        X_list <- list(X[,4:6], X[, 7:9])
      }else{
        X_list <- list(X[,1:3], X[,4:6], X[,7:9])
      }
      
      # extract the emprical data and result name
      emprical <- emprical_data[[sprintf("emprical_%s", type)]]
      proposed_name <- paste("proposed", type, sep = "_")
      alternative_name <- paste("alternative", type, sep = "_")
      
      # run the result
      switch (type,
        pairwise = {
          
          # rdhod
          statis_proposed <- computestatisticjdcov(X[,4:9], dim_list = rep(3, 2))
          
          # hod
          statis_alternative <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)$p.value
        },
        higherorder = {
          
          # rdhod
          statis_proposed <- computestatisticrdcov(X, dim_list = rep(3, 3))
          
          # hod
          statis_alternative <- boot_hodcov(X_list, 500)
        },
        full = {
          
          # rdhod
          statis_proposed <- computestatisticjdcov(X, dim_list = rep(3, 3))

          # hod
          statis_alternative <- jdcov.test(X_list, stat.type = "V", alpha = 0.05)$p.value
        }
      )
      
      # store the results 
      result[[proposed_name]][k, dist] <- length(which(emprical >= statis_proposed)) / (length(emprical) + 1)
      result[[alternative_name]][k, dist] <- statis_alternative
      
      # print the result
      print(sprintf("proposed %s_%s: %.2f", type, dist, result[[proposed_name]][k, dist]))
      print(sprintf("alternative %s_%s: %.2f", type, dist, result[[alternative_name]][k, dist])) 
    }
  }
}

# save the result
saveRDS(result, "simulation/results/dependence_results.rds")

############################### compute the rejection rate ####################
# compute the rejection rate
result <- readRDS("simulation/results/dependence_results.rds")
alpha <- 0.05
print(sapply(result, function(x){
  apply(x, 2, function(y) mean(y <= alpha))
})
)