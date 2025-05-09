library(gmm)
library(stats)
library(multivariance)
data(Finance)

# introduce the textwidth and textheight
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476

################################### pair-wise independence #####################
X_test <- c("real", "JCS", "ZOOM", "food", "manufacture","finance")

# define four industries
real_estate <- c("NHP")
JCS <- c("JCS")
ZOOM <- c("ZOOM")
food <- c("WMK", "GCO")
manufacture <- c("ROG", "MAT", "VOXX")
finance <- c("FNM")

p_rmat <- matrix(0, 6, 6)
colnames(p_rmat) <- X_test
rownames(p_rmat) <- colnames(p_rmat)
p_mat <- matrix(0, 6, 6)
colnames(p_mat) <- colnames(p_rmat)
rownames(p_mat) <- rownames(p_rmat)



# test after averaging
# average over weekday

Monday <- which(weekdays(as.Date(rownames(Finance))) == "Monday")
Friday <- which(weekdays(as.Date(rownames(Finance))) == "Friday")

wkavg_Finance <- matrix(0, nrow = length(Monday), ncol = ncol(Finance))

for (i in 1:length(Monday)) {
  start <- Monday[i]
  end <- Friday[min(which(Friday > Monday[i]))]
  wkavg_Finance[i,] <- apply(Finance[start:end,], 2, mean)
}

# extract the last 500 rows as data matrix
X <- wkavg_Finance[tail(1:nrow(wkavg_Finance), 300), ]
X_real <- X[, which(colnames(Finance) %in% real_estate)]
X_JCS <- X[, which(colnames(Finance) %in% JCS)]
X_ZOOM <- X[, which(colnames(Finance) %in% ZOOM)]
X_food <- X[, which(colnames(Finance) %in% food)]
X_manufacture <- X[, which(colnames(Finance) %in% manufacture)]
X_finance <- X[, which(colnames(Finance) %in% finance)]

# combine the list
X_list <- list(X_real = X_real, X_JCS = X_JCS, X_ZOOM = X_ZOOM, X_food = X_food, X_manufacture = X_manufacture, X_finance = X_finance)

# compute the number of samples
n <- nrow(X)

# loop over different industires
for (j in 1:6) {
  
  # stop the iteration if j is 6
  if(j == 6){
    next
  }
  for (k in (j+1):6) {
    set.seed(1)
    name_1 <- sprintf("X_%s", X_test[j])
    name_2 <- sprintf("X_%s", X_test[k])
    X1 <- as.matrix(X_list[[which(names(X_list) == name_1)]])
    X2 <- as.matrix(X_list[[which(names(X_list) == name_2)]])
    X_t <- cbind(X1, X2)
    X_tlist <- list(X1, X2)
    
    # compute empirical
    empirical <- gensamdistrjdcov(n, dim_list = c(ncol(X1), ncol(X2)), niter=2000)
    
    # rjdcov
    rstat <- computestatisticjdcov(X_t, dim_list = c(ncol(X1), ncol(X2)))
    # jdcov
    stat <- jdcov.test(X_tlist, stat.type = "V", alpha = 0.05, B = 500)
    # pvalue for rjdcov
    p_rmat[j, k] <- length(which(empirical > rstat))/length(empirical)
    # pvalue for jdcov
    p_mat[j, k] <- stat$p.value
    print(p_rmat[j, k])
    print(p_mat[j, k])
  }
}

## adjust p-value with p.adjust
p_value_rjd_p <- data.frame(p_value = numeric(15), 
                            variable_1 = numeric(15), 
                            variable_2 = numeric(15))
p_value_jd_p <- data.frame(p_value = numeric(15), 
                           variable_1 = numeric(15), 
                           variable_2 = numeric(15))
p_value_rjd_p$p_value <- c(p_rmat[1,2:6], p_rmat[2,3:6], p_rmat[3,4:6], p_rmat[4,5:6], p_rmat[5,6])
p_value_jd_p$p_value <- c(p_mat[1,2:6], p_mat[2,3:6], p_mat[3,4:6], p_mat[4,5:6], p_mat[5,6])
p_value_rjd_p$variable_1 <- c(rep(X_test[1], 5), rep(X_test[2], 4), rep(X_test[3], 3),
                              rep(X_test[4], 2), rep(X_test[5], 1))
p_value_rjd_p$variable_2 <- c(X_test[2:6], X_test[3:6], X_test[4:6], X_test[5:6],
                              X_test[6:6])
p_value_jd_p$variable_1 <- p_value_rjd_p$variable_1
p_value_jd_p$variable_2 <- p_value_rjd_p$variable_2

p_value_jd_p$p_value <- p.adjust(p_value_jd_p$p_value, method = "BH")
p_value_rjd_p$p_value <- p.adjust(p_value_rjd_p$p_value, method = "BH")

# pariwise dep table
pairwise_dep <- matrix(1, 6, 6)
pairwise_dep[1,2:6] <- p_value_rjd_p$p_value[1:5]
pairwise_dep[2,3:6] <- p_value_rjd_p$p_value[6:9]
pairwise_dep[3,4:6] <- p_value_rjd_p$p_value[9:11]
pairwise_dep[3,4:6] <- p_value_rjd_p$p_value[10:12]
pairwise_dep[4,5:6] <- p_value_rjd_p$p_value[13:14]
pairwise_dep[5,6] <- p_value_rjd_p$p_value[15]

pairwise_dep[2:6, 1] <-  round(p_value_jd_p$p_value[1:5], 3)
pairwise_dep[3:6, 2] <-  round(p_value_jd_p$p_value[6:9], 3)
pairwise_dep[4:6, 3] <-  round(p_value_jd_p$p_value[10:12], 3)
pairwise_dep[5:6, 4] <-  round(p_value_jd_p$p_value[13:14], 3)
pairwise_dep[6, 5] <-  round(p_value_jd_p$p_value[15], 3)

colnames(pairwise_dep) <- c("Real estate", "JCS", "ZOOM", "Food/Retail", "Manufacture", "Finance")
rownames(pairwise_dep) <- colnames(pairwise_dep)

# output the pairwise_dep
print(round(pairwise_dep, 3))

########################## test third order dependence #########################
JZR <- c("JCS", "ZOOM", "real")
JZF <- c("JCS", "ZOOM", "finance")
JZM <- c("JCS", "ZOOM", "manufacture")

################################## JZR #########################################
set.seed(1)
name_1 <- sprintf("X_%s", JZR[1])
name_2 <- sprintf("X_%s", JZR[2])
name_3 <- sprintf("X_%s", JZR[3])
X1 <- as.matrix(X_list[[which(names(X_list) == name_1)]])
X2 <- as.matrix(X_list[[which(names(X_list) == name_2)]])
X3 <- as.matrix(X_list[[which(names(X_list) == name_3)]])
X_t <- cbind(X1, X2, X3)
X_tlist <- list(X1, X2, X3)

# test with dcov
# bootstrap B = 500
B <- 500
boot_sample_list <- function(X, n, d){
  boot_idx <- matrix(0, n, d)
  resample_X <- list()
  for (i in 1:d) {
    boot_idx[, i] <- sample(1:n, n, replace = TRUE)
    resample_X[[i]] <- as.matrix(X[[i]][boot_idx[, i], ])
  }
  resample_X
}

boot_hodcov <- function(X, B){
  test_stat <- hodcov(X)
  boot <- numeric(B)
  n <- nrow(X[[1]])
  for (b in 1:B) {
    boot[b] <- hodcov(boot_sample_list(X, n, d = length(X)))
  }
  length(which(boot >= test_stat)) / (B + 1)
}

# test with hodcov
p.value_d_JZR <- boot_hodcov(X_tlist, B)
print(sprintf("JZR with d: %.3f", round(p.value_d_JZR, 3)))

# test with jdcov
JZR.test <- jdcov.test(X_tlist, stat.type = "V", alpha = 0.05, B = 500)
p.value_jd_JZR <- JZR.test$p.value
print(sprintf("JZR with jd: %.3f", round(p.value_jd_JZR, 3)))

## contradiction of result

# test with rdcov
empirical <- gensamdistrhodcov(n, dim_list = c(ncol(X1), ncol(X2), ncol(X3)), niter=2000)
rdstat <- computestatisticrdcov(X_t, dim_list = c(ncol(X1), ncol(X2), ncol(X3)))
p.value_rd_JZR <- length(which(empirical > rdstat))/length(empirical)
print(sprintf("JZR with rd: %.3f", round(p.value_rd_JZR, 3)))

# test with rjdcov
empirical <- gensamdistrjdcov(n, dim_list = c(ncol(X1), ncol(X2), ncol(X3)), niter=2000)
rjdstat <- computestatisticjdcov(X_t, dim_list = c(ncol(X1), ncol(X2), ncol(X3)))
p.value_rjd_JZR <- length(which(empirical > rjdstat))/length(empirical)
print(sprintf("JZR with rjd: %.3f", round(p.value_rjd_JZR, 3)))

############################### JZF ############################################
set.seed(1)
name_1 <- sprintf("X_%s", JZF[1])
name_2 <- sprintf("X_%s", JZF[2])
name_3 <- sprintf("X_%s", JZF[3])
X1 <- as.matrix(X_list[[which(names(X_list) == name_1)]])
X2 <- as.matrix(X_list[[which(names(X_list) == name_2)]])
X3 <- as.matrix(X_list[[which(names(X_list) == name_3)]])
X_t <- cbind(X1, X2, X3)
X_tlist <- list(X1, X2, X3)


# test with dcov
# bootstrap B = 500
B <- 500
p.value_d_JZF <- boot_hodcov(X_tlist, B)
round(p.value_d_JZF, 3)

# test with jdcov
JZF.test <- jdcov.test(X_tlist, stat.type = "V", alpha = 0.05, B = 500)
p.value_jd_JZF <- JZF.test$p.value
round(p.value_jd_JZF, 3)

# test with rdcov
empirical <- gensamdistrhodcov(n, dim_list = c(ncol(X1), ncol(X2), ncol(X3)), niter=2000)
rdstat <- computestatisticrdcov(X_t, dim_list = c(ncol(X1), ncol(X2), ncol(X3)))
p.value_rd_JZF <- length(which(empirical > rdstat))/length(empirical)
print(sprintf("JZF with rd: %.3f", round(p.value_rd_JZF, 3)))

# test with rjdcov
empirical <- gensamdistrjdcov(n, dim_list = c(ncol(X1), ncol(X2), ncol(X3)), niter=2000)
rjdstat <- computestatisticjdcov(X_t, dim_list = c(ncol(X1), ncol(X2), ncol(X3)))
p.value_rjd_JZF <- length(which(empirical > rjdstat))/length(empirical)
print(sprintf("JZF with rjd: %.3f", round(p.value_rjd_JZF, 3)))

################################# JZM ##########################################
set.seed(1)
name_1 <- sprintf("X_%s", JZM[1])
name_2 <- sprintf("X_%s", JZM[2])
name_3 <- sprintf("X_%s", JZM[3])
X1 <- as.matrix(X_list[[which(names(X_list) == name_1)]])
X2 <- as.matrix(X_list[[which(names(X_list) == name_2)]])
X3 <- as.matrix(X_list[[which(names(X_list) == name_3)]])
X_t <- cbind(X1, X2, X3)

# test with dcov
# bootstrap B = 500
B <- 500
p.value_d_JZM <- boot_hodcov(X_tlist, B)
round(p.value_d_JZM, 3)

# test with jdcov
JZM.test <- jdcov.test(X_tlist, stat.type = "V", alpha = 0.05, B = 500)
p.value_jd_JZM <- JZM.test$p.value
round(p.value_jd_JZM, 3)


# test with rdcov
empirical <- gensamdistrhodcov(n, dim_list = c(ncol(X1), ncol(X2), ncol(X3)), niter=2000)
rdstat <- computestatisticrdcov(X_t, dim_list = c(ncol(X1), ncol(X2), ncol(X3)))
p.value_rd_JZM <- length(which(empirical > rdstat))/length(empirical)
print(sprintf("JZM with rd: %.3f", round(p.value_rd_JZM, 3)))

# test with rjdcov
empirical <- gensamdistrjdcov(n, dim_list = c(ncol(X1), ncol(X2), ncol(X3)), niter=2000)
rjdstat <- computestatisticjdcov(X_t, dim_list = c(ncol(X1), ncol(X2), ncol(X3)))
p.value_rjd_JZM <- length(which(empirical > rjdstat))/length(empirical)
print(sprintf("JZM with rjd: %.3f", round(p.value_rjd_JZM, 3)))

######################### plot the graph #######################################
fn <- 2000
# order: 1: RE; 2: JCS; 3: Zoom; 4: Food; 5: Manufacture; 6: Finance
############################# Figure for jdcov #################################
X <- matrix(rnorm(fn*6), nrow = fn, ncol = 6)
# create pairwise dependence 
## RE
a11 <- rnorm(fn)
X[,1] <- X[, 1] + a11
X[,4] <- X[, 4] + a11
a12 <- rnorm(fn)
X[,1] <- X[, 1] + a12
X[,5] <- X[,5] + a12
a13 <- rnorm(fn)
X[,1] <- X[, 1] + a13
X[,6] <- X[, 6] + a13

## JCS
a21 <- rnorm(fn)
X[, 2] <- X[, 2] + a21
X[, 4] <- X[, 4] + a21
#jdcov
a22 <- rnorm(fn)
X[, 2] <- X[, 2] + a22
X[, 5] <- X[, 5] + a22


## ZOOM
a31 <- rnorm(fn)
X[, 3] <- X[,3]+a31
X[, 4] <- X[, 4] + a31
# jdcov
a32 <- rnorm(fn)
X[, 3] <- X[,3]+a32
X[, 5] <- X[, 5] + a32
a33 <- rnorm(fn)
X[, 3] <- X[,3]+a33
X[, 6] <- X[, 6] + a33

## F& R
a41 <- rnorm(fn)
X[,4] <- X[,4] + a41
X[,5] <- X[, 5]+ a41
a42 <- rnorm(fn)
X[,4] <- X[,4] + a42
X[,6] <- X[, 6]+ a42

## M
a51 <- rnorm(fn)
X[,5] <- X[,5] + a51
X[,6] <- X[, 6]+ a51

# store the final output
Xr <- data.frame(Real_estate = X[,1], JCS = X[,2], ZOOM = X[,3],
                 Food_Retail = X[, 4], Manufacture = X[,5], Finance = X[, 6])

## use multivariance::dependence.structure using "jdcov"
res_jdcov <- multivariance::dependence.structure(Xr, structure.type = "full")
pdf(file = "real-data/dcov_finance.pdf",   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)
plot(clean.graph(res_jdcov$graph, simplify.pairs = FALSE))
dev.off()

############################# Figure for rjdcov ################################
X <- matrix(rnorm(fn*6), nrow = fn, ncol = 6)
# create pairwise dependence 
## RE & Food; RE & Finance
a11 <- rnorm(fn)
X[,1] <- X[, 1] + a11
X[,4] <- X[, 4] + a11
a13 <- rnorm(fn)
X[,1] <- X[, 1] + a13
X[,6] <- X[, 6] + a13

## JCS & Food
a21 <- rnorm(fn)
X[, 2] <- X[, 2] + a21
X[, 4] <- X[, 4] + a21

## ZOOM & Food
a31 <- rnorm(fn)
X[, 3] <- X[,3]+a31
X[, 4] <- X[, 4] + a31

## Food & Finance
a42 <- rnorm(fn)
X[,4] <- X[,4] + a42
X[,6] <- X[, 6]+ a42

# JZM
sign_list3 <- X[, 2]*X[, 3]*X[, 5]
plus_index3 <- which(sign_list3<=0)
neg_index3 <- which(sign_list3>0)
X[, 5][plus_index3] <- -X[, 5][plus_index3]
X[, 5][neg_index3] <- X[, 5][neg_index3]

# store the final output
Xr <- data.frame(Real_estate = X[,1], JCS = X[,2], ZOOM = X[,3],
                 Food_Retail = X[, 4], Manufacture = X[,5], Finance = X[, 6])

## use multivariance::dependence.structure using "jdcov"
res_rjdcov <- multivariance::dependence.structure(Xr, structure.type = "full")
pdf(file = "real-data/rjdcov_finance.pdf",   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)
plot(clean.graph(res_rjdcov$graph, simplify.pairs = FALSE))
dev.off()
