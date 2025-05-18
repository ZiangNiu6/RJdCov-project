
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476
alpha_list <- 0.2 / 2^{seq(1, 6, 1)}

# Power computation
source("simulation/Test_function.R")
f1 <- function(p_value, alpha_list){
  sapply(alpha_list, function(x) length(which(p_value < x)) / length(p_value))
}

# create figure folder
figure_dir <- "simulation/new_figures/varying_alpha"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}

# GB
GB_power <- readRDS("simulation/results/GB_power.rds")

# plot
GB_rd_power <- apply(GB_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GB_jd_power <- apply(GB_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
GB_mat_power <- apply(GB_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
GB_hsic_power <- apply(GB_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Gaussian_banded_varying_alpha.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(alpha_list, GB_rd_power[,2], type = "o", pch = 19, col = "red",
     main="Gaussian banded", ylab="", xlab = "",ylim = c(0,1), yaxt="n",xaxt="n",
     log = "x")
lines(alpha_list, GB_jd_power[,2], type = "o", pch = 19, col = "blue")
lines(alpha_list, GB_mat_power[,2], type = "o", pch = 19, col = "purple")
lines(alpha_list, GB_hsic_power[,2], type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(round(alpha_list, 3)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Significance level", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# GT
GT_power <- readRDS("simulation/results/GT_power.rds")

# plot
GT_rd_power <- apply(GT_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GT_jd_power <- apply(GT_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
GT_mat_power <- apply(GT_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
GT_hsic_power <- apply(GT_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Gaussian_topeliz_varying_alpha.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(alpha_list, GT_rd_power[,3], type = "o", pch = 19, col = "red",
     main="Gaussian Toeplitz", ylab="", xlab = "",ylim = c(0,1), yaxt="n",xaxt="n",
     log = "x")
lines(alpha_list, GT_jd_power[,3], type = "o", pch = 19, col = "blue")
lines(alpha_list, GT_mat_power[,3], type = "o", pch = 19, col = "purple")
lines(alpha_list, GT_hsic_power[,3], type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(round(alpha_list, 3)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Significance level", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()


# Heavy tail
H_power <- readRDS("simulation/results/H_power.rds")

# plot
H_rd_power <- apply(H_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
H_jd_power <- apply(H_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
H_mat_power <- apply(H_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
H_hsic_power <- apply(H_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Cauchy_varying_alpha.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(alpha_list, H_rd_power[,4], type = "o", pch = 19, col = "red",
     main="Cauchy", ylab="", xlab = "",ylim = c(0,1), yaxt="n",xaxt="n",
     log = "x")
lines(alpha_list, H_jd_power[,4], type = "o", pch = 19, col = "blue")
lines(alpha_list, H_mat_power[,4], type = "o", pch = 19, col = "purple")
lines(alpha_list, H_hsic_power[,4], type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(round(alpha_list, 3)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Significance level", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()



# Additive sine dependence
S_power <- readRDS("simulation/results/S_power.rds")

# plot
S_rd_power <- apply(S_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
S_jd_power <- apply(S_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
S_mat_power <- apply(S_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
S_hsic_power <- apply(S_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Sin_dependence_varying_alpha.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(alpha_list, S_rd_power[,3], type = "o", pch = 19, col = "red",
     main="Sine dependence", ylab="", xlab = "",ylim = c(0,1), yaxt="n",xaxt="n",
     log = "x")
lines(alpha_list, S_jd_power[,3], type = "o", pch = 19, col = "blue")
lines(alpha_list, S_mat_power[,3], type = "o", pch = 19, col = "purple")
lines(alpha_list, S_hsic_power[,3], type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(round(alpha_list, 3)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Significance level", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()



