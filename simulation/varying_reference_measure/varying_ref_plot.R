
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476
alpha_list <- 0.05

# Power computation
f1 <- function(p_value, alpha_list){
  sapply(alpha_list, function(x) length(which(p_value < x)) / length(p_value))
}

# create figure folder
figure_dir <- "simulation/new_figures/varying_ref"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}

# GB
GB_power <- readRDS("simulation/results/GB_power_normal_ref.rds")
GB_power_unif <- readRDS("simulation/results/GB_power.rds")
GB_power_outward <- readRDS("simulation/results/GB_power_outward_ref.rds")

# plot
GB_rd_norm_power <- apply(GB_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GB_rd_unif_power <- apply(GB_power_unif[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GB_rd_outward_power <- apply(GB_power_outward[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GB_jd_power <- apply(GB_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
GB_mat_power <- apply(GB_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
GB_hsic_power <- apply(GB_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Gaussian_banded_normal_ref.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(seq(0.05,0.25,0.05), GB_rd_norm_power, type = "o", pch = 19, col = "pink",
     main="Gaussian banded", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(seq(0.05,0.25,0.05), GB_rd_unif_power, type = "o", pch = 19, col = "red")
lines(seq(0.05,0.25,0.05), GB_rd_outward_power, type = "o", pch = 19, col = "darkgreen")
lines(seq(0.05,0.25,0.05), GB_jd_power, type = "o", pch = 19, col = "blue")
lines(seq(0.05,0.25,0.05), GB_mat_power, type = "o", pch = 19, col = "purple")
lines(seq(0.05,0.25,0.05), GB_hsic_power, type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("Proposed (normal ref)", "Proposed (uniform ref)",
                               "Proposed (outward ref)", "JdCov", "MT", "dHSIC"),
       col=c("pink","red", "darkgreen", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0,1,0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(seq(0,0.25,0.05)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# GT
GT_power <- readRDS("simulation/results/GT_power_normal_ref.rds")
GT_power_unif <- readRDS("simulation/results/GT_power.rds")
GT_power_outward <- readRDS("simulation/results/GT_power_outward_ref.rds")

# plot
GT_rd_norm_power <- apply(GT_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GT_rd_unif_power <- apply(GT_power_unif[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GT_rd_outward_power <- apply(GT_power_outward[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
GT_jd_power <- apply(GT_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
GT_mat_power <- apply(GT_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
GT_hsic_power <- apply(GT_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Gaussian_topeliz_normal_ref.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(seq(0.05,0.25,0.05), GT_rd_norm_power, type = "o", pch = 19, col = "pink",
     main="Gaussian Toeplitz", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(seq(0.05,0.25,0.05), GT_rd_unif_power, type = "o", pch = 19, col = "red")
lines(seq(0.05,0.25,0.05), GT_rd_outward_power, type = "o", pch = 19, col = "darkgreen")
lines(seq(0.05,0.25,0.05), GT_jd_power, type = "o", pch = 19, col = "blue")
lines(seq(0.05,0.25,0.05), GT_mat_power, type = "o", pch = 19, col = "purple")
lines(seq(0.05,0.25,0.05), GT_hsic_power, type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("Proposed (normal ref)", "Proposed (uniform ref)",
                               "Proposed (outward ref)", "JdCov", "MT", "dHSIC"),
       col=c("pink", "red", "darkgreen", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis = 0.7, mgp = c(1,0.6,0))
axis(side = 1, at = c(seq(0, 0.25, 0.05)), cex.axis = 0.7, mgp = c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()


# Heavy tail
H_power <- readRDS("simulation/results/H_power_normal_ref.rds")
H_power_unif <- readRDS("simulation/results/H_power.rds")
H_power_outward <- readRDS("simulation/results/H_power_outward_ref.rds")

# plot
H_rd_norm_power <- apply(H_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
H_rd_unif_power <- apply(H_power_unif[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
H_rd_outward_power <- apply(H_power_outward[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
H_jd_power <- apply(H_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
H_mat_power <- apply(H_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
H_hsic_power <- apply(H_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Cauchy_normal_ref.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)


plot(seq(0.02,0.1,0.02), H_rd_norm_power, type = "o", pch = 19, col = "pink",
     main="Cauchy", ylab="", xlab = "",ylim = c(0.05, 1), yaxt="n",xaxt="n")
lines(seq(0.02,0.1,0.02), H_rd_unif_power, type = "o", pch = 19, col = "red")
lines(seq(0.02,0.1,0.02), H_rd_outward_power, type = "o", pch = 19, col = "darkgreen")
lines(seq(0.02,0.1,0.02), H_jd_power, type = "o", pch = 19, col = "blue")
lines(seq(0.02,0.1,0.02), H_mat_power, type = "o", pch = 19, col = "purple")
lines(seq(0.02,0.1,0.02), H_hsic_power, type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("Proposed (normal ref)", "Proposed (uniform ref)",
                               "Proposed (outward ref)", "JdCov", "MT", "dHSIC"),
       col=c("pink", "red", "darkgreen", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(seq(0, 0.1, 0.02)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()



# Additive sine dependence
S_power <- readRDS("simulation/results/S_power_normal_ref.rds")
S_power_unif <- readRDS("simulation/results/S_power.rds")
S_power_outward <- readRDS("simulation/results/S_power_outward_ref.rds")

# plot
S_rd_norm_power <- apply(S_power[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
S_rd_unif_power <- apply(S_power_unif[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
S_rd_outward_power <- apply(S_power_outward[, "rdhdcov", ], 2, function(x) f1(x, alpha_list))
S_jd_power <- apply(S_power[, "hdcov", ], 2, function(x) f1(x, alpha_list))
S_mat_power <- apply(S_power[, "matteson", ], 2, function(x) f1(x, alpha_list))
S_hsic_power <- apply(S_power[, "dhsic", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Sin_dependence_normal_ref.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)


plot(seq(0.1, 0.5, length.out = 5), S_rd_norm_power, type = "o", pch = 19, col = "pink",
     main="Sine dependence", ylab="", xlab = "",ylim = c(0.05, 1), yaxt="n",xaxt="n")
lines(seq(0.1, 0.5, length.out = 5), S_rd_unif_power, type = "o", pch = 19, col = "red")
lines(seq(0.1, 0.5, length.out = 5), S_rd_outward_power, type = "o", pch = 19, col = "darkgreen")
lines(seq(0.1, 0.5, length.out = 5), S_jd_power, type = "o", pch = 19, col = "blue")
lines(seq(0.1, 0.5, length.out = 5), S_mat_power, type = "o", pch = 19, col = "purple")
lines(seq(0.1, 0.5, length.out = 5), S_hsic_power, type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("Proposed (normal ref)", "Proposed (uniform ref)",
                           "Proposed (outward ref)", "JdCov", "MT", "dHSIC"),
       col=c("pink", "red", "darkgreen", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(seq(0, 0.5, 0.1)),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)


dev.off()


