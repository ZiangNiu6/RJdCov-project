
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476

# Power computation
f1 <- function(p_value, alpha_list){
  sapply(alpha_list, function(x) length(which(p_value < x)) / length(p_value))
}

# create figure folder
figure_dir <- "simulation/new_figures/low_dimension"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}

################################## plot for varying dimension ##################

# Spherical setting
SP_power <- readRDS("simulation/results/SP_power.rds")

# plot for sphereical setting
SP_rd_power <- apply(SP_power[, "rdhdcov", ], 2, function(x) f1(x, 0.05))
SP_jd_power <- apply(SP_power[, "hdcov", ], 2, function(x) f1(x, 0.05))
SP_mat_power <- apply(SP_power[, "matteson", ], 2, function(x) f1(x, 0.05))
SP_hsic_power <- apply(SP_power[, "dhsic", ], 2, function(x) f1(x, 0.05))
signal_list <- seq(0, 0.25, length.out = 6)

pdf(file = sprintf("%s/sphere_case.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(signal_list, SP_rd_power, type = "o", pch = 19, col = "red",
     main="Spherical", ylab="", xlab = "",ylim = c(0, 1), yaxt="n",xaxt="n")
lines(signal_list, SP_jd_power, type = "o", pch = 19, col = "blue")
lines(signal_list, SP_mat_power, type = "o", pch = 19, col = "purple")
lines(signal_list, SP_hsic_power, type = "o", pch = 19, col = "orange")
abline(h = 0.05, lty = 2, col = "red")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at = signal_list, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# Low-dimensional Gaussian setting
L_power <- readRDS("simulation/results/L_power.rds")

# plot for L
L_rd_power <- apply(L_power[, "rdhdcov", ], 2, function(x) f1(x, 0.05))
L_jd_power <- apply(L_power[, "hdcov", ], 2, function(x) f1(x, 0.05))
L_mat_power <- apply(L_power[, "matteson", ], 2, function(x) f1(x, 0.05))
L_hsic_power <- apply(L_power[, "dhsic", ], 2, function(x) f1(x, 0.05))
signal_list <- seq(0, 0.15, length.out = 6)

pdf(file = sprintf("%s/low_dim_gaussian_case.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(signal_list, L_rd_power, type = "o", pch = 19, col = "red",
     main="Degenerate Gaussian (uniform reference)", ylab="", xlab = "",ylim = c(0, 1), yaxt="n",xaxt="n")
lines(signal_list, L_jd_power, type = "o", pch = 19, col = "blue")
lines(signal_list, L_mat_power, type = "o", pch = 19, col = "purple")
lines(signal_list, L_hsic_power, type = "o", pch = 19, col = "orange")
abline(h = 0.05, lty = 2, col = "red")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at = signal_list, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()


# Low-dimensional Gaussian setting
L_power <- readRDS("simulation/results/L_power_normal_ref.rds")

# plot for L
L_rd_power <- apply(L_power[, "rdhdcov", ], 2, function(x) f1(x, 0.05))
L_jd_power <- apply(L_power[, "hdcov", ], 2, function(x) f1(x, 0.05))
L_mat_power <- apply(L_power[, "matteson", ], 2, function(x) f1(x, 0.05))
L_hsic_power <- apply(L_power[, "dhsic", ], 2, function(x) f1(x, 0.05))
signal_list <- seq(0, 0.15, length.out = 6)

pdf(file = sprintf("%s/low_dim_gaussian_case_normal_ref.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(signal_list, L_rd_power, type = "o", pch = 19, col = "red",
     main="Degenerate Gaussian (normal reference)", ylab="", xlab = "",ylim = c(0, 1), yaxt="n",xaxt="n")
lines(signal_list, L_jd_power, type = "o", pch = 19, col = "blue")
lines(signal_list, L_mat_power, type = "o", pch = 19, col = "purple")
lines(signal_list, L_hsic_power, type = "o", pch = 19, col = "orange")
abline(h = 0.05, lty = 2, col = "red")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at = signal_list, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# Low-dimensional mutlivariate T setting
T_power <- readRDS("simulation/results/T_power.rds")

# plot for L
T_rd_power <- apply(T_power[, "rdhdcov", ], 2, function(x) f1(x, 0.05))
T_jd_power <- apply(T_power[, "hdcov", ], 2, function(x) f1(x, 0.05))
T_mat_power <- apply(T_power[, "matteson", ], 2, function(x) f1(x, 0.05))
T_hsic_power <- apply(T_power[, "dhsic", ], 2, function(x) f1(x, 0.05))
signal_list <- seq(0, 0.15, length.out = 6)

pdf(file = sprintf("%s/low_dim_t.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(signal_list, T_rd_power, type = "o", pch = 19, col = "red",
     main="Degenerate Multivariate t-distribution", ylab="", xlab = "",ylim = c(0, 1), yaxt="n",xaxt="n")
lines(signal_list, T_jd_power, type = "o", pch = 19, col = "blue")
lines(signal_list, T_mat_power, type = "o", pch = 19, col = "purple")
lines(signal_list, T_hsic_power, type = "o", pch = 19, col = "orange")
abline(h = 0.05, lty = 2, col = "red")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at = signal_list, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()
