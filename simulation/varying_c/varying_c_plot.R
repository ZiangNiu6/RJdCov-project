
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476
alpha_list <- 0.05

# Power computation
f1 <- function(p_value, alpha_list){
  sapply(alpha_list, function(x) length(which(p_value < x)) / length(p_value))
}

# create figure folder
figure_dir <- "simulation/new_figures/varying_c"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}
signal_list <- seq(0, 0.25, length.out = 6)

# GB
GT_power <- readRDS("simulation/results/GT_power_varying_c.rds")

# plot
GT_1_power <- apply(GT_power[, "0.2", ], 2, function(x) f1(x, alpha_list))
GT_2_power <- apply(GT_power[, "1", ], 2, function(x) f1(x, alpha_list))
GT_3_power <- apply(GT_power[, "5", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Gaussian_toeplitz_varying_c.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(signal_list, GT_1_power, type = "o", pch = 19, col = "red",
     main="Gaussian Toeplitz", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(signal_list, GT_2_power, type = "o", pch = 19, col = "blue")
lines(signal_list, GT_3_power, type = "o", pch = 19, col = "purple")
legend("topleft", legend=c("c = 0.2", "c = 1", "c = 5"),
       col=c("red", "blue", "purple"), lty = 1, cex=0.7)
abline(h = alpha_list, lty = 2, col = "red")
axis(side = 2, at=c(seq(0,1,0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(signal_list),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# Heavy tail
H_power <- readRDS("simulation/results/H_power_varying_c.rds")
signal_list <- seq(0, 0.1, 0.02)

# plot
H_1_power <- apply(H_power[, "0.2", ], 2, function(x) f1(x, alpha_list))
H_2_power <- apply(H_power[, "1", ], 2, function(x) f1(x, alpha_list))
H_3_power <- apply(H_power[, "5", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Cauchy_varying_c.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)


plot(signal_list, H_1_power, type = "o", pch = 19, col = "red",
     main="Cauchy", ylab="", xlab = "",ylim = c(0.05, 1), yaxt="n",xaxt="n")
lines(signal_list, H_2_power, type = "o", pch = 19, col = "blue")
lines(signal_list, H_3_power, type = "o", pch = 19, col = "purple")
abline(h = alpha_list, lty = 2, col = "red")
legend("topleft", legend=c("c = 0.2", "c = 1", "c = 5"),
       col=c("red", "blue", "purple"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=c(signal_list),cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()


