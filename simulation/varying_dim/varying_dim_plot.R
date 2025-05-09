
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476
alpha_list <- 0.2 / 2^{seq(1, 6, 1)}

# Power computation
f1 <- function(p_value, alpha_list){
  sapply(alpha_list, function(x) length(which(p_value < x)) / length(p_value))
}

# create figure folder
figure_dir <- "simulation/new_figures/varying_dim"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}

################################## plot for varying dimension ##################

# GB
GB_power <- readRDS("simulation/results/GB_power_varying_dim.rds")

# plot for GB
GB_rd_power <- apply(GB_power[, "rdhdcov", ], 2, function(x) f1(x, 0.05))
GB_jd_power <- apply(GB_power[, "hdcov", ], 2, function(x) f1(x, 0.05))
GB_mat_power <- apply(GB_power[, "matteson", ], 2, function(x) f1(x, 0.05))
GB_hsic_power <- apply(GB_power[, "dhsic", ], 2, function(x) f1(x, 0.05))
dim_list <- as.numeric(names(GB_rd_power))

pdf(file = sprintf("%s/Gaussian_banded_varying_dim.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(dim_list, GB_rd_power, type = "o", pch = 19, col = "red",
     main="Gaussian banded", ylab="", xlab = "",ylim = c(0,0.7), yaxt="n",xaxt="n")
lines(dim_list, GB_jd_power, type = "o", pch = 19, col = "blue")
lines(dim_list, GB_mat_power, type = "o", pch = 19, col = "purple")
lines(dim_list, GB_hsic_power, type = "o", pch = 19, col = "orange")
legend("topright", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 0.75, 0.15)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at = dim_list, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Dimension", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# Heavy tail
H_power <- readRDS("simulation/results/H_power_varying_dim.rds")

# plot for H
H_rd_power <- apply(H_power[, "rdhdcov", ], 2, function(x) f1(x, 0.05))
H_jd_power <- apply(H_power[, "hdcov", ], 2, function(x) f1(x, 0.05))
H_mat_power <- apply(H_power[, "matteson", ], 2, function(x) f1(x, 0.05))
H_hsic_power <- apply(H_power[, "dhsic", ], 2, function(x) f1(x, 0.05))
dim_list <- as.numeric(names(H_rd_power))

pdf(file = sprintf("%s/Cauchy_varying_dim.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.45*TEXTHEIGHT)

plot(dim_list, H_rd_power, type = "o", pch = 19, col = "red",
     main="Cauchy", ylab="", xlab = "",ylim = c(0,1), yaxt="n",xaxt="n")
lines(dim_list, H_jd_power, type = "o", pch = 19, col = "blue")
lines(dim_list, H_mat_power, type = "o", pch = 19, col = "purple")
lines(dim_list, H_hsic_power, type = "o", pch = 19, col = "orange")
legend("topleft", inset = c(0, 0.3),
       legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = c(seq(0, 1, 0.2)), cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at = dim_list, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Dimension", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()



