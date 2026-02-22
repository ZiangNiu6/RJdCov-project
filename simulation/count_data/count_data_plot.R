
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476
alpha_list <- 0.05

# Power computation
f1 <- function(p_value, alpha_list){
  sapply(alpha_list, function(x) length(which(p_value < x)) / length(p_value))
}

# create figure folder
figure_dir <- "simulation/new_figures/count_data"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}

# ============================================================
# Cauchy discretized
# ============================================================
cauchy_res <- readRDS("simulation/results/cauchy_discretized_n200.rds")
signal_list <- c(0, 0.2, 0.4, 0.6, 0.8)

cauchy_rd_power  <- apply(cauchy_res[, "RJdCov", ], 2, function(x) f1(x, alpha_list))
cauchy_jd_power  <- apply(cauchy_res[, "JdCov", ],  2, function(x) f1(x, alpha_list))
cauchy_mat_power <- apply(cauchy_res[, "Matteson",], 2, function(x) f1(x, alpha_list))
cauchy_hsic_power <- apply(cauchy_res[, "dHSIC", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Cauchy_discretized.pdf", figure_dir),
    width =  0.8*TEXTWIDTH,
    height = 0.45*TEXTHEIGHT)

plot(signal_list, cauchy_rd_power, type = "o", pch = 19, col = "red",
     main="Cauchy discretized", ylab="", xlab = "", ylim = c(0, 1), yaxt="n", xaxt="n")
abline(h = 0.05, col = "red", lty = 2)
lines(signal_list, cauchy_rd_power, type = "o", pch = 19, col = "red")
lines(signal_list, cauchy_jd_power, type = "o", pch = 19, col = "blue")
lines(signal_list, cauchy_mat_power, type = "o", pch = 19, col = "purple")
lines(signal_list, cauchy_hsic_power, type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = seq(0, 1, 0.2), cex.axis = 0.7, mgp = c(1, 0.6, 0))
axis(side = 1, at = signal_list, cex.axis = 0.7, mgp = c(1, 0.6, 0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()

# ============================================================
# Poisson discretized
# ============================================================
poisson_res <- readRDS("simulation/results/poisson_discretized_n200.rds")
a_list <- c(0, 0.05, 0.1, 0.15, 0.2)

poisson_rd_power  <- apply(poisson_res[, "RJdCov", ], 2, function(x) f1(x, alpha_list))
poisson_jd_power  <- apply(poisson_res[, "JdCov", ],  2, function(x) f1(x, alpha_list))
poisson_mat_power <- apply(poisson_res[, "Matteson",], 2, function(x) f1(x, alpha_list))
poisson_hsic_power <- apply(poisson_res[, "dHSIC", ], 2, function(x) f1(x, alpha_list))

pdf(file = sprintf("%s/Poisson_discretized.pdf", figure_dir),
    width =  0.8*TEXTWIDTH,
    height = 0.45*TEXTHEIGHT)

plot(a_list, poisson_rd_power, type = "o", pch = 19, col = "red",
     main="Poisson discretized", ylab="", xlab = "", ylim = c(0, 1), yaxt="n", xaxt="n")
abline(h = 0.05, col = "red", lty = 2)
lines(a_list, poisson_rd_power, type = "o", pch = 19, col = "red")
lines(a_list, poisson_jd_power, type = "o", pch = 19, col = "blue")
lines(a_list, poisson_mat_power, type = "o", pch = 19, col = "purple")
lines(a_list, poisson_hsic_power, type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("Proposed", "JdCov", "MT", "dHSIC"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at = seq(0, 1, 0.2), cex.axis = 0.7, mgp = c(1, 0.6, 0))
axis(side = 1, at = a_list, cex.axis = 0.7, mgp = c(1, 0.6, 0))
mtext("Magnitude of signal", side=1, line = 1.6, cex=0.9)
mtext("Power", side=2, line = 1.6, cex=0.9)

dev.off()
