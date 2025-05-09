
TEXTWIDTH = 6.0689
TEXTHEIGHT = 9.33476

# specify the result directory
summary_dir <- "dependence_recovery/results/summary"

# create plotting directory
figure_dir <- "dependence_recovery/figure"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir)
  cat("Directory created:", figure_dir, "\n")
} else {
  cat("Directory already exists:", figure_dir, "\n")
}


## Gaussian setting
n <- seq(40, 200, by = 40)
correlation <- seq(0.1, 0.25, by = 0.05)
Gaussian_power <- readRDS(sprintf("%s/Gaussian_diff_set.rds", summary_dir))
Gaussian_consistency <- readRDS(sprintf("%s/Gaussian_consistency.rds", summary_dir))
Gaussian_fdr <- readRDS(sprintf("%s/Gaussian_fdr.rds", summary_dir))

## plot
pdf(file = sprintf("%s/Gaussian_diff_set.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.7*TEXTHEIGHT)
plot(n, Gaussian_power[1, ], type = "o", pch = 19, col = "red",
     main="Gaussian dependence set recovery", ylab="", xlab = "",ylim = c(0,5), yaxt="n",xaxt="n")
lines(n, Gaussian_power[2, ], type = "o", pch = 19, col = "blue")
lines(n, Gaussian_power[3, ], type = "o", pch = 19, col = "purple")
lines(n, Gaussian_power[4, ], type = "o", pch = 19, col = "orange")
legend("bottomleft", legend=c("a = 0.10", "a = 0.15", "a = 0.20", "a = 0.25"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 5, 1)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Sample size", side=1, line = 1.6, cex=0.9)
mtext("Set difference", side=2, line = 1.6, cex=0.9)

dev.off()

# Gaussian_consistency 

pdf(file = sprintf("%s/Gaussian_consisteny.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.7*TEXTHEIGHT)


plot(n, Gaussian_consistency[1, ], type = "o", pch = 19, col = "red",
     main="Gaussian dependence set recovery", ylab="", xlab = "",ylim = c(0,1), yaxt="n",xaxt="n")
lines(n, Gaussian_consistency[2, ], type = "o", pch = 19, col = "blue")
lines(n, Gaussian_consistency[3, ], type = "o", pch = 19, col = "purple")
lines(n, Gaussian_consistency[4, ], type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("a = 0.10", "a = 0.15", "a = 0.20", "a = 0.25"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0,1,0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Sample size", side=1, line = 1.6, cex=0.9)
mtext("Consistency rate", side=2, line = 1.6, cex=0.9)

dev.off()


# FDR

pdf(file = sprintf("%s/Gaussian_fdr.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.7*TEXTHEIGHT)


plot(n, Gaussian_fdr[1, ], type = "o", pch = 19, col = "red",
     main="Gaussian dependence set recovery", ylab="", xlab = "",ylim = c(0.04,0.16), yaxt="n",xaxt="n")
lines(n, Gaussian_fdr[2, ], type = "o", pch = 19, col = "blue")
lines(n, Gaussian_fdr[3, ], type = "o", pch = 19, col = "purple")
lines(n, Gaussian_fdr[4, ], type = "o", pch = 19, col = "orange")
legend("topleft", legend=c("a = 0.10", "a = 0.15", "a = 0.20", "a = 0.25"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0,.2,0.02)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Sample size", side=1, line = 1.6, cex=0.9)
mtext("False discovery rate", side=2, line = 1.6, cex=0.9)

dev.off()


# TEXTWIDTH = 6.0689
# TEXTHEIGHT = 9.33476
# 
# 
# ## Student setting
# n <- seq(40, 200, by = 40)
# df <- c(1, 2, 5, 10, Inf)
# Student_power <- readRDS(sprintf("%s/Student_diff_set.rds", summary_dir))
# Student_consistency <- readRDS(sprintf("%s/Student_consistency.rds", summary_dir))
# Student_fdr <- readRDS(sprintf("%s/Student_fdr.rds", summary_dir))
# 
# ## plot
# pdf(file = sprintf("%s/Student_diff_set.pdf", figure_dir),   # The directory you want to save the file in
#     width =  0.8*TEXTWIDTH, # The width of the plot in inches
#     height = 0.7*TEXTHEIGHT)
# plot(n, Student_power[1, ], type = "o", pch = 19, col = "red",
#      main="Higher order dependence set recovery", ylab="", xlab = "",ylim = c(0, 1.5), yaxt="n",xaxt="n")
# lines(n, Student_power[2, ], type = "o", pch = 19, col = "blue")
# lines(n, Student_power[3, ], type = "o", pch = 19, col = "purple")
# lines(n, Student_power[4, ], type = "o", pch = 19, col = "orange")
# lines(n, Student_power[5, ], type = "o", pch = 19, col = "green")
# legend("bottomright", legend=c("df = 1", "df = 2", "df = 5", "df = 10", "df = Inf"),
#        col=c("red", "blue", "purple", "orange", "green"), lty = 1, cex=0.7)
# axis(side = 2, at=c(seq(0, 1.5, 0.3)),cex.axis=0.7, mgp=c(1,0.6,0))
# axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
# mtext("Sample size", side=1, line = 1.6, cex=0.9)
# mtext("Size difference between discovery and oracle sets", side=2, line = 1.6, cex=0.9)
# 
# dev.off()
# 
# # Student_consistency 
# 
# pdf(file = sprintf("%s/Student_consisteny.pdf", figure_dir),   # The directory you want to save the file in
#     width =  0.8*TEXTWIDTH, # The width of the plot in inches
#     height = 0.7*TEXTHEIGHT)
# 
# 
# plot(n, Student_consistency[1, ], type = "o", pch = 19, col = "red",
#      main="Higher order dependence set recovery", ylab="", xlab = "",ylim = c(0.65,1), yaxt="n",xaxt="n")
# lines(n, Student_consistency[2, ], type = "o", pch = 19, col = "blue")
# lines(n, Student_consistency[3, ], type = "o", pch = 19, col = "purple")
# lines(n, Student_consistency[4, ], type = "o", pch = 19, col = "orange")
# lines(n, Student_consistency[5, ], type = "o", pch = 19, col = "green")
# legend("bottomright", legend=c("df = 1", "df = 2", "df = 5", "df = 10", "df = Inf"),
#        col=c("red", "blue", "purple", "orange", "green"), lty = 1, cex=0.7)
# axis(side = 2, at=c(seq(0,1,0.05)),cex.axis=0.7, mgp=c(1,0.6,0))
# axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
# mtext("Sample size", side=1, line = 1.6, cex=0.9)
# mtext("Consistency rate", side=2, line = 1.6, cex=0.9)
# 
# dev.off()
# 
# 
# # FDR
# 
# pdf(file = sprintf("%s/Student_fdr.pdf", figure_dir),   # The directory you want to save the file in
#     width =  0.8*TEXTWIDTH, # The width of the plot in inches
#     height = 0.7*TEXTHEIGHT)
# 
# 
# plot(n, Student_fdr[1, ], type = "o", pch = 19, col = "red",
#      main="Higher order dependence set recovery", ylab="", xlab = "",ylim = c(0.1,0.26), yaxt="n",xaxt="n")
# lines(n, Student_fdr[2, ], type = "o", pch = 19, col = "blue")
# lines(n, Student_fdr[3, ], type = "o", pch = 19, col = "purple")
# lines(n, Student_fdr[4, ], type = "o", pch = 19, col = "orange")
# lines(n, Student_fdr[5, ], type = "o", pch = 19, col = "green")
# legend("bottomright", legend=c("df = 1", "df = 2", "df = 5", "df = 10", "df = Inf"),
#        col=c("red", "blue", "purple", "orange", "green"), lty = 1, cex=0.7)
# axis(side = 2, at=c(seq(0.1,.3,0.04)),cex.axis=0.7, mgp=c(1,0.6,0))
# axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
# mtext("Sample size", side=1, line = 1.6, cex=0.9)
# mtext("False discovery rate", side=2, line = 1.6, cex=0.9)
# 
# dev.off()

## Student setting
n <- seq(40, 200, by = 40)
noise_level <- seq(0, 1.5, length.out = 4)
Student_power <- readRDS(sprintf("%s/Student_diff_set_noise.rds", summary_dir))
Student_consistency <- readRDS(sprintf("%s/Student_consistency_noise.rds", summary_dir))
Student_fdr <- readRDS(sprintf("%s/Student_fdr_noise.rds", summary_dir))

## plot
pdf(file = sprintf("%s/Student_diff_set_noise.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.7*TEXTHEIGHT)
plot(n, Student_power[1, ], type = "o", pch = 19, col = "red",
     main="Higher order dependence set recovery", ylab="", xlab = "",ylim = c(0, 4), yaxt="n",xaxt="n")
lines(n, Student_power[2, ], type = "o", pch = 19, col = "blue")
lines(n, Student_power[3, ], type = "o", pch = 19, col = "purple")
lines(n, Student_power[4, ], type = "o", pch = 19, col = "orange")
legend("topright", legend=c("sd = 0", "sd = 0.5", "sd = 1", "sd = 1.5"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 4, 1)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Sample size", side=1, line = 1.6, cex=0.9)
mtext("Set difference", side=2, line = 1.6, cex=0.9)

dev.off()

# Student_consistency 

pdf(file = sprintf("%s/Student_consisteny_noise.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.7*TEXTHEIGHT)


plot(n, Student_consistency[1, ], type = "o", pch = 19, col = "red",
     main="Higher order dependence set recovery", ylab="", xlab = "",ylim = c(0.05,1), yaxt="n",xaxt="n")
lines(n, Student_consistency[2, ], type = "o", pch = 19, col = "blue")
lines(n, Student_consistency[3, ], type = "o", pch = 19, col = "purple")
lines(n, Student_consistency[4, ], type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("sd = 0", "sd = 0.5", "sd = 1", "sd = 1.5"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0, 1, 0.2)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=n,cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Sample size", side=1, line = 1.6, cex=0.9)
mtext("Consistency rate", side=2, line = 1.6, cex=0.9)

dev.off()


# FDR

pdf(file = sprintf("%s/Student_fdr_noise.pdf", figure_dir),   # The directory you want to save the file in
    width =  0.8*TEXTWIDTH, # The width of the plot in inches
    height = 0.7*TEXTHEIGHT)


plot(n, Student_fdr[1, ], type = "o", pch = 19, col = "red",
     main="Higher order dependence set recovery", ylab="", xlab = "",ylim = c(0.04, 0.12), yaxt="n",xaxt="n")
lines(n, Student_fdr[2, ], type = "o", pch = 19, col = "blue")
lines(n, Student_fdr[3, ], type = "o", pch = 19, col = "purple")
lines(n, Student_fdr[4, ], type = "o", pch = 19, col = "orange")
legend("bottomright", legend=c("sd = 0", "sd = 0.5", "sd = 1", "sd = 1.5"),
       col=c("red", "blue", "purple", "orange"), lty = 1, cex=0.7)
axis(side = 2, at=c(seq(0.02, 0.12, 0.02)),cex.axis=0.7, mgp=c(1,0.6,0))
axis(side = 1, at=n, cex.axis=0.7, mgp=c(1,0.6,0))
mtext("Sample size", side=1, line = 1.6, cex=0.9)
mtext("False discovery rate", side=2, line = 1.6, cex=0.9)

dev.off()
