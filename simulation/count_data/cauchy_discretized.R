source("simulation/Test_function.R")

log_file <- "simulation/results/log_cauchy_discretized.txt"
con <- file(log_file, open="wt")
sink(con); sink(con, type="message")

n <- 200; B <- 500; alpha <- 0.05
group <- c(1,1,1,2,2,2,3,3,3)
signal_list <- c(0, 0.2, 0.4, 0.6, 0.8)

set.seed(1)
cat("Precomputing null distribution (n=200, dim_list=c(3,3,3))...\n")
null_jdcov <- gensamdistrjdcov(n, dim_list=rep(3,3), niter=1000)
cat("Done.\n\n")

joint_methods <- c("RJdCov","JdCov","dHSIC","Matteson")
results <- array(NA, dim=c(B, length(joint_methods), length(signal_list)),
                 dimnames=list(NULL, joint_methods, paste0("s=",signal_list)))

for (si in seq_along(signal_list)) {
  signal <- signal_list[si]
  cat(sprintf("=== signal=%.1f ===\n", signal))
  t0 <- Sys.time()

  for (b in 1:B) {
    X_cont <- matrix(rcauchy(n*9), nrow=n, ncol=9) + signal * rcauchy(n)
    X <- round(abs(X_cont))
    X_list <- list(X[,1:3], X[,4:6], X[,7:9])

    stat <- computestatisticjdcov(X, dim_list=rep(3,3))
    results[b,"RJdCov",si] <- length(which(null_jdcov >= stat))/(length(null_jdcov)+1)

    results[b,"JdCov",si] <- jdcov.test(X_list, stat.type="V", B=500, alpha=0.05)$p.value

    results[b,"dHSIC",si] <- dhsic.test(X_list, B=500)$p.value

    results[b,"Matteson",si] <- boot_matteson(X, B=500, group=group)

    if (b %% 50 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
      cat(sprintf("  rep %d/%d (%.1f min)\n", b, B, elapsed))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  cat(sprintf("  Done signal=%.1f in %.1f min\n\n", signal, elapsed))
}

cat("\n=== Rejection rates (alpha=0.05) ===\n")
cat(sprintf("%-10s", "Method"))
for (s in signal_list) cat(sprintf("  s=%.1f", s))
cat("\n")
for (m in joint_methods) {
  cat(sprintf("%-10s", m))
  for (si in seq_along(signal_list)) {
    pv <- results[,m,si]; pv <- pv[!is.na(pv)]
    cat(sprintf("  %.3f", mean(pv < alpha)))
  }
  cat("\n")
}

saveRDS(results, "simulation/results/cauchy_discretized_n200.rds")
cat("\nSaved to simulation/results/cauchy_discretized_n200.rds\n")

sink(type="message"); sink()
close(con)
