source("simulation/Test_function.R")

log_file <- "simulation/results/log_poisson_discretized.txt"
con <- file(log_file, open="wt")
sink(con); sink(con, type="message")

n <- 200; B <- 500; alpha <- 0.05
group <- c(1,1,2,2,3,3)
mu <- 20
a_list <- c(0, 0.05, 0.1, 0.15, 0.2)

set.seed(1)
cat("Precomputing null distribution (n=200, dim_list=c(2,2,2))...\n")
null_jdcov <- gensamdistrjdcov(n, dim_list=c(2,2,2), niter=1000)
cat("Done.\n\n")

joint_methods <- c("RJdCov","JdCov","dHSIC","Matteson")
results <- array(NA, dim=c(B, length(joint_methods), length(a_list)),
                 dimnames=list(NULL, joint_methods, paste0("a=",a_list)))

for (ai in seq_along(a_list)) {
  a <- a_list[ai]
  cat(sprintf("=== a=%.2f ===\n", a))
  t0 <- Sys.time()

  for (b in 1:B) {
    X1 <- matrix(rpois(n*2, lambda=mu), nrow=n, ncol=2)
    X2 <- matrix(rpois(n*2, lambda=mu), nrow=n, ncol=2)
    X3_prime <- matrix(rpois(n*2, lambda=mu), nrow=n, ncol=2)
    X3 <- round(a*(X1 + X2) + X3_prime)

    X <- cbind(X1, X2, X3)
    X_list <- list(X1, X2, X3)

    stat <- computestatisticjdcov(X, dim_list=c(2,2,2))
    results[b,"RJdCov",ai] <- length(which(null_jdcov >= stat))/(length(null_jdcov)+1)

    results[b,"JdCov",ai] <- jdcov.test(X_list, stat.type="V", B=500, alpha=0.05)$p.value

    results[b,"dHSIC",ai] <- dhsic.test(X_list, B=500)$p.value

    results[b,"Matteson",ai] <- boot_matteson(X, B=500, group=group)

    if (b %% 50 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
      cat(sprintf("  rep %d/%d (%.1f min)\n", b, B, elapsed))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  cat(sprintf("  Done a=%.2f in %.1f min\n\n", a, elapsed))
}

cat("\n=== Rejection rates (alpha=0.05) ===\n")
cat(sprintf("%-10s", "Method"))
for (a in a_list) cat(sprintf("  a=%.2f", a))
cat("\n")
for (m in joint_methods) {
  cat(sprintf("%-10s", m))
  for (ai in seq_along(a_list)) {
    pv <- results[,m,ai]; pv <- pv[!is.na(pv)]
    cat(sprintf("  %.3f", mean(pv < alpha)))
  }
  cat("\n")
}

saveRDS(results, "simulation/results/poisson_discretized_n200.rds")
cat("\nSaved to simulation/results/poisson_discretized_n200.rds\n")

sink(type="message"); sink()
close(con)
