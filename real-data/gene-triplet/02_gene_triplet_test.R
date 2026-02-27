###############################################################################
# 02_gene_triplet_test.R
# Test all 216 single-gene triplets (one from each module: CYTO, ECM, PROLIF)
# for pairwise independence and 3-way dependence.
#
# Per triplet:
#   Pairwise: RdCov (null precomputed) + dCov (bootstrap, B=500)
#   3-way joint: RJdCov (null precomputed) + JdCov (B=500)
#   3-way total: RdCov (null precomputed) + HodCov (bootstrap, B=500)
#
# Flags triplets with "higher-order dependence only":
#   all pairwise Holm > 0.05, but joint test <= 0.05
###############################################################################

source("real-data/function.R")
library(jdcov)

# --- Load individual gene Z-scores ---
dat <- readRDS("real-data/gene-triplet/TCGA_SARC_gene_zscores.rds")
gene_zscores <- dat$gene_zscores
modules <- dat$modules
n <- nrow(gene_zscores)

cat(sprintf("Loaded gene Z-scores: %d samples x %d genes\n",
            nrow(gene_zscores), ncol(gene_zscores)))
cat(sprintf("CYTO genes:   %s\n", paste(modules$CYTO, collapse = ", ")))
cat(sprintf("ECM genes:    %s\n", paste(modules$ECM, collapse = ", ")))
cat(sprintf("PROLIF genes: %s\n", paste(modules$PROLIF, collapse = ", ")))

B <- 500
B_null <- 2000

# --- Helpers ---
p_from_null <- function(stat, nullvec) {
  (sum(nullvec >= stat) + 1) / (length(nullvec) + 1)
}

boot_sample_list <- function(Xlist, n, d) {
  resample_X <- list()
  for (i in 1:d) {
    boot_idx <- sample(1:n, n, replace = TRUE)
    resample_X[[i]] <- as.matrix(Xlist[[i]][boot_idx, , drop = FALSE])
  }
  resample_X
}

boot_hodcov <- function(Xlist, B) {
  test_stat <- hodcov(Xlist)
  boot <- numeric(B)
  n <- nrow(Xlist[[1]])
  for (b in 1:B) {
    boot[b] <- hodcov(boot_sample_list(Xlist, n, d = length(Xlist)))
  }
  (sum(boot >= test_stat) + 1) / (B + 1)
}

# --- Precompute null distributions (depend only on n, reused for all 216) ---
set.seed(42)
cat("\nPrecomputing null distributions...\n")

cat("  Pairwise RdCov null (dim_list=c(1,1))...\n")
null2_rdcov <- gensamdistrhodcov(N = n, dim_list = c(1, 1), niter = B_null)

cat("  3-way RdCov null (dim_list=c(1,1,1))...\n")
null3_rdcov <- gensamdistrhodcov(N = n, dim_list = c(1, 1, 1), niter = B_null)

cat("  3-way RJdCov null (dim_list=c(1,1,1))...\n")
null3_rjdcov <- gensamdistrjdcov(N = n, dim_list = c(1, 1, 1), niter = B_null)

cat("Null distributions ready.\n\n")

# --- Test one triplet ---
test_one_triplet <- function(x1, x2, x3, B = 500) {
  X_pair <- cbind(x1, x2, x3)

  # Pairwise RdCov
  stat12_rd <- computestatisticrdcov(cbind(x1, x2), dim_list = c(1, 1), n = n)
  stat13_rd <- computestatisticrdcov(cbind(x1, x3), dim_list = c(1, 1), n = n)
  stat23_rd <- computestatisticrdcov(cbind(x2, x3), dim_list = c(1, 1), n = n)
  p12_rd <- p_from_null(stat12_rd, null2_rdcov)
  p13_rd <- p_from_null(stat13_rd, null2_rdcov)
  p23_rd <- p_from_null(stat23_rd, null2_rdcov)

  # Pairwise dCov (bootstrap)
  p12_dcov <- boot_hodcov(list(as.matrix(x1), as.matrix(x2)), B = B)
  p13_dcov <- boot_hodcov(list(as.matrix(x1), as.matrix(x3)), B = B)
  p23_dcov <- boot_hodcov(list(as.matrix(x2), as.matrix(x3)), B = B)

  # Holm correction (3 pairwise tests per method)
  rd_holm   <- p.adjust(c(p12_rd, p13_rd, p23_rd), method = "holm")
  dcov_holm <- p.adjust(c(p12_dcov, p13_dcov, p23_dcov), method = "holm")

  # 3-way RJdCov
  stat_rjdcov <- computestatisticjdcov(X_pair, dim_list = c(1, 1, 1), n = n)
  p_rjdcov <- p_from_null(stat_rjdcov, null3_rjdcov)

  # 3-way JdCov
  X_list <- list(as.matrix(x1), as.matrix(x2), as.matrix(x3))
  jdcov_res <- jdcov.test(X_list, stat.type = "V", B = B, alpha = 0.05)
  p_jdcov <- jdcov_res$p.value

  # 3-way RdCov (higher-order)
  stat_rdcov <- computestatisticrdcov(X_pair, dim_list = c(1, 1, 1), n = n)
  p_rdcov <- p_from_null(stat_rdcov, null3_rdcov)

  # 3-way HodCov (bootstrap)
  p_hodcov <- boot_hodcov(X_list, B = B)

  data.frame(
    p12_rd = p12_rd, p13_rd = p13_rd, p23_rd = p23_rd,
    p12_rd_holm = rd_holm[1], p13_rd_holm = rd_holm[2], p23_rd_holm = rd_holm[3],
    p12_dcov = p12_dcov, p13_dcov = p13_dcov, p23_dcov = p23_dcov,
    p12_dcov_holm = dcov_holm[1], p13_dcov_holm = dcov_holm[2], p23_dcov_holm = dcov_holm[3],
    p_rjdcov = p_rjdcov, p_jdcov = p_jdcov, p_rdcov = p_rdcov, p_hodcov = p_hodcov,
    stringsAsFactors = FALSE
  )
}

# --- Triple loop over CYTO x ECM x PROLIF ---
cyto_genes   <- intersect(modules$CYTO, colnames(gene_zscores))
ecm_genes    <- intersect(modules$ECM, colnames(gene_zscores))
prolif_genes <- intersect(modules$PROLIF, colnames(gene_zscores))

total_combos <- length(cyto_genes) * length(ecm_genes) * length(prolif_genes)
cat(sprintf("Testing %d triplet combinations (%d x %d x %d)\n\n",
            total_combos, length(cyto_genes), length(ecm_genes), length(prolif_genes)))

results_list <- vector("list", total_combos)
combo_idx <- 0

for (g1 in cyto_genes) {
  for (g2 in ecm_genes) {
    for (g3 in prolif_genes) {
      combo_idx <- combo_idx + 1

      if (combo_idx %% 10 == 1) {
        cat(sprintf("[%d/%d] %s - %s - %s\n", combo_idx, total_combos, g1, g2, g3))
      }

      x1 <- as.matrix(gene_zscores[, g1])
      x2 <- as.matrix(gene_zscores[, g2])
      x3 <- as.matrix(gene_zscores[, g3])

      row <- test_one_triplet(x1, x2, x3, B = B)
      row <- cbind(data.frame(gene_CYTO = g1, gene_ECM = g2, gene_PROLIF = g3,
                              stringsAsFactors = FALSE), row)
      results_list[[combo_idx]] <- row
    }
  }
}

results <- do.call(rbind, results_list)

# --- Flag higher-order only triplets ---
results$higher_order_only_rd <- with(results,
  p12_rd_holm > 0.05 & p13_rd_holm > 0.05 & p23_rd_holm > 0.05 & p_rjdcov <= 0.05
)
results$higher_order_only_dcov <- with(results,
  p12_dcov_holm > 0.05 & p13_dcov_holm > 0.05 & p23_dcov_holm > 0.05 & p_jdcov <= 0.05
)

###############################################################################
# === SUMMARY ===
###############################################################################
cat("\n===============================================\n")
cat(sprintf("  Combinatorial test complete: %d triplets\n", nrow(results)))
cat("===============================================\n")

n_ho_rd   <- sum(results$higher_order_only_rd)
n_ho_dcov <- sum(results$higher_order_only_dcov)

cat(sprintf("\nHigher-order only (RdCov-based): %d / %d triplets\n", n_ho_rd, nrow(results)))
cat(sprintf("Higher-order only (dCov-based):  %d / %d triplets\n", n_ho_dcov, nrow(results)))

if (n_ho_rd > 0) {
  cat("\n--- Higher-order only triplets (RdCov-based) ---\n")
  ho_rd <- results[results$higher_order_only_rd, ]
  print(ho_rd[, c("gene_CYTO", "gene_ECM", "gene_PROLIF",
                   "p12_rd_holm", "p13_rd_holm", "p23_rd_holm", "p_rjdcov")],
        row.names = FALSE, digits = 4)
}

if (n_ho_dcov > 0) {
  cat("\n--- Higher-order only triplets (dCov-based) ---\n")
  ho_dcov <- results[results$higher_order_only_dcov, ]
  print(ho_dcov[, c("gene_CYTO", "gene_ECM", "gene_PROLIF",
                     "p12_dcov_holm", "p13_dcov_holm", "p23_dcov_holm", "p_jdcov")],
        row.names = FALSE, digits = 4)
}

# --- Save ---
outfile <- "real-data/gene-triplet/gene_combinatorial_results.rds"
saveRDS(results, file = outfile)
cat(sprintf("\nFull results saved to %s\n", outfile))
