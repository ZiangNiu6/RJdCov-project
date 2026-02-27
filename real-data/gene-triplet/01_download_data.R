###############################################################################
# 01_download_data.R
# Download TCGA-SARC RNA-seq counts, normalize, extract gene Z-scores, save RDS
###############################################################################

# --- Install / load packages ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (pkg in c("TCGAbiolinks", "SummarizedExperiment", "edgeR")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

# --- 1. Query & download TCGA-SARC STAR-Counts (primary tumors, open) ---
query <- GDCquery(
  project          = "TCGA-SARC",
  data.category    = "Transcriptome Profiling",
  data.type        = "Gene Expression Quantification",
  workflow.type    = "STAR - Counts",
  sample.type      = "Primary Tumor",
  access           = "open"
)

GDCdownload(query)
se <- GDCprepare(query)

# --- 2. Extract count matrix & map gene IDs ---
counts <- assay(se, "unstranded")
rd <- rowData(se)

# Use gene symbols from rowData; fall back to Ensembl IDs if unavailable
if ("gene_name" %in% colnames(rd)) {
  gene_symbols <- rd$gene_name
} else {
  gene_symbols <- rownames(counts)
}

# Remove duplicated / NA gene symbols (keep first occurrence)
keep <- !is.na(gene_symbols) & !duplicated(gene_symbols)
counts <- counts[keep, ]
rownames(counts) <- gene_symbols[keep]

# --- 3. Filter low-expression genes ---
counts <- counts[rowSums(counts >= 10) >= 10, ]

# --- 4. TMM normalization + logCPM ---
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE)

# --- 5. Gene-wise Z-scores ---
zscores <- t(scale(t(logcpm)))

# --- 6. Define 3 gene modules (6 genes each) ---
modules <- list(
  CYTO   = c("TRAC", "NKG7", "KLRD1", "PRF1", "GZMB", "GNLY"),
  ECM    = c("TGFB1", "SERPINE1", "COL1A1", "FN1", "ACTA2", "TAGLN"),
  PROLIF = c("MKI67", "TOP2A", "CDK1", "CCNB1", "MCM2", "UBE2C")
)

# --- 7. Save individual gene Z-scores for combinatorial analysis ---
all_genes <- unlist(modules)
gene_zscores <- t(zscores[intersect(all_genes, rownames(zscores)), ])
cat(sprintf("Gene Z-score matrix: %d samples x %d genes\n",
            nrow(gene_zscores), ncol(gene_zscores)))
saveRDS(list(gene_zscores = gene_zscores, modules = modules),
        file = "real-data/gene-triplet/TCGA_SARC_gene_zscores.rds")
cat("Saved gene Z-scores to real-data/gene-triplet/TCGA_SARC_gene_zscores.rds\n")
