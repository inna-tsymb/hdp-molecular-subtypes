#!/usr/bin/env Rscript
# 08_dea_stratified.R
# -------------------------------------------------------------------
# Branch 1, Option B: STRATIFIED DEA.
# Each dataset is analyzed SEPARATELY (its own batch -> ComBat_seq is NOT needed).
#   GSE203507  EOPE  vs control   design ~ category
#   GSE306864  LOPE  vs control   design ~ sex + category   (sex is an effect modifier!)
#   GSE190971  mixed vs control   design ~ category
# Then: compare DEG lists (shared core vs EOPE/LOPE-specific).
#
# Input:  data/dea_counts.csv (18325 x 91), data/dea_metadata.csv
# Output: results/deg_<GSE>.csv (full),
#         results/deg_<GSE>_sig.csv (padj<0.05),
#         results/deg_comparison.csv, results/deg_summary.txt
#
# Run:  Rscript scripts/08_dea_stratified.R
# Requires: DESeq2, apeglm  (BiocManager::install(c("DESeq2","apeglm")))
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
})

PADJ <- 0.05
LFC  <- 1.0                      # |log2FC| > 1 -> at least 2x change
dir.create("results", showWarnings = FALSE)

counts <- as.matrix(read.csv("data/dea_counts.csv", row.names = 1, check.names = FALSE))
meta   <- read.csv("data/dea_metadata.csv", stringsAsFactors = FALSE)
stopifnot(all(meta$sample == colnames(counts)))

# dataset-specific config: design formula
configs <- list(
  GSE203507 = list(design = ~ category,        onset = "EOPE"),
  GSE306864 = list(design = ~ sex + category,  onset = "LOPE"),
  GSE190971 = list(design = ~ category,        onset = "mixed")
)

run_one <- function(gse, cfg) {
  cat("\n==============================\n", gse, "(", cfg$onset, ")\n==============================\n")
  idx <- meta$dataset == gse
  cd  <- counts[, idx, drop = FALSE]
  md  <- meta[idx, , drop = FALSE]

  # Control must be the reference so log2FC is "Case vs Control"
  md$category <- factor(md$category, levels = c("Control", "Case"))
  if ("sex" %in% all.vars(cfg$design)) md$sex <- factor(md$sex)

  # low-count filter: gene must have >=10 total counts and be present
  # in at least the size of the smaller group
  min_grp <- min(table(md$category))
  keep <- rowSums(cd >= 10) >= min_grp
  cd <- cd[keep, ]
  cat("Genes after filter:", nrow(cd), "| samples:", ncol(cd),
      "| balance:", paste(names(table(md$category)), table(md$category), collapse=" / "), "\n")

  dds <- DESeqDataSetFromMatrix(cd, colData = md, design = cfg$design)
  dds <- DESeq(dds, quiet = TRUE)

  res <- results(dds, name = "category_Case_vs_Control", alpha = PADJ)
  # shrink for correct ranking by effect size
  res <- lfcShrink(dds, coef = "category_Case_vs_Control", type = "apeglm", res = res)

  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[order(df$padj), c("gene","baseMean","log2FoldChange","lfcSE","pvalue","padj")]

  sig <- subset(df, !is.na(padj) & padj < PADJ & abs(log2FoldChange) > LFC)
  cat("DEG (padj<", PADJ, " & |log2FC|>", LFC, "): ", nrow(sig), "\n", sep="")

  write.csv(df,  sprintf("results/deg_%s.csv", gse), row.names = FALSE)
  write.csv(sig, sprintf("results/deg_%s_sig.csv", gse), row.names = FALSE)
  list(all = df, sig = sig, onset = cfg$onset)
}

out <- lapply(names(configs), function(g) run_one(g, configs[[g]]))
names(out) <- names(configs)

# ── Compare DEG lists ───────────────────────────────────────────
eope <- out$GSE203507$sig$gene
lope <- out$GSE306864$sig$gene
mix  <- out$GSE190971$sig$gene

core       <- Reduce(intersect, list(eope, lope, mix))
eope_lope  <- intersect(eope, lope)
eope_only  <- setdiff(eope, union(lope, mix))
lope_only  <- setdiff(lope, union(eope, mix))

comp <- data.frame(
  set = c("EOPE (203507)","LOPE (306864)","mixed (190971)",
          "EOPE∩LOPE","core (all 3)","EOPE-only","LOPE-only"),
  n   = c(length(eope), length(lope), length(mix),
          length(eope_lope), length(core), length(eope_only), length(lope_only))
)
write.csv(comp, "results/deg_comparison.csv", row.names = FALSE)

# shared effect direction for EOPE∩LOPE (subtype consistency)
if (length(eope_lope) > 0) {
  e <- out$GSE203507$all; l <- out$GSE306864$all
  rownames(e) <- e$gene; rownames(l) <- l$gene
  same_dir <- sign(e[eope_lope,"log2FoldChange"]) == sign(l[eope_lope,"log2FoldChange"])
  cat("\nEOPE∩LOPE:", length(eope_lope), "genes, concordant direction:",
      sum(same_dir, na.rm=TRUE), "\n")
  write.csv(data.frame(gene=eope_lope, same_direction=same_dir),
            "results/deg_shared_direction.csv", row.names = FALSE)
}

sink("results/deg_summary.txt")
cat("STRATIFIED DEA — summary (Option B)\n")
cat("padj <", PADJ, " & |log2FC| >", LFC, "\n\n")
print(comp)
cat("\nInterpretation:\n")
cat("  core        -> shared PE signature (onset-independent)\n")
cat("  EOPE-only   -> early-specific genes\n")
cat("  LOPE-only   -> late-specific genes\n")
cat("  EOPE∩LOPE   -> subtype overlap (check direction concordance!)\n")
sink()

cat("\n== SUMMARY ==\n"); print(comp)
cat("\nSaved in results/:\n  deg_<GSE>.csv / _sig.csv, deg_comparison.csv,\n  deg_shared_direction.csv, deg_summary.txt\n")

suppressPackageStartupMessages(library(DESeq2))
counts <- as.matrix(read.csv("data/dea_counts.csv", row.names=1, check.names=FALSE))
meta   <- read.csv("data/dea_metadata.csv", stringsAsFactors=FALSE)
idx <- meta$dataset=="GSE306864" & meta$sex=="male"
cd  <- counts[, idx]; md <- meta[idx,]
md$category <- factor(md$category, levels=c("Control","Case"))
cat("male-only balance:", paste(names(table(md$category)), table(md$category)), "\n")
keep <- rowSums(cd>=10) >= min(table(md$category)); cd <- cd[keep,]
dds <- DESeq(DESeqDataSetFromMatrix(cd, md, ~category), quiet=TRUE)
res <- lfcShrink(dds, coef="category_Case_vs_Control", type="apeglm")
sig <- subset(as.data.frame(res), !is.na(padj) & padj<0.05 & abs(log2FoldChange)>1)
cat("male-only LOPE DEG:", nrow(sig), "\n")
