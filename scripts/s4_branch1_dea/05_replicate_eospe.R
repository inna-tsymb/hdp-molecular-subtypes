#!/usr/bin/env Rscript
# s4_branch1_dea/05_replicate_eospe.R
# -------------------------------------------------------------------
# INDEPENDENT REPLICATION of the EOPE signature on GSE148241 (EOSPE, early-onset
# severe). Question: does the hypoxic HIF-1 signature of early PE
# replicate in an independent early cohort?
#
# 1) own DESeq2 (~ category) on the full gene space of 148241;
# 2) DEG overlap with EOPE (GSE203507) + direction concordance;
# 3) enrichment (GO+KEGG) -> do hypoxia/HIF-1 terms replicate.
#
# Input:  data/gse148241_counts.csv, data/gse148241_meta.csv
# Output: results/deg_GSE148241_sig.csv, results/replication_eospe.txt,
#        results/enrich_GSE148241_*.csv
#
# Run:  Rscript scripts/s4_branch1_dea/05_replicate_eospe.R
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2); library(clusterProfiler); library(org.Hs.eg.db)
})
PADJ <- 0.05; LFC <- 1.0
dir.create("results", showWarnings = FALSE)

counts <- as.matrix(read.csv("data/gse148241_counts.csv", row.names = 1, check.names = FALSE))
meta   <- read.csv("data/gse148241_meta.csv", stringsAsFactors = FALSE)
stopifnot(all(meta$sample == colnames(counts)))
meta$category <- factor(meta$category, levels = c("Control", "Case"))

min_grp <- min(table(meta$category))
keep <- rowSums(counts >= 10) >= min_grp
counts <- counts[keep, ]
cat("Genes after filter:", nrow(counts), "| balance:",
    paste(names(table(meta$category)), table(meta$category), collapse = " / "), "\n")

dds <- DESeqDataSetFromMatrix(counts, meta, ~ category)
# sparse matrix: try poscounts first; if it fails, fall back to library-size normalization
dds <- tryCatch({
  d <- estimateSizeFactors(dds, type = "poscounts")
  if (any(is.na(sizeFactors(d)))) stop("NA size factors")
  d
}, error = function(e) {
  message("poscounts failed (", conditionMessage(e), ") -> library-size normalization")
  sf <- colSums(counts)
  sizeFactors(dds) <- sf / exp(mean(log(sf)))
  dds
})
dds <- DESeq(dds, quiet = TRUE)
res <- lfcShrink(dds, coef = "category_Case_vs_Control", type = "apeglm")
df <- as.data.frame(res); df$gene <- rownames(df)
sig <- subset(df, !is.na(padj) & padj < PADJ & abs(log2FoldChange) > LFC)
sig <- sig[order(sig$padj), ]
write.csv(sig, "results/deg_GSE148241_sig.csv", row.names = FALSE)
cat("DEG GSE148241 (EOSPE):", nrow(sig), "\n")

# ── replication vs EOPE (GSE203507) ───────────────────────────────────
eope <- read.csv("results/deg_GSE203507_sig.csv")
ov <- intersect(sig$gene, eope$gene)
cat("Overlap with EOPE (203507):", length(ov), "genes\n")

conc <- NA
if (length(ov) > 0) {
  rownames(eope) <- eope$gene
  s_dir <- sign(sig$log2FoldChange[match(ov, sig$gene)])
  e_dir <- sign(eope[ov, "log2FoldChange"])
  conc <- sum(s_dir == e_dir)
  cat("Concordant direction:", conc, "/", length(ov), "\n")
}

# ── enrichment: does hypoxia replicate? ────────────────────────────
universe <- rownames(counts)
go <- tryCatch(enrichGO(sig$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                        ont = "BP", universe = universe, pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10), error = function(e) NULL)
kegg <- tryCatch({
  ent <- bitr(sig$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  uni <- bitr(universe, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  enrichKEGG(ent, organism = "hsa", universe = uni, pvalueCutoff = 0.05)
}, error = function(e) NULL)

hyp_terms <- character(0)
if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
  gd <- as.data.frame(go); write.csv(gd, "results/enrich_GSE148241_GO.csv", row.names = FALSE)
  hyp_terms <- c(hyp_terms, grep("hypox|oxygen|HIF", gd$Description,
                                 ignore.case = TRUE, value = TRUE))
}
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  kd <- as.data.frame(kegg); write.csv(kd, "results/enrich_GSE148241_KEGG.csv", row.names = FALSE)
  hyp_terms <- c(hyp_terms, grep("HIF|hypox", kd$Description,
                                 ignore.case = TRUE, value = TRUE))
}

# ── replication summary ──────────────────────────────────────────────
sink("results/replication_eospe.txt")
cat("EOPE signature replication on GSE148241 (EOSPE)\n")
cat("=============================================\n")
cat("DEG 148241:", nrow(sig), " | DEG EOPE(203507):", nrow(eope), "\n")
cat("Overlap:", length(ov), " | concordant direction:", conc, "/", length(ov), "\n\n")
cat("Hypoxia/HIF terms in 148241:\n")
if (length(hyp_terms) > 0) {
  cat(paste0("  - ", unique(hyp_terms), "\n"))
} else {
  cat("  (no obvious hypoxia/HIF terms)\n")
}
cat("\nConclusion: if overlap is present with concordant direction and hypoxia terms\n")
cat("replicate -> the EOPE hypoxia signature is supported independently.\n")
sink()

cat("\nHypoxia terms:", if (length(hyp_terms)) paste(unique(hyp_terms), collapse="; ") else "none", "\n")
cat("Saved: results/deg_GSE148241_sig.csv, replication_eospe.txt, enrich_GSE148241_*.csv\n")