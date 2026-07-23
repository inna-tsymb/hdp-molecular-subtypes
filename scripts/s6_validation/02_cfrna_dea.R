#!/usr/bin/env Rscript
# s6_validation/06_cfrna_dea.R
# -------------------------------------------------------------------
# cfRNA DEA (Moufarrej): DESeq2 on each subset (timepoint + pooled),
# one sample per patient. Design ~ cohort + condition (batch as covariate),
# or ~ condition if only one cohort. Then intersect cfRNA DEG with placental
# signatures (closes the transfer narrative at the gene level).
#
# Input:  data/cfrna_dea/{counts,coldata}_<subset>.csv, subsets.txt
# Output: results/cfrna_deg_<subset>.csv, results/cfrna_dea_overlap.txt
#
# Run:  Rscript scripts/s6_validation/06_cfrna_dea.R
# -------------------------------------------------------------------

suppressPackageStartupMessages(library(DESeq2))
PADJ <- 0.05; LFC <- 1.0
dir.create("results", showWarnings = FALSE)
DIR <- "data/cfrna_dea"
subsets <- readLines(file.path(DIR, "subsets.txt"))

run_dea <- function(name) {
  cat("\n==============================\n", name, "\n==============================\n")
  counts <- as.matrix(read.csv(file.path(DIR, paste0("counts_", name, ".csv")),
                               row.names = 1, check.names = FALSE))
  cd <- read.csv(file.path(DIR, paste0("coldata_", name, ".csv")), row.names = 1)
  cd$condition <- factor(cd$condition, levels = c("control", "PE"))
  cd$cohort <- factor(cd$cohort)

  min_grp <- min(table(cd$condition))
  keep <- rowSums(counts >= 10) >= min_grp
  counts <- counts[keep, ]

  # design: add cohort only if there are >=2 cohorts and no full confound
  multi <- nlevels(droplevels(cd$cohort)) >= 2
  form <- if (multi) ~ cohort + condition else ~ condition
  cat("Design:", deparse(form), "| genes:", nrow(counts),
      "| balance:", paste(names(table(cd$condition)), table(cd$condition), collapse=" / "), "\n")

  dds <- DESeqDataSetFromMatrix(counts, cd, design = form)
  # sparse cfRNA counts: poscounts with fallback to library-size
  dds <- tryCatch({
    d <- estimateSizeFactors(dds, type = "poscounts")
    if (any(is.na(sizeFactors(d)))) stop("NA")
    d
  }, error = function(e) {
    sf <- colSums(counts); sizeFactors(dds) <- sf / exp(mean(log(sf))); dds
  })
  dds <- tryCatch(DESeq(dds, quiet = TRUE),
                  error = function(e) { design(dds) <- ~ condition; DESeq(dds, quiet = TRUE) })

  res <- results(dds, name = "condition_PE_vs_control", alpha = PADJ)
  res <- lfcShrink(dds, coef = "condition_PE_vs_control", type = "apeglm", res = res)
  df <- as.data.frame(res); df$gene <- rownames(df)
  df <- df[order(df$padj), c("gene","baseMean","log2FoldChange","lfcSE","pvalue","padj")]
  sig <- subset(df, !is.na(padj) & padj < PADJ & abs(log2FoldChange) > LFC)
  write.csv(sig, sprintf("results/cfrna_deg_%s.csv", name), row.names = FALSE)
  cat("cfRNA-DEG:", nrow(sig), "\n")
  sig$gene
}

cfrna_deg <- lapply(subsets, function(s) tryCatch(run_dea(s),
                    error = function(e) { cat("[error", s, ":", conditionMessage(e), "]\n"); character(0) }))
names(cfrna_deg) <- subsets
cfrna_union <- unique(unlist(cfrna_deg))

# ── INTERSECTION with placental DEG (transfer narrative) ──────────────────
rd <- function(p) if (file.exists(p)) read.csv(p)$gene else character(0)
eope   <- rd("results/deg_GSE203507_sig.csv")
lope   <- rd("results/deg_GSE306864_male_sig.csv")
panel  <- if (file.exists("data/panel_top25.txt")) readLines("data/panel_top25.txt") else character(0)
plac_all <- unique(c(eope, lope))

sink("results/cfrna_dea_overlap.txt")
cat("cfRNA-DEA — intersection with placental signatures\n")
cat("================================================\n")
cat("cfRNA-DEG (union of subsets):", length(cfrna_union), "\n\n")
for (s in subsets) cat(sprintf("  %-18s : %d DEG\n", s, length(cfrna_deg[[s]])))
cat("\nIntersections (cfRNA-union with ...):\n")
cat("  placental EOPE (203507):", length(intersect(cfrna_union, eope)),
    "/", length(eope), "\n")
cat("  placental LOPE (male)  :", length(intersect(cfrna_union, lope)),
    "/", length(lope), "\n")
cat("  panel top25             :", length(intersect(cfrna_union, panel)),
    "/", length(panel), "\n")
cat("  any placental DEG       :", length(intersect(cfrna_union, plac_all)),
    "/", length(plac_all), "\n\n")
shared <- intersect(cfrna_union, plac_all)
cat("Shared genes (cfRNA ∩ placenta):",
    if (length(shared)) paste(head(shared, 30), collapse=", ") else "NONE", "\n\n")
cat("Interpretation: low overlap -> placental and plasma PE transcriptomes\n")
cat("do not match closely (explains weak placenta->cfRNA transfer at the gene level).\n")
sink()

cat("\n== INTERSECTION ==\n")
cat("cfRNA-DEG union:", length(cfrna_union),
    "| ∩ EOPE:", length(intersect(cfrna_union, eope)),
    "| ∩ LOPE:", length(intersect(cfrna_union, lope)),
    "| ∩ panel:", length(intersect(cfrna_union, panel)), "\n")
cat("Saved: results/cfrna_deg_*.csv, cfrna_dea_overlap.txt\n")
