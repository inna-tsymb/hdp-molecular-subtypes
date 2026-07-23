#!/usr/bin/env Rscript
# 09_enrichment.R
# -------------------------------------------------------------------
# Functional enrichment (KEGG + GO BP) for PE subtypes.
# Compares EOPE vs LOPE(male) at the PATHWAY level -> mechanistic validation
# of the hypothesis "early and late PE are different diseases".
#
# Gene sets:
#   EOPE       = DEG GSE203507                (283)
#   LOPE_male  = DEG GSE306864 male-only      (~117; computed and saved here)
#   shared     = intersection EOPE ∩ LOPE_male
# Universe    = 18325 tested genes (master_gene_list) — CRITICAL.
#
# Output: results/enrich_<set>_<GO|KEGG>.csv + dotplot PDFs,
#        results/enrich_pathway_compare.csv (shared vs specific pathways)
#
# Run:  Rscript scripts/09_enrichment.R
# Requires: clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2, DESeq2
#   BiocManager::install(c("clusterProfiler","org.Hs.eg.db","enrichplot"))
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2); library(DESeq2)
})

PADJ <- 0.05; LFC <- 1.0
dir.create("results", showWarnings = FALSE)

universe_sym <- readLines("data/master_gene_list.txt")            # 18325 genes
universe_ent <- bitr(universe_sym, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

# ── LOPE male-only: compute and save as a standalone result ─────
counts <- as.matrix(read.csv("data/dea_counts.csv", row.names=1, check.names=FALSE))
meta   <- read.csv("data/dea_metadata.csv", stringsAsFactors=FALSE)
idx <- meta$dataset=="GSE306864" & meta$sex=="male"
cd  <- counts[, idx]; md <- meta[idx, ]
md$category <- factor(md$category, levels=c("Control","Case"))
keep <- rowSums(cd>=10) >= min(table(md$category)); cd <- cd[keep, ]
dds <- DESeq(DESeqDataSetFromMatrix(cd, md, ~category), quiet=TRUE)
res <- lfcShrink(dds, coef="category_Case_vs_Control", type="apeglm")
res <- as.data.frame(res); res$gene <- rownames(res)
lope_male_sig <- subset(res, !is.na(padj) & padj<PADJ & abs(log2FoldChange)>LFC)
write.csv(lope_male_sig[order(lope_male_sig$padj),],
          "results/deg_GSE306864_male_sig.csv", row.names=FALSE)

# ── gene sets ─────────────────────────────────────────────────────
eope <- read.csv("results/deg_GSE203507_sig.csv")$gene
lope <- lope_male_sig$gene
gene_sets <- list(
  EOPE      = eope,
  LOPE_male = lope,
  shared    = intersect(eope, lope)
)
cat("Sets:", sapply(gene_sets, length), "\n")

# ── enrichment functions with universe ────────────────────────────────────
run_go <- function(genes) {
  enrichGO(genes, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="BP",
           universe=universe_sym, pAdjustMethod="BH",
           pvalueCutoff=0.05, qvalueCutoff=0.10, readable=FALSE)
}
run_kegg <- function(genes) {
  ent <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  tryCatch(
    enrichKEGG(ent, organism="hsa", universe=universe_ent,
               pAdjustMethod="BH", pvalueCutoff=0.05),
    error=function(e){ message("KEGG skipped (internet required): ", e$message); NULL }
  )
}

save_enrich <- function(obj, tag) {
  if (is.null(obj) || nrow(as.data.frame(obj))==0) {
    cat("  ", tag, ": empty\n"); return(invisible(NULL))
  }
  df <- as.data.frame(obj)
  write.csv(df, sprintf("results/enrich_%s.csv", tag), row.names=FALSE)
  n <- min(15, nrow(df))
  ggsave(sprintf("results/enrich_%s_dotplot.pdf", tag),
         dotplot(obj, showCategory=n) + ggtitle(tag),
         width=8, height=max(4, n*0.35))
  cat("  ", tag, ":", nrow(df), "terms\n")
}

for (nm in names(gene_sets)) {
  g <- gene_sets[[nm]]
  cat("\n==", nm, "(", length(g), "genes) ==\n")
  if (length(g) < 5) { cat("  too few genes for enrichment\n"); next }
  save_enrich(run_go(g),   paste0(nm, "_GO"))
  save_enrich(run_kegg(g), paste0(nm, "_KEGG"))
}

# ── pathway comparison: EOPE vs LOPE (shared vs specific) ───────────
top_terms <- function(tag) {
  f <- sprintf("results/enrich_%s.csv", tag)
  if (file.exists(f)) as.data.frame(read.csv(f))$Description else character(0)
}
for (db in c("GO","KEGG")) {
  e <- top_terms(paste0("EOPE_", db)); l <- top_terms(paste0("LOPE_male_", db))
  if (length(e)==0 && length(l)==0) next
  comp <- data.frame(
    term = union(e, l),
    in_EOPE = union(e,l) %in% e,
    in_LOPE_male = union(e,l) %in% l
  )
  comp$class <- ifelse(comp$in_EOPE & comp$in_LOPE_male, "shared",
                ifelse(comp$in_EOPE, "EOPE-specific", "LOPE-specific"))
  write.csv(comp, sprintf("results/enrich_compare_%s.csv", db), row.names=FALSE)
  cat("\n[", db, "] pathways: EOPE-spec=", sum(comp$class=="EOPE-specific"),
      " LOPE-spec=", sum(comp$class=="LOPE-specific"),
      " shared=", sum(comp$class=="shared"), "\n", sep="")
}

cat("\nDone. See results/enrich_*.csv and *_dotplot.pdf\n")
