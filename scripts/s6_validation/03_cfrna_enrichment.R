#!/usr/bin/env Rscript
# s6_validation/07_cfrna_enrichment.R
# -------------------------------------------------------------------
# Functional enrichment of the NATIVE cfRNA PE signature (13-20 wk, 76 DEG).
# Question: which biology is active in plasma and how does it differ from placenta
# (hypoxia/HIF-1 for EOPE)? Provides symmetry: "placenta = X, plasma = Y".
#
# Universe = genes tested in the same subset (not the whole genome!).
# Additionally: check against known cfRNA PE markers from literature
# (Del Vecchio 2020: S100A8, MS4A3, MMP8, BCL2L15, ALPL).
#
# Input:  results/cfrna_deg_mid_13_20.csv, data/cfrna_dea/counts_mid_13_20.csv
# Output: results/enrich_cfrna_mid_{GO,KEGG}.csv + dotplot, cfrna_biology.txt
#
# Run:  Rscript scripts/s6_validation/07_cfrna_enrichment.R
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2)
})
dir.create("results", showWarnings = FALSE)

sig <- read.csv("results/cfrna_deg_mid_13_20.csv")
genes <- sig$gene
# universe = genes that passed the filter in the same subset
cnt <- read.csv("data/cfrna_dea/counts_mid_13_20.csv", row.names = 1, check.names = FALSE)
universe <- rownames(cnt)
cat("cfRNA-DEG (13-20 wk):", length(genes), "| universe:", length(universe), "\n")

up   <- sig$gene[sig$log2FoldChange > 0]
down <- sig$gene[sig$log2FoldChange < 0]
cat("Up in PE:", length(up), "| Down in PE:", length(down), "\n")

save_e <- function(obj, tag) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    cat("  ", tag, ": empty\n"); return(invisible(NULL))
  }
  df <- as.data.frame(obj)
  write.csv(df, sprintf("results/enrich_cfrna_%s.csv", tag), row.names = FALSE)
  n <- min(15, nrow(df))
  ggsave(sprintf("results/enrich_cfrna_%s_dotplot.pdf", tag),
         dotplot(obj, showCategory = n) + ggtitle(paste("cfRNA 13-20wk:", tag)),
         width = 8, height = max(4, n * 0.35))
  cat("  ", tag, ":", nrow(df), "terms | top:",
      paste(head(df$Description, 3), collapse = "; "), "\n")
  df$Description
}

go <- tryCatch(enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                        universe = universe, pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, qvalueCutoff = 0.10),
               error = function(e) NULL)
go_terms <- save_e(go, "mid_GO")

kegg <- tryCatch({
  ent <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  uni <- bitr(universe, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  enrichKEGG(ent, organism = "hsa", universe = uni, pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
}, error = function(e) { message("KEGG: ", conditionMessage(e)); NULL })
kegg_terms <- save_e(kegg, "mid_KEGG")

# ── compare with placental biology ──────────────────────────────
plac_go <- if (file.exists("results/enrich_EOPE_GO.csv"))
  read.csv("results/enrich_EOPE_GO.csv")$Description else character(0)
shared_terms <- intersect(c(go_terms, kegg_terms), plac_go)
hyp <- grep("hypox|oxygen|HIF", c(go_terms, kegg_terms), ignore.case = TRUE, value = TRUE)

# ── check against literature cfRNA markers ───────────────────────────
lit <- c("S100A8", "MS4A3", "MMP8", "BCL2L15", "ALPL")   # Del Vecchio 2020
found <- intersect(genes, lit)

sink("results/cfrna_biology.txt")
cat("NATIVE cfRNA PE SIGNATURE (13-20 wk)\n")
cat("=======================================\n")
cat("DEG:", length(genes), " (up:", length(up), " down:", length(down), ")\n\n")
cat("GO BP terms:", length(go_terms), "\n")
if (length(go_terms)) cat(paste0("  - ", head(go_terms, 15), "\n"))
cat("\nKEGG pathways:", length(kegg_terms), "\n")
if (length(kegg_terms)) cat(paste0("  - ", head(kegg_terms, 10), "\n"))
cat("\nShared terms with placental EOPE enrichment:",
    if (length(shared_terms)) paste(shared_terms, collapse = "; ") else "NONE", "\n")
cat("Hypoxia/HIF terms in plasma:",
    if (length(hyp)) paste(hyp, collapse = "; ") else "NONE", "\n\n")
cat("Literature cfRNA PE markers (Del Vecchio 2020) among our DEG:",
    if (length(found)) paste(found, collapse = ", ") else "none", "\n\n")
cat("TOP-15 DEG by padj:\n")
print(head(sig[order(sig$padj), c("gene", "log2FoldChange", "padj")], 15))
sink()

cat("\nShared terms with placenta:", length(shared_terms),
    "| hypoxia terms in plasma:", length(hyp),
    "| lit. markers found:", length(found), "\n")
cat("Saved: results/enrich_cfrna_*.csv, cfrna_biology.txt\n")
