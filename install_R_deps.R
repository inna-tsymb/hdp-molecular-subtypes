#!/usr/bin/env Rscript
# install_R_deps.R — install R dependencies for Branch 1.
# Run: Rscript install_R_deps.R
# After installation, capture the session info to a file for reproducibility.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc <- c("DESeq2", "apeglm", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
cran <- c("ggplot2")

for (p in cran)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

BiocManager::install(setdiff(bioc, rownames(installed.packages())),
                     update = FALSE, ask = FALSE)

cat("\nInstalled R dependencies. Verifying load:\n")
invisible(lapply(c(bioc, cran), function(p)
  cat(" ", p, ":", requireNamespace(p, quietly = TRUE), "\n")))

writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")
cat("\nsessionInfo() -> R_sessionInfo.txt\n")
