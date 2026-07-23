#!/usr/bin/env bash
# s1_acquisition/download_series.sh
# ---------------------------------------------------------------
# Downloads series count matrices that GEOparse skips because those files
# are stored as series-level supplements instead of per-GSM files.
# Also downloads cfRNA Moufarrej matrices for validation.
#
# Run from the project root:
#   bash scripts/s1_acquisition/download_series.sh
# ---------------------------------------------------------------
set -euo pipefail
B="https://ftp.ncbi.nlm.nih.gov/geo/series"

mkdir -p data/series_counts
echo "== Series count matrices (Branch 1 DEA) =="
curl -sS -o data/series_counts/GSE306864_rawCounts.txt.gz \
  "$B/GSE306nnn/GSE306864/suppl/GSE306864_rawCounts.txt.gz"
curl -sS -o data/series_counts/GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz \
  "$B/GSE190nnn/GSE190971/suppl/GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz"
# GSE186257 is DESeq2-normalized (ML branch, not DEA):
curl -sS -o data/series_counts/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz \
  "$B/GSE186nnn/GSE186257/suppl/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz"

mkdir -p data/moufarrej
echo "== cfRNA validation matrices (Moufarrej GSE192902) =="
for f in Discovery Validation1 Validation2; do
  curl -sS -o "data/moufarrej/GSE192902_counts_${f}_postQC.csv.gz" \
    "$B/GSE192nnn/GSE192902/suppl/GSE192902_counts_${f}_postQC.csv.gz"
done

echo "== Done =="
ls -lh data/series_counts data/moufarrej
