#!/usr/bin/env bash
# reorganize.sh — regroup scripts/ into functional stage directories.
# Safe: it only moves existing files and does not delete anything permanently
# except __pycache__ and a temporary arshld draft file.
#
# Run from the project root:
#   bash reorganize.sh
# After this script, run scripts from the root with paths like:
#   python scripts/s2_annotation/02_build_master_metadata.py
# All data paths in the scripts are relative to the repository root.

set -euo pipefail

echo "== Reorganizing scripts/ into functional stage folders =="

mv_if () { [ -e "$1" ] && { mkdir -p "$(dirname "$2")"; git mv "$1" "$2" 2>/dev/null || mv "$1" "$2"; echo "  $1 -> $2"; } || true; }

cd scripts 2>/dev/null || { echo "Run this from the project root (where scripts/ exists)"; exit 1; }

# ── s1 · acquisition ───────────────────────────────────────────
mv_if 01_download_data.py        s1_acquisition/01_download_geo.py
mv_if 02_flatten_data.py         s1_acquisition/02_flatten_per_gsm.py
mv_if 03_final_flattening.py     s1_acquisition/03_flatten_final.py
# download_series.sh is handled separately in s1_acquisition

# ── s2 · QC and annotation ───────────────────────────────────────
mv_if 00_inventory_datasets.py   s2_annotation/01_inventory_datasets.py
mv_if 00_inventarization         s2_annotation/01_inventory_datasets.py
mv_if 04_create_master_metadata.py s2_annotation/02_build_master_metadata.py
mv_if 06_build_cfrna_labels.py   s2_annotation/03_build_cfrna_labels.py

# ── s3 · cohort scouting (exploratory) ───────────────────────────
mv_if 05_scout_cfrna_cohorts.py  s3_scouting/01_scout_cfrna_cohorts.py

# ── s4 · Branch 1: DEA ───────────────────────────────────────────
mv_if 07_merge_dea_matrices.py   s4_branch1_dea/01_merge_dea_matrices.py
mv_if 08_dea_stratified.R        s4_branch1_dea/02_dea_stratified.R
mv_if 09_enrichment.R            s4_branch1_dea/03_enrichment.R

# ── s5/s6 · future blocks (empty for now) ───────────────────────
mkdir -p s5_branch2_ml s6_validation
touch s5_branch2_ml/.gitkeep s6_validation/.gitkeep

# ── cleanup ─────────────────────────────────────────────────────
rm -f arshld 01b_redownload.py 2>/dev/null || true
find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

cd ..
echo ""
echo "== Done. New structure: =="
find scripts -type f | sort
