#!/usr/bin/env python3
"""
06_build_cfrna_labels.py
------------------------
Builds a single label table for cfRNA validation (Moufarrej / GSE192902),
checks that counts align with labels across ALL three matrices, and
writes the ready-to-use cfrna_labels.csv — the single source of truth for validation.

Run from project root:
    python 06_build_cfrna_labels.py
"""

import os
import glob
import pandas as pd
import GEOparse

MOUF_DIR = "data/moufarrej"
META_DIR = "data/metadata"
OUT_CSV = "data/cfrna_labels.csv"

# matrix file -> cohort name mapping (for matrix_file)
MATRICES = {
    "GSE192902_counts_Discovery_postQC.csv.gz":   "Discovery",
    "GSE192902_counts_Validation1_postQC.csv.gz": "Validation 1",
    "GSE192902_counts_Validation2_postQC.csv.gz": "Validation 2",
}


def load_labels():
    gse = GEOparse.get_GEO("GSE192902", destdir=META_DIR,
                           include_data=False, silent=True)
    ph = gse.phenotype_data
    lab = ph[[
        "title",
        "characteristics_ch1.0.disease",
        "characteristics_ch1.1.sampling time group",
        "characteristics_ch1.2.cohort",
    ]].copy()
    lab.columns = ["sample", "disease", "timepoint", "cohort"]
    lab["patient"] = lab["sample"].str.extract(r"^(\d+)_")
    lab["label"] = lab["disease"].map(
        lambda d: "control" if d == "control" else "PE")
    lab["is_severe"] = lab["disease"].eq("severe pre-eclampsia")
    return lab.reset_index(drop=True)


def header_samples(path):
    """Extract sample names from the matrix header (skip gene_name, gene_num)."""
    cols = pd.read_csv(path, nrows=0).columns.tolist()[2:]
    clean = [c.split(".")[0] for c in cols]  # remove pandas suffix .1
    return cols, clean


def main():
    lab = load_labels()
    known = set(lab["sample"])

    print("=" * 60)
    print("Samples by diagnosis:")
    print(lab["disease"].value_counts().to_string())
    print("\nUnique patients (not samples):")
    print(lab.groupby("patient")["label"].first().value_counts().to_string())
    # merge information for each matrix
    print("\n" + "=" * 60)
    matrix_of = {}   # sample -> matrix file
    for fname, coh in MATRICES.items():
        path = os.path.join(MOUF_DIR, fname)
        if not os.path.exists(path):
            print(f"[!] missing file: {path} — skipping")
            continue
        raw_cols, clean = header_samples(path)
        matched, unmatched = [], []
        for c, cc in zip(raw_cols, clean):
            (matched if cc in known else unmatched).append((c, cc))
        for c, cc in matched:
            matrix_of[cc] = fname
        print(f"\n{fname}  [{coh}]")
        print(f"  columns: {len(raw_cols)} | matched: {len(matched)} | "
              f"unmatched: {len(unmatched)}")
        if unmatched:
            dupes = [c for c, _ in unmatched if "." in c]
            real = [c for c, _ in unmatched if "." not in c]
            if dupes:
                print(f"    - {len(dupes)} duplicates (.N suffix) -> ignored: "
                      f"{dupes[:5]}{' ...' if len(dupes) > 5 else ''}")
            if real:
                print(f"    - {len(real)} unlabeled columns (dropped): {real}")

    # final table: only samples present in both labels and matrices
    lab["matrix_file"] = lab["sample"].map(matrix_of)
    final = lab.dropna(subset=["matrix_file"]).copy()

    print("\n" + "=" * 60)
    print(f"Valid samples (label + matrix): {len(final)}")
    print("\nPatients × timepoint (unique patient count):")
    piv = (final.groupby(["timepoint", "label"])["patient"]
           .nunique().unstack(fill_value=0))
    print(piv.to_string())

    final.to_csv(OUT_CSV, index=False)
    print(f"\nSaved -> {OUT_CSV}")
    print("Columns:", list(final.columns))


if __name__ == "__main__":
    main()
