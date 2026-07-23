#!/usr/bin/env python3
"""
07_merge_dea_matrices.py
------------------------
Collects THREE DEA datasets with different structures and ID spaces into a SINGLE
raw count matrix genes(symbol) x 91, aligned with master_metadata,
and prepares input for the R script (ComBat_seq -> DESeq2).

Structures handled:
  GSE203507 — 20 per-GSM files (symbol in quotes, 2 columns)   -> collect
  GSE306864 — series matrix, Ensembl ID, columns SCP/STP         -> Ensembl->symbol
  GSE190971 — series matrix, gene symbol, columns V..._PE/NORMAL

Target ID space: GENE SYMBOL (2 of 3 are already there; only GSE306864 is mapped).
HTSeq/STAR tails (__no_feature etc.) are trimmed.

Output:
  data/dea_counts.csv     genes x 91  (raw integer counts)
  data/dea_metadata.csv   91 x [sample,dataset,category,onset,sex]  (in column order)
  data/master_gene_list.txt

Run from project root:
    python scripts/07_merge_dea_matrices.py
Requires: pip install mygene   (for Ensembl->symbol)
"""

import os
import re
import gzip
import glob
import sys
import pandas as pd

MASTER = "data/master_metadata.csv"
OUT_COUNTS = "data/dea_counts.csv"
OUT_META = "data/dea_metadata.csv"
OUT_GENES = "data/master_gene_list.txt"


def strip_htseq_tails(df):
    mask = ~df.index.astype(str).str.startswith("__")
    return df[mask]


def collapse_symbols(df):
    """Sum counts for duplicate gene symbols (after mapping, duplicates may occur)."""
    df = df[df.index.notna()]
    df.index = df.index.astype(str)
    return df.groupby(level=0).sum()


# ── Dataset loaders ────────────────────────────────────────
def load_gse203507(meta_rows):
    """20 per-GSM files -> genes(symbol) x samples."""
    series = {}
    for _, r in meta_rows.iterrows():
        path = r["matrix_file"]
        if not isinstance(path, str) or not os.path.exists(path):
            print(f"   [!] missing file for {r['sample_col']}: {path}")
            continue
        s = pd.read_csv(path, sep="\t", header=None, index_col=0,
                        compression="gzip")
        s.index = s.index.astype(str).str.strip('"')      # strip quotes
        series[r["sample_col"]] = s.iloc[:, 0]
    df = pd.DataFrame(series)
    return collapse_symbols(strip_htseq_tails(df))


def load_gse190971(path):
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return collapse_symbols(strip_htseq_tails(df))


def load_gse306864(path, ens2sym):
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    df = strip_htseq_tails(df)
    # Ensembl IDs may have version suffixes (ENSG....1) — strip them
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    df.index = df.index.map(lambda e: ens2sym.get(e))     # -> symbol or None
    return collapse_symbols(df)


def ensembl_to_symbol(ensembl_ids):
    try:
        import mygene
    except ImportError:
        sys.exit("mygene is required: pip install mygene")
    mg = mygene.MyGeneInfo()
    res = mg.querymany(ensembl_ids, scopes="ensembl.gene",
                       fields="symbol", species="human",
                       as_dataframe=True, verbose=False)
    res = res[~res.index.duplicated(keep="first")]
    return res["symbol"].dropna().to_dict()


# ── Main ──────────────────────────────────────────────────────────
def main():
    master = pd.read_csv(MASTER)
    dea = master[master["branch"] == "DEA"].copy()
    print(f"DEA samples in master: {len(dea)}")

    frames = {}

    # GSE203507
    m = dea[dea.dataset == "GSE203507"]
    frames["GSE203507"] = load_gse203507(m)
    print(f"GSE203507: {frames['GSE203507'].shape} (genes x samples)")

    # GSE190971
    p = dea[dea.dataset == "GSE190971"]["matrix_file"].dropna().iloc[0]
    frames["GSE190971"] = load_gse190971(p)
    print(f"GSE190971: {frames['GSE190971'].shape}")

    # GSE306864 — first map Ensembl->symbol
    p = dea[dea.dataset == "GSE306864"]["matrix_file"].dropna().iloc[0]
    ens_ids = pd.read_csv(p, sep="\t", index_col=0, compression="gzip",
                          usecols=[0]).index.astype(str)
    ens_ids = [e.split(".")[0] for e in ens_ids if not e.startswith("__")]
    print(f"GSE306864: mapping {len(ens_ids)} Ensembl IDs -> symbol ...")
    ens2sym = ensembl_to_symbol(ens_ids)
    print(f"   mapped: {len(ens2sym)}/{len(ens_ids)}")
    frames["GSE306864"] = load_gse306864(p, ens2sym)
    print(f"GSE306864: {frames['GSE306864'].shape}")

    # ── master_gene_list = intersection of symbols ──
    gene_sets = [set(f.index) for f in frames.values()]
    common = sorted(set.intersection(*gene_sets))
    print(f"\nShared genes (master_gene_list): {len(common)}")
    for name, f in frames.items():
        print(f"   {name}: {f.shape[0]} genes, {len(set(f.index) & set(common))} in intersection")

    # ── assemble a single genes x samples matrix ──
    parts = [f.loc[common] for f in frames.values()]
    counts = pd.concat(parts, axis=1)
    counts = counts.round().astype(int)          # ComBat_seq requires int

    # ── align metadata with matrix columns ──
    meta = (dea.set_index("sample_col")
            .loc[counts.columns, ["dataset", "category", "onset", "sex"]]
            .rename_axis("sample")
            .reset_index())

    assert list(meta["sample"]) == list(counts.columns), "mismatch!"

    counts.to_csv(OUT_COUNTS)
    meta.to_csv(OUT_META, index=False)
    with open(OUT_GENES, "w") as fh:
        fh.write("\n".join(common))

    print("\n" + "=" * 55)
    print(f"Matrix: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print("Sample distribution (dataset x category):")
    print(pd.crosstab(meta["dataset"], meta["category"]).to_string())
    print(f"\nSaved:\n  {OUT_COUNTS}\n  {OUT_META}\n  {OUT_GENES}")


if __name__ == "__main__":
    main()
