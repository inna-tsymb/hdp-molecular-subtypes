#!/usr/bin/env python3
"""
s5_branch2_ml/01_build_ml_matrix.py
-----------------------------------
Builds the ML training pool (5 RNA-seq datasets) into a normalized
samples x genes matrix for predictive panel selection.

Key principles (leak avoidance):
  * per-dataset NORMALIZATION: log1p -> z-score each gene WITHIN each dataset.
    Makes raw/RPKM/DESeq2-norm comparable, acts like batch correction, and
    prevents information leakage across LODO folds (each dataset is standardized independently).
  * Genes are reduced to a COMMON symbol space (intersection across 5 datasets).
  * Branch 1 DEG lists are NOT included here (to avoid circularity);
    feature selection happens inside CV folds (next script).

Sources:
  GSE203507/306864/190971 — from ready data/dea_counts.csv (symbol, raw)
  GSE186257 — series DESeq2-norm matrix (detect ID -> symbol)
  GSE148241 — per-GSM RPKM files (detect ID -> symbol)

Output:
  data/ml_expr.csv      samples x genes (normalized)
  data/ml_metadata.csv  sample,dataset,category,onset
  data/ml_gene_space.txt

Run from project root:  python scripts/s5_branch2_ml/01_build_ml_matrix.py
"""

import os
import re
import numpy as np
import pandas as pd

MASTER = "data/master_metadata.csv"
DEA_COUNTS = "data/dea_counts.csv"
SERIES_DIR = "data/series_counts"
OUT_EXPR = "data/ml_expr.csv"
OUT_META = "data/ml_metadata.csv"
OUT_GENES = "data/ml_gene_space.txt"

ML_DATASETS = ["GSE203507", "GSE306864", "GSE190971", "GSE186257", "GSE148241"]
ENSEMBL_RE = re.compile(r"^ENSG\d+")


def looks_ensembl(index):
    idx = list(index)
    if not idx:
        return False
    frac = sum(bool(ENSEMBL_RE.match(str(i))) for i in idx) / len(idx)
    return frac > 0.3          # fraction of ENSG IDs across the whole index, not just the first 50


def ensembl_to_symbol(ids):
    import mygene
    mg = mygene.MyGeneInfo()
    r = mg.querymany([i.split(".")[0] for i in ids], scopes="ensembl.gene",
                     fields="symbol", species="human",
                     as_dataframe=True, verbose=False)
    r = r[~r.index.duplicated(keep="first")]
    if "symbol" not in r.columns:
        print("   [!] mygene did not map any IDs — check the Ensembl format")
        return {}
    return r["symbol"].dropna().to_dict()


def to_symbol(df):
    df.index = df.index.astype(str).str.strip().str.strip('"').str.strip()
    # composite IDs like ENSG00000000003$TSPAN6 -> keep the symbol part (after $ or |)
    if df.index.str.contains(r"[$|]").mean() > 0.5:
        df.index = df.index.str.split(r"[$|]").str[-1]
    df = df[~df.index.str.startswith("__")]
    if looks_ensembl(df.index):
        df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
        m = ensembl_to_symbol(list(df.index))
        n_before = len(df)
        df.index = df.index.map(lambda e: m.get(e))
        df = df[df.index.notna()]
        print(f"   Ensembl->symbol: {len(df)}/{n_before} mapped")
    return df.groupby(level=0).sum()


def load_from_dea(dataset, meta):
    dea = pd.read_csv(DEA_COUNTS, index_col=0)
    cols = [c for c in meta[meta.dataset == dataset]["sample_col"] if c in dea.columns]
    return dea[cols]


def load_series_186257(meta):
    f = os.path.join(SERIES_DIR, "GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz")
    df = to_symbol(pd.read_csv(f, sep="\t", index_col=0, compression="gzip"))
    keep = [c for c in df.columns
            if c in set(meta[meta.dataset == "GSE186257"]["sample_col"])]
    return df[keep]


def load_per_gsm_148241(meta):
    cols = {}
    for _, r in meta[meta.dataset == "GSE148241"].iterrows():
        p = r["matrix_file"]
        if not isinstance(p, str) or not os.path.exists(p):
            continue
        # files: header ENSEMBL_ID\traw.counts, 2 columns, Ensembl without versions
        s = pd.read_csv(p, sep="\t", index_col=0,
                        compression="gzip" if p.endswith(".gz") else None,
                        usecols=[0, 1], low_memory=False)
        vals = pd.to_numeric(s.iloc[:, 0], errors="coerce")
        # clean index: remove NaN/empty values and collapse duplicates
        vals.index = vals.index.astype(str).str.strip()
        vals = vals[vals.index.notna() & (vals.index != "") & (vals.index != "nan")]
        vals = vals.groupby(level=0).sum()
        cols[r["sample_col"]] = vals
    # assemble by outer join reindex over all unique genes
    all_genes = sorted(set().union(*[set(v.index) for v in cols.values()]))
    df = pd.DataFrame({k: v.reindex(all_genes) for k, v in cols.items()},
                      index=all_genes).fillna(0)
    return to_symbol(df)


def norm_per_dataset(df):
    """genes x samples -> log1p -> z-score each gene (row)."""
    x = np.log1p(df.astype(float))
    mu, sd = x.mean(axis=1), x.std(axis=1).replace(0, np.nan)
    return x.sub(mu, axis=0).div(sd, axis=0).dropna(how="all")


def main():
    master = pd.read_csv(MASTER)
    ml = master[master.dataset.isin(ML_DATASETS)].copy()
    print("ML samples in master:", len(ml))

    frames = {
        "GSE203507": load_from_dea("GSE203507", ml),
        "GSE306864": load_from_dea("GSE306864", ml),
        "GSE190971": load_from_dea("GSE190971", ml),
    }
    print("Loading GSE186257 (series, DESeq2-norm)...")
    frames["GSE186257"] = load_series_186257(ml)
    print("Loading GSE148241 (per-GSM, RPKM)...")
    frames["GSE148241"] = load_per_gsm_148241(ml)

    for k, v in frames.items():
        print(f"  {k}: {v.shape[0]} genes x {v.shape[1]} samples")

    common = sorted(set.intersection(*[set(f.index) for f in frames.values()]))
    print(f"\nShared genes (5 datasets): {len(common)}")

    parts = [norm_per_dataset(f.loc[common]).T for f in frames.values()]
    expr = pd.concat(parts, axis=0).dropna(axis=1)

    meta = (ml.set_index("sample_col")
            .loc[expr.index, ["dataset", "category", "onset"]]
            .rename_axis("sample").reset_index())
    assert list(meta["sample"]) == list(expr.index)

    expr.to_csv(OUT_EXPR)
    meta.to_csv(OUT_META, index=False)
    with open(OUT_GENES, "w") as fh:
        fh.write("\n".join(expr.columns))

    print("\n" + "=" * 55)
    print(f"ML matrix: {expr.shape[0]} samples x {expr.shape[1]} genes")
    print("\nSamples by dataset × category:")
    print(pd.crosstab(meta["dataset"], meta["category"]).to_string())
    print(f"\nSaved:\n  {OUT_EXPR}\n  {OUT_META}\n  {OUT_GENES}")


if __name__ == "__main__":
    main()