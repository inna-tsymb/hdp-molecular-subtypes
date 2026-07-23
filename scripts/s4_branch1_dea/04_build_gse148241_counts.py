#!/usr/bin/env python3
"""
s4_branch1_dea/04_build_gse148241_counts.py
-------------------------------------------
Collects RAW counts for GSE148241 (EOSPE, placenta-only) into a
genes(symbol) x samples matrix for independent replication of the EOPE signature.

GSE148241 per-GSM files: header ENSEMBL_ID\\traw.counts, integer counts.
-> suitable for DESeq2 (unlike the normalized series matrix).

Output:
  data/gse148241_counts.csv   genes(symbol) x samples  (raw int)
  data/gse148241_meta.csv     sample, category

Run from project root:  python scripts/s4_branch1_dea/04_build_gse148241_counts.py
"""

import os
import pandas as pd

MASTER = "data/master_metadata.csv"


def ensembl_to_symbol(ids):
    import mygene
    mg = mygene.MyGeneInfo()
    r = mg.querymany([i.split(".")[0] for i in ids], scopes="ensembl.gene",
                     fields="symbol", species="human",
                     as_dataframe=True, verbose=False)
    r = r[~r.index.duplicated(keep="first")]
    return r["symbol"].dropna().to_dict() if "symbol" in r.columns else {}


def main():
    m = pd.read_csv(MASTER)
    rows = m[m.dataset == "GSE148241"]
    print(f"GSE148241 samples in master (placenta): {len(rows)}")

    import glob
    # GSE148241 has multiple supplementary files per GSM (raw.counts, editing.ratio, ...);
    # take ONLY raw.counts by matching GSM in raw_counts/
    raw_files = glob.glob("data/raw_counts/*raw.counts*")
    cols = {}
    for _, r in rows.iterrows():
        gsm = r["sample_col"]
        cand = [f for f in raw_files if gsm in f]
        if not cand:
            print(f"   [!] missing raw.counts for {gsm}")
            continue
        p = cand[0]
        s = pd.read_csv(p, sep="\t", index_col=0,
                        compression="gzip" if p.endswith(".gz") else None,
                        usecols=[0, 1], low_memory=False)
        vals = pd.to_numeric(s.iloc[:, 0], errors="coerce")
        vals.index = vals.index.astype(str).str.strip().str.replace(r"\.\d+$", "", regex=True)
        vals = vals[vals.index.notna() & (vals.index != "") & (vals.index != "nan")]
        vals = vals.groupby(level=0).sum()
        cols[r["sample_col"]] = vals

    all_ens = sorted(set().union(*[set(v.index) for v in cols.values()]))
    df = pd.DataFrame({k: v.reindex(all_ens) for k, v in cols.items()},
                      index=all_ens).fillna(0)
    print(f"Collected (Ensembl): {df.shape[0]} genes x {df.shape[1]} samples")

    # Ensembl -> symbol
    print("Mapping Ensembl -> symbol ...")
    e2s = ensembl_to_symbol(list(df.index))
    df.index = df.index.map(lambda e: e2s.get(e))
    df = df[df.index.notna()]
    df = df.groupby(level=0).sum()          # collapse duplicate symbols (sum -> int)
    df = df.round().astype(int)
    print(f"In symbols: {df.shape[0]} genes x {df.shape[1]} samples")

    meta = (rows.set_index("sample_col").loc[df.columns, ["category"]]
            .rename_axis("sample").reset_index())
    assert list(meta["sample"]) == list(df.columns)

    df.to_csv("data/gse148241_counts.csv")
    meta.to_csv("data/gse148241_meta.csv", index=False)
    print("\nBalance:", meta["category"].value_counts().to_dict())
    print("Saved: data/gse148241_counts.csv, data/gse148241_meta.csv")


if __name__ == "__main__":
    main()
