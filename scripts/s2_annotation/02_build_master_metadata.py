#!/usr/bin/env python3
"""
04_create_master_metadata.py  (robust, final version)
-------------------------------------------------------
Single master_metadata for the whole project. Reconciles TWO data structures:

  A) per-GSM files in data/raw_counts/ (flattened):  GSE203507, GSE148241
  B) series matrices in data/series_counts/:          GSE306864, GSE190971, GSE186257

+ microarray branch (GSE75010) and cfRNA validation (separate, in cfrna_labels.csv).

Each dataset has an EXPLICIT map: where the diagnosis is, how to read it, where it goes
(branch), which tissue filter and covariates apply. No automatic column searching.

Output: data/master_metadata.csv with columns:
  GSM_ID, dataset, sample_col, raw_label, subtype, category,
  onset, tissue, sex, branch, data_type, matrix_file

Run from project root:
    python scripts/04_create_master_metadata.py
"""

import os
import re
import gzip
import csv
import glob
import pandas as pd
import GEOparse

RAW_DIR = "data/raw_counts"
SERIES_DIR = "data/series_counts"
META_DIR = "data/metadata"
OUT_CSV = "data/master_metadata.csv"

# ── Normalize diagnosis -> subtype + binary category ────────────────
SUBTYPE_RULES = [
    # (regex over raw lowercased label, subtype, category)
    # IMPORTANT: the control guard goes first, otherwise 'non-PE' is incorrectly
    # matched by the PE rule (because it contains the substring 'pe').
    (r"non.?pe|non.?hypertensive|^control$|^normal$|^np$|^tb$",
     "Control", "Control"),
    (r"eope.*fgr|fgr.*eope", "PE_Early_FGR", "Case"),
    (r"eospe|early[- ]onset severe", "PE_Early_Severe", "Case"),
    (r"eope|early[- ]onset", "PE_Early", "Case"),
    (r"severe pre.?eclampsia", "PE_Severe", "Case"),
    (r"pre.?eclampsia|(^|[^a-z])pe([^a-z]|$)", "PE_General", "Case"),
    (r"normal|control", "Control", "Control"),   # catches values like Non-hypertensive Control
]


def map_subtype(raw):
    v = str(raw).lower().strip()
    for pat, sub, cat in SUBTYPE_RULES:
        if re.search(pat, v):
            return sub, cat
    return "Unknown", "Unknown"


# ── Dataset configuration (single source of truth) ───────────────────
# structure: 'per_gsm' | 'series'
# label_col: name of the phenotype column containing diagnosis labels
#            (for per_gsm or series-via-meta datasets)
# join:      how to merge series matrix columns with metadata labels
#            'name_suffix' | 'via_meta' (through which metadata field)
DATASETS = {
    "GSE203507": dict(
        structure="per_gsm", branch="DEA", data_type="raw_counts",
        label_col="characteristics_ch1.1.clinical group",
        tissue_col="characteristics_ch1.0.tissue",
        tissue_keep="villous",                       # keep only villous tissue
        onset="EOPE",
        drop_labels={"fgr", "iptb"},                 # non-hypertensive controls -> drop
    ),
    "GSE148241": dict(
        structure="per_gsm", branch="ML", data_type="rpkm",
        label_col="characteristics_ch1.0.subject status",
        tissue_col="characteristics_ch1.1.tissue",
        tissue_keep="placenta",
        onset="EOPE",
    ),
    "GSE306864": dict(
        structure="series", branch="DEA", data_type="raw_counts",
        matrix="GSE306864_rawCounts.txt.gz",
        join="via_meta", meta_key="description", meta_key_regex=r"((?:SCP|STP)\d+)",
        label_col="characteristics_ch1.2.treatment",
        sex_col="characteristics_ch1.1.genotype",
        onset="LOPE",
    ),
    "GSE190971": dict(
        structure="series", branch="DEA", data_type="raw_counts",
        matrix="GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz",
        join="name_suffix",                          # label encoded in the column name
        onset="mixed",
    ),
    "GSE186257": dict(
        structure="series", branch="ML", data_type="deseq2_norm",
        matrix="GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz",
        join="via_meta", meta_key="title", meta_key_regex=r"(P\d+)",
        label_col="characteristics_ch1.0.disease",
        onset="severe",
    ),
    "GSE75010": dict(
        structure="per_gsm", branch="microarray", data_type="microarray_norm",
        label_col="characteristics_ch1.0.diagnosis",
        tissue_col="characteristics_ch1.1.tissue",
        tissue_keep="placenta",
        onset="mixed",
    ),
}


def load_pheno(acc):
    gse = GEOparse.get_GEO(acc, destdir=META_DIR, include_data=False, silent=True)
    return gse.phenotype_data


def matrix_columns(path):
    with gzip.open(path, "rt") as fh:
        return next(csv.reader(fh, delimiter="\t"))[1:]  # skip gene ID column


# ── Structure handlers ──────────────────────────────────────────────
def handle_per_gsm(acc, cfg, raw_files):
    ph = load_pheno(acc)
    rows = []
    for gsm, r in ph.iterrows():
        raw = r.get(cfg["label_col"], "")
        # tissue filter
        if cfg.get("tissue_col"):
            tis = str(r.get(cfg["tissue_col"], "")).lower()
            if cfg["tissue_keep"] not in tis:
                continue
        # drop labels (for example FGR/iPTB in GSE203507)
        if str(raw).lower().strip() in cfg.get("drop_labels", set()):
            continue
        sub, cat = map_subtype(raw)
        path = next((os.path.join(RAW_DIR, f) for f in raw_files if gsm in f), None)
        rows.append(dict(
            GSM_ID=gsm, dataset=acc, sample_col=gsm, raw_label=raw,
            subtype=sub, category=cat, onset=cfg["onset"],
            tissue=str(r.get(cfg.get("tissue_col", ""), "")),
            sex="", branch=cfg["branch"], data_type=cfg["data_type"],
            matrix_file=path,
        ))
    return rows


def handle_series(acc, cfg):
    ph = load_pheno(acc)
    mpath = os.path.join(SERIES_DIR, cfg["matrix"])
    cols = matrix_columns(mpath)
    rows = []

    if cfg["join"] == "name_suffix":
        # label is directly in the column name (GSE190971: ..._PE / ..._NORMAL)
        for c in cols:
            cat = "Case" if c.upper().endswith("_PE") else "Control"
            sub = "PE_General" if cat == "Case" else "Control"
            rows.append(dict(
                GSM_ID="", dataset=acc, sample_col=c, raw_label=c,
                subtype=sub, category=cat, onset=cfg["onset"], tissue="Placenta",
                sex="", branch=cfg["branch"], data_type=cfg["data_type"],
                matrix_file=mpath,
            ))
        return rows

    if cfg["join"] == "via_meta":
        # build mapping from internal ID -> metadata row
        rx = re.compile(cfg["meta_key_regex"])
        key2meta = {}
        for gsm, r in ph.iterrows():
            m = rx.search(str(r.get(cfg["meta_key"], "")))
            if m:
                key2meta[m.group(1)] = (gsm, r)
        for c in cols:
            m = rx.search(c)
            key = m.group(1) if m else c
            hit = key2meta.get(key)
            if hit is None:
                rows.append(dict(
                    GSM_ID="", dataset=acc, sample_col=c, raw_label="",
                    subtype="Unknown", category="Unknown", onset=cfg["onset"],
                    tissue="", sex="", branch=cfg["branch"],
                    data_type=cfg["data_type"], matrix_file=mpath))
                continue
            gsm, r = hit
            raw = r.get(cfg["label_col"], "")
            sub, cat = map_subtype(raw)
            sex = str(r.get(cfg.get("sex_col", ""), "")) if cfg.get("sex_col") else ""
            rows.append(dict(
                GSM_ID=gsm, dataset=acc, sample_col=c, raw_label=raw,
                subtype=sub, category=cat, onset=cfg["onset"], tissue="",
                sex=sex, branch=cfg["branch"], data_type=cfg["data_type"],
                matrix_file=mpath))
        return rows

    raise ValueError(f"unknown join type: {cfg['join']}")


def main():
    raw_files = [os.path.basename(f) for f in glob.glob(os.path.join(RAW_DIR, "*"))]
    all_rows = []
    for acc, cfg in DATASETS.items():
        if cfg["structure"] == "per_gsm":
            rows = handle_per_gsm(acc, cfg, raw_files)
        else:
            rows = handle_series(acc, cfg)
        all_rows.extend(rows)
        cats = pd.Series([r["category"] for r in rows]).value_counts().to_dict()
        print(f"{acc:12s} [{cfg['branch']:10s}] {len(rows):3d} samples  {cats}")

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)

    print("\n" + "=" * 60)
    print("MASTER METADATA — summary")
    print("=" * 60)
    print("\nBy branch × category:")
    print(pd.crosstab(df["branch"], df["category"]).to_string())
    print("\nDEA branch by onset × category:")
    dea = df[df["branch"] == "DEA"]
    print(pd.crosstab(dea["onset"], dea["category"]).to_string())
    unk = df[df["category"] == "Unknown"]
    if len(unk):
        print(f"\n[!] {len(unk)} samples labeled Unknown — please review:")
        print(unk[["dataset", "sample_col", "raw_label"]].to_string(index=False))
    else:
        print("\n✓ No Unknown labels — all samples are annotated.")
    missing = df[(df["branch"] != "microarray") & (df["matrix_file"].isna())]
    if len(missing):
        print(f"\n[!] {len(missing)} samples missing a matrix file (per-GSM not found).")
    print(f"\nSaved -> {OUT_CSV}  ({len(df)} rows)")


if __name__ == "__main__":
    main()