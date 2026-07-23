#!/usr/bin/env python3
"""
00_inventory_datasets.py
------------------------
Inventory GEO datasets BEFORE building the pipeline.

For each GSE it collects three independent signals:
  1) platform technology (RNA-seq vs microarray)  -> from GPL
  2) data format (raw counts / FPKM / normalized)  -> from the data_processing field in SOFT
  3) integer vs float in actual local files                -> heuristic sniff
Plus: diagnosis columns + their unique values, sample counts,
supplement format (one matrix per series vs file-per-GSM).

It does not reload anything heavy: reads cached .soft.gz from data/metadata.
Run:  python 00_inventory_datasets.py
Output:   console report + data/inventory_summary.csv
"""

import os
import re
import glob
import gzip
import GEOparse
import pandas as pd

# ── Paths (adjust if structure differs) ─────────────────────────────
METADATA_DIR = "data/metadata"
RAW_DIR      = "data/raw_counts"
NORM_DIR     = "data/normalized_data"
OUT_CSV      = "data/inventory_summary.csv"

ACCESSIONS = [
    "GSE203507", "GSE148241", "GSE306864", "GSE190971",
    "GSE186257", "GSE75010", "GSE97320", "GSE154377", "GSE192902",
]

# Keywords used to find the diagnosis column in phenotype_data
DIAG_KEYWORDS = [
    "disease", "clinical", "diagnos", "status", "condition",
    "group", "phenotype", "subtype", "onset", "characteristics",
]


# ── Helper classifiers ──────────────────────────────────────────
def _first(meta, key):
    """GEOparse stores values as lists; safely take the first element."""
    v = meta.get(key, [""])
    return (v[0] if isinstance(v, (list, tuple)) and v else (v or "")) or ""


def classify_platform(gpl):
    tech = _first(gpl.metadata, "technology").lower()
    title = _first(gpl.metadata, "title")
    tl = title.lower()
    if "sequencing" in tech or "seq" in tl:
        kind = "RNA-seq"
    elif "oligonucleotide" in tech or "array" in tech or "beadchip" in tl \
            or "affymetrix" in tl or "genechip" in tl or "array" in tl:
        kind = "microarray"
    else:
        kind = "?"
    return kind, title


def guess_format_from_processing(text):
    """Classify the format based on what the authors wrote in data_processing."""
    t = (text or "").lower()
    if any(k in t for k in ["raw count", "read count", "htseq", "featurecounts",
                            "raw read", "count matrix", "star ", "counts were"]):
        return "RAW COUNTS"
    if any(k in t for k in ["fpkm", "rpkm", "tpm"]):
        return "FPKM/RPKM/TPM (normalized)"
    if any(k in t for k in ["deseq", "tmm", "vst", "rlog", "cpm",
                            "quantile", "rma", "normalized", "log2", "log-transformed"]):
        return "normalized/transformed"
    return "unclear"


def sniff_int_vs_float(paths, gsms):
    """Find the local dataset file and inspect whether it contains integers or floats."""
    matched = None
    for p in paths:
        base = os.path.basename(p)
        if any(g in base for g in gsms):
            matched = p
            break
    if matched is None:
        return "n/a (local file not found)"
    try:
        opener = gzip.open if matched.endswith(".gz") else open
        with opener(matched, "rt", errors="ignore") as fh:
            lines = [next(fh) for _ in range(50)]
    except Exception as e:  # noqa: BLE001
        return f"n/a (could not read: {e})"

    nums = re.findall(r"(?<![A-Za-z_])-?\d+\.?\d*(?:[eE][-+]?\d+)?", " ".join(lines))
    vals = []
    for n in nums:
        try:
            vals.append(float(n))
        except ValueError:
            pass
    vals = [v for v in vals if abs(v) > 1e-9][:200]  # drop zeros/id
    if not vals:
        return "n/a (no numeric values found)"
    frac_non_integer = sum(1 for v in vals if abs(v - round(v)) > 1e-6) / len(vals)
    if frac_non_integer > 0.05:
        return f"FLOAT -> likely normalized ({frac_non_integer:.0%} fractional values)"
    return "INTEGER -> likely raw counts"


def find_diag_columns(pheno):
    cols = [c for c in pheno.columns
            if any(k in c.lower() for k in DIAG_KEYWORDS)]
    out = []
    for c in cols:
        uniq = pheno[c].dropna().astype(str).unique().tolist()
        preview = "; ".join(uniq[:6]) + (" ..." if len(uniq) > 6 else "")
        out.append((c, len(uniq), preview))
    return out


def supp_summary(gse):
    series_supp = gse.metadata.get("supplementary_file", []) or []
    per_gsm = 0
    for gsm in gse.gsms.values():
        if any(k.startswith("supplementary_file") for k in gsm.metadata):
            per_gsm += 1
    if per_gsm > 0:
        return f"per-GSM files ({per_gsm} samples have supp files)"
    if series_supp:
        names = ", ".join(os.path.basename(s) for s in series_supp[:3])
        return f"series-level matrix: {names}"
    return "no supplementary files found in SOFT"


# ── Main loop ────────────────────────────────────────────────────
def inventory():
    local_files = glob.glob(os.path.join(RAW_DIR, "*")) + \
                  glob.glob(os.path.join(NORM_DIR, "*"))
    rows = []

    for acc in ACCESSIONS:
        print("\n" + "=" * 70)
        print(f"  {acc}")
        print("=" * 70)
        try:
            gse = GEOparse.get_GEO(geo=acc, destdir=METADATA_DIR,
                                   include_data=False, silent=True)
        except Exception as e:  # noqa: BLE001
            print(f"  [ERROR] failed to download SOFT: {e}")
            rows.append({"GSE": acc, "platform": "ERROR"})
            continue

        # 1) platform
        plats = [classify_platform(g) for g in gse.gpls.values()]
        kind = "/".join(sorted({k for k, _ in plats})) or "?"
        titles = "; ".join(t for _, t in plats)

        # 2) format from data_processing (take from the first GSM)
        dp = ""
        if gse.gsms:
            first_gsm = next(iter(gse.gsms.values()))
            dp = _first(first_gsm.metadata, "data_processing")
        fmt_proc = guess_format_from_processing(dp)

        # 3) integer vs float in local files
        gsms = list(gse.gsms.keys())
        fmt_sniff = sniff_int_vs_float(local_files, gsms)

        pheno = gse.phenotype_data
        n = len(pheno)
        diag_cols = find_diag_columns(pheno)
        supp = supp_summary(gse)

        # ── print report ──
        print(f"  Technology      : {kind}")
        print(f"  Platform(s)     : {titles}")
        print(f"  Samples (GSM)   : {n}")
        print(f"  Format (SOFT)   : {fmt_proc}")
        print(f"    data_processing (excerpt): {dp[:140]}")
        print(f"  Format (sniff)  : {fmt_sniff}")
        print(f"  Supplementary   : {supp}")
        print(f"  Diagnosis cols :")
        if diag_cols:
            for c, nuniq, preview in diag_cols:
                print(f"     - {c}  ({nuniq} unique): {preview}")
        else:
            print("     [!] no diagnosis columns found by keywords")

        # mismatch flags
        flags = []
        if kind == "microarray":
            flags.append("MICROARRAY -> not suitable for ComBat-seq/DESeq2")
        if "FPKM" in fmt_proc or "normalized" in fmt_proc or "FLOAT" in fmt_sniff:
            flags.append("NOT raw counts -> ML/limma branch only")
        if flags:
            print("  ⚠ WARNING      : " + " | ".join(flags))

        rows.append({
            "GSE": acc,
            "technology": kind,
            "n_samples": n,
            "format_data_processing": fmt_proc,
            "format_sniff": fmt_sniff,
            "supplementary": supp,
            "diag_columns": " || ".join(c for c, _, _ in diag_cols),
            "flags": " | ".join(flags),
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print("\n" + "=" * 70)
    print(f"Summary saved -> {OUT_CSV}")
    print("=" * 70)
    with pd.option_context("display.max_colwidth", 40, "display.width", 200):
        print(df[["GSE", "technology", "n_samples",
                  "format_data_processing", "flags"]].to_string(index=False))


if __name__ == "__main__":
    inventory()