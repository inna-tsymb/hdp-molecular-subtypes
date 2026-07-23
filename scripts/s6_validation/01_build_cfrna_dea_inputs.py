#!/usr/bin/env python3
"""
s6_validation/05_build_cfrna_dea_inputs.py
------------------------------------------
Prepares inputs for cfRNA DEA (Moufarrej): raw counts + coldata for subsets.
Each subset uses ONE sample per patient (avoids longitudinal pseudo-replication).

Subsets:
  early_le12, mid_13_20, late_ge23   — by timepoint
  pooled_antenatal                    — earliest antenatal sample per patient

Output: data/cfrna_dea/counts_<name>.csv, coldata_<name>.csv, subsets.txt

Run from project root:  python scripts/s6_validation/05_build_cfrna_dea_inputs.py
"""

import os
import glob
import numpy as np
import pandas as pd

OUT = "data/cfrna_dea"
os.makedirs(OUT, exist_ok=True)

# ── Moufarrej counts ─────────────────────────────────────────────────
mats = []
for f in sorted(glob.glob("data/moufarrej/GSE192902_counts_*postQC.csv.gz")):
    df = pd.read_csv(f, index_col=0)
    df = df.drop(columns=[c for c in ["gene_num"] if c in df.columns])
    mats.append(df)
counts = pd.concat(mats, axis=1)
counts = counts[~counts.index.duplicated(keep="first")]
counts.columns = [str(c).split(".")[0] for c in counts.columns]
counts = counts.loc[:, ~counts.columns.duplicated(keep="first")]

lab = pd.read_csv("data/cfrna_labels.csv")
lab["sample"] = lab["sample"].astype(str)
lab = lab.drop_duplicates("sample").set_index("sample")
common = [c for c in counts.columns if c in lab.index]
counts = counts[common]
info = lab.loc[common].copy()
info["condition"] = np.where(info["label"] == "PE", "PE", "control")
print(f"Moufarrej: {counts.shape[0]} genes x {counts.shape[1]} samples")

TP = {
    "early_le12": ["≤12 weeks gestation", "<=12 weeks gestation"],
    "mid_13_20": ["13-20 weeks gestation"],
    "late_ge23": ["≥23 weeks gestation"],
}
TP_ORDER = ["≤12 weeks gestation", "<=12 weeks gestation",
            "13-20 weeks gestation", "≥23 weeks gestation"]


def emit(name, sub):
    sub = sub[~sub["patient"].duplicated(keep="first")]        # one sample per patient
    if sub["condition"].nunique() < 2 or len(sub) < 10:
        print(f"  {name}: skip (n={len(sub)}, classes={sub['condition'].nunique()})")
        return None
    cd = counts[sub.index].round().astype(int)
    cd.to_csv(f"{OUT}/counts_{name}.csv")
    sub[["condition", "cohort", "timepoint", "patient"]].to_csv(f"{OUT}/coldata_{name}.csv")
    bal = sub["condition"].value_counts().to_dict()
    coh = sub["cohort"].nunique()
    print(f"  {name}: {cd.shape[1]} samples {bal}  cohorts={coh}")
    return name


subsets = []
for name, tps in TP.items():
    s = info[info["timepoint"].isin(tps)]
    if emit(name, s):
        subsets.append(name)

# pooled: earliest antenatal sample per patient
ant = info[info["timepoint"].isin(sum(TP.values(), []))].copy()
ant["tp_rank"] = ant["timepoint"].map({t: i for i, t in enumerate(TP_ORDER)})
ant = ant.sort_values("tp_rank")
if emit("pooled_antenatal", ant):
    subsets.append("pooled_antenatal")

with open(f"{OUT}/subsets.txt", "w") as fh:
    fh.write("\n".join(subsets))
print(f"\nSubsets: {subsets}")
print(f"Saved in {OUT}/")
