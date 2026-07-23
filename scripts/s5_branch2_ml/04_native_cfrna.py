#!/usr/bin/env python3
"""
s5_branch2_ml/04_native_cfrna.py
--------------------------------
Native cfRNA PE predictor: the model is TRAINED on plasma (LODO across three
Moufarrej cohorts), not transferred from placenta. Provides a contrast:
   transfer (placenta->cfRNA, weak)  vs native cfRNA (model sees plasma).

Design (leak control, same as the panel):
  * 3 Moufarrej cohorts (Discovery/Val1/Val2) = 3 LODO folds.
  * Per-cohort normalization (log+z) — the held-out cohort is standardized on its own.
  * Feature selection (pre-filter + LASSO) INSIDE the fold.
  * Evaluation at the PATIENT level (one sample per patient within a timepoint).
  * Pooled out-of-fold AUC + bootstrap CI (resample patients).
  * Feature spaces: transferable (6184, main contrast) + full Moufarrej (ceiling).
  * Focus: early plasma (<=20 weeks). Co-temporal (>=23) is secondary (only 2 cohorts).

Output: results/native_cfrna.csv

Run from project root:  python scripts/s5_branch2_ml/04_native_cfrna.py
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)
RNG = np.random.default_rng(0)
B = 2000
K_PREFILTER = 500
C_GRID = np.logspace(-2.0, 0.3, 8)

TRANSFER_REF = {"early (<=20wk)": 0.657, "co-temporal (>=23wk)": 0.544}  # transfer top25 reference

# ── Moufarrej: count matrix + labels ───────────────────────────────
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
info = lab.loc[common, ["patient", "label", "timepoint", "cohort"]].copy()
info["y"] = (info["label"] == "PE").astype(int)
print(f"Moufarrej: {counts.shape[0]} genes x {counts.shape[1]} samples")

FEATURE_SPACES = {
    "transferable": [g for g in
                     open("data/transferable_genes.txt").read().split()
                     if g in counts.index],
    "full": list(counts.index),
}
TP_GROUPS = {
    "early (<=20wk)": ["<=12 weeks gestation", "13-20 weeks gestation",
                       "≤12 weeks gestation"],
    "co-temporal (>=23wk)": ["≥23 weeks gestation"],
}


def make_pipe(C):
    return Pipeline([
        ("sel", SelectKBest(f_classif, k=K_PREFILTER)),
        ("clf", LogisticRegression(penalty="l1", solver="saga", C=C,
                                   class_weight="balanced", max_iter=10000,
                                   random_state=0)),
    ])


def fit_l1(Xtr, ytr):
    inner = StratifiedKFold(4, shuffle=True, random_state=0)
    best_C, best_s = C_GRID[0], -1.0
    for C in C_GRID:
        s = cross_val_score(make_pipe(C), Xtr, ytr, cv=inner,
                            scoring="roc_auc").mean()
        if s > best_s:
            best_s, best_C = s, C
    p = make_pipe(best_C)
    p.fit(Xtr, ytr)
    return p


def norm_cohort(sub):
    """genes x samples -> log+z by gene within the cohort."""
    x = np.log1p(sub.astype(float))
    return x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1).replace(0, np.nan), axis=0)


def boot_auc(y, s, groups):
    base = roc_auc_score(y, s)
    uniq = np.array(sorted(set(groups)))
    g2i = {g: np.where(groups == g)[0] for g in uniq}
    aucs = []
    for _ in range(B):
        samp = RNG.choice(uniq, len(uniq), replace=True)
        idx = np.concatenate([g2i[g] for g in samp])
        if len(set(y[idx])) == 2:
            aucs.append(roc_auc_score(y[idx], s[idx]))
    lo, hi = np.percentile(aucs, [2.5, 97.5]) if aucs else (np.nan, np.nan)
    return base, lo, hi


rows = []
for space_name, genes in FEATURE_SPACES.items():
    genes = [g for g in genes if g in counts.index]
    for tp_name, tps in TP_GROUPS.items():
        sel = info[info["timepoint"].isin(tps)].copy()
        # one sample per patient (keep sample index)
        sel = sel[~sel["patient"].duplicated(keep="first")]
        cohorts = [c for c in ["Discovery", "Validation 1", "Validation 2"]
                   if (sel["cohort"] == c).sum() >= 5]
        if len(cohorts) < 2:
            continue
        oof_prob, oof_y, oof_pat = [], [], []
        for held in cohorts:
            tr_idx = sel.index[sel["cohort"] != held]
            te_idx = sel.index[sel["cohort"] == held]
            if info.loc[te_idx, "y"].nunique() < 2:
                continue
            # normalize per cohort separately for train (each cohort) and test
            def build(idxs):
                parts = []
                for c in sel.loc[idxs, "cohort"].unique():
                    cidx = [i for i in idxs if sel.loc[i, "cohort"] == c]
                    parts.append(norm_cohort(counts.loc[genes, cidx]).T)
                return pd.concat(parts).fillna(0.0)
            Xtr = build(tr_idx); Xte = build(te_idx)
            Xtr = Xtr.reindex(columns=genes).fillna(0.0)
            Xte = Xte.reindex(columns=genes).fillna(0.0)
            ytr = info.loc[Xtr.index, "y"].values
            model = fit_l1(Xtr.values, ytr)
            prob = model.predict_proba(Xte.values)[:, 1]
            oof_prob += list(prob)
            oof_y += list(info.loc[Xte.index, "y"].values)
            oof_pat += list(sel.loc[Xte.index, "patient"].values)
        if len(set(oof_y)) < 2:
            continue
        auc, lo, hi = boot_auc(np.array(oof_y), np.array(oof_prob),
                               np.array(oof_pat))
        rows.append({"space": space_name, "timepoint": tp_name,
                     "n_folds": len(cohorts), "n_patients": len(oof_y),
                     "n_PE": int(sum(oof_y)), "AUC": round(auc, 3),
                     "CI_low": round(lo, 3), "CI_high": round(hi, 3)})
        print(f"  {space_name:12s} | {tp_name:20s} | n={len(oof_y):3d} "
              f"(PE={int(sum(oof_y)):2d}) | AUC={auc:.3f} [{lo:.3f}-{hi:.3f}]")

res = pd.DataFrame(rows)
res.to_csv("results/native_cfrna.csv", index=False)

print("\n" + "=" * 60)
print("CONTRAST: transfer (placenta->cfRNA) vs native cfRNA")
print("=" * 60)
for tp, ref in TRANSFER_REF.items():
    nat = res[(res.space == "transferable") & (res.timepoint == tp)]
    if len(nat):
        n = nat.iloc[0]
    print(f"  {tp:20s}: transfer={ref:.3f}  ->  native={n['AUC']:.3f} "
              f"[{n['CI_low']:.3f}-{n['CI_high']:.3f}]  (same 6184-feature space)")
print("\nSaved -> results/native_cfrna.csv")
