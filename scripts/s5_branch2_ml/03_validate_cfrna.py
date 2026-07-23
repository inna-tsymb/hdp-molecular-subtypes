#!/usr/bin/env python3
"""
s5_branch2_ml/03_validate_cfrna.py
----------------------------------
Branch 2 finale: frozen placental panel -> external test on cfRNA
(Moufarrej, GSE192902). Answers the project's central question:
Does the placental PE signature transfer to NON-INVASIVE plasma?

Leak controls:
  * The model is trained on ALL placenta samples (panel genes), frozen,
    and only then applied to cfRNA. cfRNA is not touched during panel building.
  * Per-dataset normalization (log+z) — Moufarrej is standardized on its own.
  * Evaluation at the PATIENT level (one sample per patient within a timepoint)
    to avoid longitudinal pseudo-replication. AUC with bootstrap CI (resample patients).
  * Stratification by timepoint: ≥23 weeks (co-temporal) vs early plasma (predictive).

Panels: top25 (primary) + stable (core, appendix).
SHAP shows gene contributions to the final top25 model.

Output: results/cfrna_validation.csv, results/shap_top25_*.pdf

Run from project root:  python scripts/s5_branch2_ml/03_validate_cfrna.py
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)
RNG = np.random.default_rng(0)
B = 2000                       # bootstrap iterations

# ── placental training data ─────────────────────────────────────
expr = pd.read_csv("data/ml_expr.csv", index_col=0)
meta = pd.read_csv("data/ml_metadata.csv")
y_train = (meta["category"].values == "Case").astype(int)

panels = {
    "top25": open("data/panel_top25.txt").read().split(),
    "stable": open("data/panel_stable.txt").read().split(),
}


def train_frozen(panel):
    Xtr = expr[panel].values
    m = LogisticRegression(penalty="l2", class_weight="balanced",
                           max_iter=5000, random_state=0)
    m.fit(Xtr, y_train)
    return m


# ── cfRNA Moufarrej: matrix + labels ─────────────────────────────────
def load_moufarrej(panel):
    mats = []
    for f in sorted(glob.glob("data/moufarrej/GSE192902_counts_*postQC.csv.gz")):
        df = pd.read_csv(f, index_col=0)
        df = df.drop(columns=[c for c in ["gene_num"] if c in df.columns])
        mats.append(df)
    counts = pd.concat(mats, axis=1)
    counts = counts[~counts.index.duplicated(keep="first")]
    # panel genes (all transferable to Moufarrej by design)
    counts = counts.reindex(panel).dropna(how="all")
    # log + z-score per gene (Moufarrej is one dataset)
    x = np.log1p(counts.astype(float))
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1).replace(0, np.nan), axis=0)
    z = z.dropna(how="all").fillna(0.0)
    return z.T                                     # samples x genes


def attach_labels(Xcf):
    """Attach label metadata to cfRNA samples."""
    lab = pd.read_csv("data/cfrna_labels.csv")
    lab["key"] = lab["sample"].astype(str)
    Xcf = Xcf.copy()
    Xcf.index = [str(s).split(".")[0] for s in Xcf.index]   # drop .N duplicates
    Xcf = Xcf[~Xcf.index.duplicated(keep="first")]
    keep = [s for s in Xcf.index if s in set(lab["key"])]
    Xcf = Xcf.loc[keep]
    info = lab.set_index("key").loc[keep, ["patient", "label", "timepoint"]]
    return Xcf, info


def boot_auc(yv, sv, groups):
    """AUC + 95% CI, resample at the patient level (groups)."""
    base = roc_auc_score(yv, sv)
    uniq = np.array(sorted(set(groups)))
    g2idx = {g: np.where(groups == g)[0] for g in uniq}
    aucs = []
    for _ in range(B):
        samp = RNG.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([g2idx[g] for g in samp])
        yy, ss = yv[idx], sv[idx]
        if len(set(yy)) == 2:
            aucs.append(roc_auc_score(yy, ss))
    lo, hi = np.percentile(aucs, [2.5, 97.5]) if aucs else (np.nan, np.nan)
    return base, lo, hi


def one_per_patient(info, scores):
    """One sample per patient within each timepoint (first sample)."""
    d = info.copy()
    d["score"] = scores
    d["y"] = (d["label"] == "PE").astype(int)
    d = d.reset_index().drop_duplicates(subset=["patient"], keep="first")
    return d


# ── validation ────────────────────────────────────────────────────────
rows = []
TIMEPOINTS = ["≥23 weeks gestation", "13-20 weeks gestation", "≤12 weeks gestation"]

for pname, panel in panels.items():
    model = train_frozen(panel)
    Xcf = load_moufarrej(panel)
    Xcf, info = attach_labels(Xcf)
    # align gene order to the model
    Xcf = Xcf.reindex(columns=panel).fillna(0.0)
    prob = model.predict_proba(Xcf.values)[:, 1]

    for tp in TIMEPOINTS + ["ALL antenatal"]:
        if tp == "ALL antenatal":
            mask = info["timepoint"].isin(TIMEPOINTS).values
        else:
            mask = (info["timepoint"] == tp).values
        if mask.sum() < 6:
            continue
        d = one_per_patient(info[mask], prob[mask])
        if d["y"].nunique() < 2:
            continue
        auc, lo, hi = boot_auc(d["y"].values, d["score"].values, d["patient"].values)
        rows.append({"panel": pname, "timepoint": tp,
                     "n_patients": len(d), "n_PE": int(d["y"].sum()),
                     "AUC": round(auc, 3), "CI_low": round(lo, 3), "CI_high": round(hi, 3)})
        print(f"  {pname:7s} | {tp:22s} | n={len(d):3d} (PE={int(d['y'].sum()):2d}) "
              f"| AUC={auc:.3f} [{lo:.3f}-{hi:.3f}]")

res = pd.DataFrame(rows)
res.to_csv("results/cfrna_validation.csv", index=False)
print("\nSaved -> results/cfrna_validation.csv")

# ── SHAP for top25 ───────────────────────────────────────────────────
try:
    import shap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panel = panels["top25"]
    model = train_frozen(panel)
    Xtr = expr[panel]
    expl = shap.LinearExplainer(model, Xtr)
    sv = expl.shap_values(Xtr)
    shap.summary_plot(sv, Xtr, show=False, max_display=25)
    plt.tight_layout(); plt.savefig("results/shap_top25_summary.pdf"); plt.close()
    shap.summary_plot(sv, Xtr, plot_type="bar", show=False, max_display=25)
    plt.tight_layout(); plt.savefig("results/shap_top25_bar.pdf"); plt.close()
    print("SHAP -> results/shap_top25_summary.pdf, shap_top25_bar.pdf")
except Exception as e:
    print(f"[SHAP skipped: {e}]  (pip install shap matplotlib)")
