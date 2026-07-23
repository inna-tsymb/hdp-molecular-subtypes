#!/usr/bin/env python3
"""
s5_branch2_ml/02_select_panel.py  (v2 — stabilized selection)
--------------------------------------------------------------
Selects a predictive PE vs control panel (Branch 2), leak-resistant pipeline.

Changes from v1 (fixing unstable regularization):
  * IN-FOLD prefilter: top-K_PREFILTER genes by univariate f_classif on
    TRAINING data, then LASSO. Removes noise before regularization.
  * RESTRICTED C grid (prevent near-zero penalty) -> models are
    consistently sparse across folds.
  * Explicit LogisticRegression(penalty='l1') + manual inner CV for C.

Design (unchanged):
  * Feature space = ML ∩ Moufarrej (transferable).
  * LODO-CV, feature selection inside the fold -> honest evaluation.
  * Panels: stable (>= STABLE_K folds) and top25. External test is Moufarrej.

Output: data/panel_stable.txt, panel_top25.txt, transferable_genes.txt
       results/ml_lodo_auc.csv, ml_feature_stability.csv

Run from project root:  python scripts/s5_branch2_ml/02_select_panel.py
"""

import os
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore", category=FutureWarning)

K_PREFILTER = 500                       # top univariate genes before LASSO
C_GRID = np.logspace(-2.0, 0.3, 8)      # limited grid -> sparsity
STABLE_K = 3
TOP_N = 25
os.makedirs("results", exist_ok=True)

expr = pd.read_csv("data/ml_expr.csv", index_col=0)
meta = pd.read_csv("data/ml_metadata.csv")
assert list(meta["sample"]) == list(expr.index)
y = (meta["category"].values == "Case").astype(int)
datasets = meta["dataset"].values

mou = set(pd.read_csv("data/moufarrej/GSE192902_counts_Discovery_postQC.csv.gz",
                      usecols=["gene_name"])["gene_name"].astype(str))
transfer = [g for g in expr.columns if g in mou]
expr = expr[transfer].fillna(0.0)
with open("data/transferable_genes.txt", "w") as fh:
    fh.write("\n".join(transfer))
print(f"Transferable genes (ML ∩ Moufarrej): {len(transfer)}")
print(f"Matrix: {expr.shape[0]} samples x {expr.shape[1]} genes\n")

X = expr.values
gene_names = np.array(expr.columns)
N = len(gene_names)


def make_pipe(C):
    return Pipeline([
        ("sel", SelectKBest(f_classif, k=min(K_PREFILTER, N))),
        ("clf", LogisticRegression(penalty="l1", solver="saga", C=C,
                                   class_weight="balanced", max_iter=10000,
                                   random_state=0)),
    ])


def fit_l1(Xtr, ytr):
    """Pre-filter + L1, C tuned with inner CV on training data."""
    inner = StratifiedKFold(4, shuffle=True, random_state=0)
    best_C, best_s = C_GRID[0], -1.0
    for C in C_GRID:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = cross_val_score(make_pipe(C), Xtr, ytr, cv=inner,
                                scoring="roc_auc").mean()
        if s > best_s:
            best_s, best_C = s, C
    pipe = make_pipe(best_C)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pipe.fit(Xtr, ytr)
    return pipe


def selected_mask(pipe):
    support = pipe.named_steps["sel"].get_support()
    coef = pipe.named_steps["clf"].coef_[0]
    full = np.zeros(N)
    full[support] = coef
    return full


# ── LODO-CV ──────────────────────────────────────────────────────────
fold_auc = []
selection_count = np.zeros(N, dtype=int)
coef_accum = np.zeros(N)

for ds in sorted(set(datasets)):
    te, tr = datasets == ds, datasets != ds
    if len(set(y[te])) < 2:
        print(f"  {ds}: skipping (only one class in test set)")
        continue
    pipe = fit_l1(X[tr], y[tr])
    prob = pipe.predict_proba(X[te])[:, 1]
    auc = roc_auc_score(y[te], prob)
    coef = selected_mask(pipe)
    nsel = int((coef != 0).sum())
    fold_auc.append({"held_out": ds, "n_test": int(te.sum()),
                     "n_features": nsel, "auc": round(auc, 3)})
    selection_count += (coef != 0).astype(int)
    coef_accum += np.abs(coef)
    print(f"  {ds:11s}: AUC={auc:.3f}  features={nsel}")

auc_df = pd.DataFrame(fold_auc)
auc_df.to_csv("results/ml_lodo_auc.csv", index=False)
print(f"\nLODO-AUC mean: {auc_df['auc'].mean():.3f} (±{auc_df['auc'].std():.3f})")

# ── stability + panels ────────────────────────────────────────────
stab = (pd.DataFrame({"gene": gene_names, "n_folds_selected": selection_count,
                      "mean_abs_coef": np.round(coef_accum / max(len(fold_auc), 1), 4)})
        .sort_values(["n_folds_selected", "mean_abs_coef"], ascending=False))
stab.to_csv("results/ml_feature_stability.csv", index=False)
panel_stable = stab[stab["n_folds_selected"] >= STABLE_K]["gene"].tolist()

full = fit_l1(X, y)
coef_full = pd.Series(np.abs(selected_mask(full)), index=gene_names)
panel_top25 = coef_full[coef_full > 0].sort_values(ascending=False).head(TOP_N).index.tolist()

with open("data/panel_stable.txt", "w") as fh:
    fh.write("\n".join(panel_stable))
with open("data/panel_top25.txt", "w") as fh:
    fh.write("\n".join(panel_top25))

def read_deg(p):
    return set(pd.read_csv(p)["gene"]) if os.path.exists(p) else set()
deg_all = read_deg("results/deg_GSE203507_sig.csv") | read_deg("results/deg_GSE306864_male_sig.csv")

print("\n" + "=" * 55)
print(f"STABLE panel (>= {STABLE_K} folds): {len(panel_stable)} genes")
print(f"TOP-{TOP_N} panel: {len(panel_top25)} genes")
print(f"Overlap stable ∩ top25: {len(set(panel_stable) & set(panel_top25))}")
print(f"\nIn Branch 1 DEG (EOPE∪LOPE):")
print(f"  stable: {len(set(panel_stable) & deg_all)}/{len(panel_stable)}")
print(f"  top25 : {len(set(panel_top25) & deg_all)}/{len(panel_top25)}")
print("\nTOP-15 stable features:")
print(stab.head(15).to_string(index=False))
print("\nSaved: data/panel_*.txt, results/ml_lodo_auc.csv, ml_feature_stability.csv")
