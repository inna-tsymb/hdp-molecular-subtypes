#!/usr/bin/env python3
"""
s5_branch2_ml/05_trees_shap.py
------------------------------
Nonlinear models (RandomForest + XGBoost) with LODO-CV and tree SHAP.
Covers planned items: RF/XGBoost, SHAP summary, SHAP force/waterfall
for example patients (hypoxic EOPE-like vs immune LOPE-like).

  1) LODO-CV: RF and XGB on transferable space -> AUC (compare with LASSO 0.86).
  2) SHAP summary (global) on the final RF (top25 panel).
  3) SHAP waterfall for 2 example Case patients:
       - driven by HK2 (hypoxic -> EOPE-like)
       - driven by CX3CR1 (immune -> LOPE-like)

Output: results/trees_lodo_auc.csv, results/shap_rf_summary.pdf,
       results/shap_waterfall_hypoxic.pdf, results/shap_waterfall_immune.pdf

Run from project root:  python scripts/s5_branch2_ml/05_trees_shap.py
Requires: shap, matplotlib, (optional) xgboost
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)

expr = pd.read_csv("data/ml_expr.csv", index_col=0)
meta = pd.read_csv("data/ml_metadata.csv")
y = (meta["category"].values == "Case").astype(int)
datasets = meta["dataset"].values

transfer = [g for g in open("data/transferable_genes.txt").read().split()
            if g in expr.columns]
panel = [g for g in open("data/panel_top25.txt").read().split()
         if g in expr.columns]
Xt = expr[transfer].fillna(0.0).values

# ── models ───────────────────────────────────────────────────────────
def rf_model():
    return RandomForestClassifier(n_estimators=500, class_weight="balanced",
                                  max_features="sqrt", random_state=0, n_jobs=-1)

def xgb_model(scale_pos_weight=1.0):
    import xgboost as xgb
    return xgb.XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.1,
                             subsample=0.8, colsample_bytree=0.8,
                             eval_metric="logloss", scale_pos_weight=scale_pos_weight,
                             random_state=0, n_jobs=-1)

HAVE_XGB = True
try:
    import xgboost  # noqa  (on macOS this may fail without libomp.dylib)
except Exception as e:
    HAVE_XGB = False
    print(f"[XGBoost unavailable -> only RandomForest. Reason: {type(e).__name__}]")
    print(" macOS: `brew install libomp` or `conda install -c conda-forge libomp`")

# ── LODO-CV ──────────────────────────────────────────────────────────
def lodo(make_model, is_xgb=False):
    aucs = []
    for ds in sorted(set(datasets)):
        te, tr = datasets == ds, datasets != ds
        if len(set(y[te])) < 2:
            continue
        if is_xgb:
            spw = (y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1)
            m = make_model(spw)
        else:
            m = make_model()
        m.fit(Xt[tr], y[tr])
        p = m.predict_proba(Xt[te])[:, 1]
        aucs.append({"held_out": ds, "auc": round(roc_auc_score(y[te], p), 3)})
    return pd.DataFrame(aucs)

rows = []
rf_auc = lodo(rf_model)
rf_auc["model"] = "RandomForest"
rows.append(rf_auc)
print("RandomForest LODO:")
print(rf_auc.to_string(index=False))
print(f"  mean AUC: {rf_auc['auc'].mean():.3f}\n")

if HAVE_XGB:
    xgb_auc = lodo(xgb_model, is_xgb=True)
    xgb_auc["model"] = "XGBoost"
    rows.append(xgb_auc)
    print("XGBoost LODO:")
    print(xgb_auc.to_string(index=False))
    print(f"  mean AUC: {xgb_auc['auc'].mean():.3f}\n")

pd.concat(rows).to_csv("results/trees_lodo_auc.csv", index=False)

# ── SHAP on final RF (top25 panel, all data) ────────────────────
import shap
Xp = expr[panel].fillna(0.0)
rf = rf_model()
rf.fit(Xp.values, y)

explainer = shap.TreeExplainer(rf)
sv = explainer.shap_values(Xp.values)
# binary RF: use class 1 (Case)
if isinstance(sv, list):
    sv1, base1 = sv[1], explainer.expected_value[1]
elif np.ndim(sv) == 3:
    sv1, base1 = sv[:, :, 1], np.ravel(explainer.expected_value)[1]
else:
    sv1, base1 = sv, np.ravel(explainer.expected_value)[0]

# global summary
shap.summary_plot(sv1, Xp, feature_names=panel, show=False, max_display=25)
plt.tight_layout(); plt.savefig("results/shap_rf_summary.pdf"); plt.close()

# ── waterfall: hypoxic (HK2) vs immune (CX3CR1) patient ──────────
sv_df = pd.DataFrame(sv1, columns=panel, index=Xp.index)
case_mask = y == 1

def pick_patient(gene):
    if gene not in sv_df.columns:
        return None
    s = sv_df.loc[case_mask, gene]
    return s.idxmax()                       # Case where the gene contributes most to PE

def waterfall(sample_id, tag, title):
    if sample_id is None:
        print(f"  [{tag}: gene missing from panel, skipping]")
        return
    i = list(Xp.index).index(sample_id)
    try:
        expl = shap.Explanation(values=sv1[i], base_values=base1,
                                data=Xp.iloc[i].values, feature_names=panel)
        shap.plots.waterfall(expl, max_display=15, show=False)
    except Exception:
        # fallback: barh of SHAP contributions for this patient
        row = sv_df.loc[sample_id].sort_values()
        plt.figure(figsize=(7, 6))
        plt.barh(row.index, row.values,
                 color=["#c0392b" if v > 0 else "#2980b9" for v in row.values])
        plt.axvline(0, color="k", lw=0.5); plt.xlabel("SHAP (contribution to P(PE))")
    plt.title(title)
    plt.tight_layout()
    out = f"results/shap_waterfall_{tag}.pdf"
    plt.savefig(out); plt.close()
    print(f"  {tag}: patient {sample_id} -> {out}")

print("Waterfall examples:")
waterfall(pick_patient("HK2"), "hypoxic",
      "Hypoxic (EOPE-like) profile: prediction driven by HK2")
waterfall(pick_patient("CX3CR1"), "immune",
      "Immune (LOPE-like) profile: prediction driven by CX3CR1")

print("\nSaved: results/trees_lodo_auc.csv, shap_rf_summary.pdf,")
print("           shap_waterfall_hypoxic.pdf, shap_waterfall_immune.pdf")
