#!/usr/bin/env python3
"""
s4_branch1_dea/03_nmf_consensus.py  (uncontrolled subtype check)
-----------------------------------------------------------------------
Consensus NMF on PE samples: does unsupervised learning recover
an EOPE/LOPE split without labels?

IMPORTANT (confound): onset ≈ dataset (EOPE=203507, LOPE=306864), so
clusters may reflect BATCH rather than biology. Therefore:
  * compute ARI vs onset AND vs dataset (to distinguish them);
  * inspect where mixed samples (190971) fall — a biological test.

Method: consensus clustering (Monti) — 100 NMF subsample runs per k,
consensus matrix, cophenetic correlation for k selection, ARI + silhouette.

Input:  data/dea_counts.csv, data/dea_metadata.csv
Output: results/nmf_consensus.csv, results/nmf_consensus_heatmap_k*.pdf

Run from project root:  python scripts/s4_branch1_dea/03_nmf_consensus.py
"""

import os
import warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import adjusted_rand_score, silhouette_score
from scipy.cluster.hierarchy import linkage, cophenet, fcluster
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)
RNG = np.random.default_rng(0)
N_HVG = 2000
N_RUNS = 100
KS = [2, 3, 4]

dea = pd.read_csv("data/dea_counts.csv", index_col=0)
meta = pd.read_csv("data/dea_metadata.csv")
pe = meta[meta["category"] == "Case"].reset_index(drop=True)
print(f"PE samples: {len(pe)}  ({pe['onset'].value_counts().to_dict()})")

# ── non-negative input: log-CPM + top HVG ──────────────────────────────
X = dea[pe["sample"]].T.astype(float)                 # samples x genes (raw)
cpm = X.div(X.sum(axis=1), axis=0) * 1e6
logx = np.log1p(cpm)
hvg = logx.var().sort_values(ascending=False).head(N_HVG).index
M = logx[hvg].values                                  # >= 0, samples x HVG
n = M.shape[0]

onset = pe["onset"].values
dataset = pe["dataset"].values


def consensus_matrix(k):
    conn = np.zeros((n, n)); cnt = np.zeros((n, n))
    for _ in range(N_RUNS):
        idx = RNG.choice(n, int(0.8 * n), replace=False)
        W = NMF(n_components=k, init="nndsvda", max_iter=500,
                random_state=int(RNG.integers(1e6))).fit_transform(M[idx])
        lab = W.argmax(axis=1)
        for a in range(len(idx)):
            for b in range(len(idx)):
                conn[idx[a], idx[b]] += (lab[a] == lab[b])
                cnt[idx[a], idx[b]] += 1
    with np.errstate(invalid="ignore"):
        C = np.divide(conn, cnt, out=np.zeros_like(conn), where=cnt > 0)
    return C


rows = []
for k in KS:
    C = consensus_matrix(k)
    dist = 1 - C
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist, checks=False), method="average")
    coph = cophenet(Z, squareform(dist, checks=False))[0]
    clusters = fcluster(Z, k, criterion="maxclust")
    sil = silhouette_score(M, clusters) if len(set(clusters)) > 1 else np.nan
    ari_onset = adjusted_rand_score(onset, clusters)
    ari_dataset = adjusted_rand_score(dataset, clusters)
    rows.append({"k": k, "cophenetic": round(coph, 3),
                 "silhouette": round(sil, 3),
                 "ARI_onset": round(ari_onset, 3),
                 "ARI_dataset": round(ari_dataset, 3)})
    print(f"\nk={k}: cophenetic={coph:.3f} silhouette={sil:.3f} "
          f"ARI(onset)={ari_onset:.3f} ARI(dataset)={ari_dataset:.3f}")
    ct = pd.crosstab(pd.Series(onset, name="onset"),
                     pd.Series(clusters, name="cluster"))
    print(ct.to_string())

    # consensus heatmap, ordered by clusters
    order = np.argsort(clusters)
    plt.figure(figsize=(6, 5))
    plt.imshow(C[np.ix_(order, order)], cmap="RdBu_r", vmin=0, vmax=1)
    plt.colorbar(label="consensus"); plt.title(f"Consensus NMF k={k}")
    plt.xticks([]); plt.yticks([])
    plt.tight_layout(); plt.savefig(f"results/nmf_consensus_heatmap_k{k}.pdf"); plt.close()

res = pd.DataFrame(rows)
res.to_csv("results/nmf_consensus.csv", index=False)

# ── mixed test: where do GSE190971 samples fall? ─────────────────────────
k_best = res.loc[res["cophenetic"].idxmax(), "k"]
C = consensus_matrix(int(k_best))
dist = 1 - C; np.fill_diagonal(dist, 0)
Z = linkage(squareform(dist, checks=False), method="average")
clusters = fcluster(Z, int(k_best), criterion="maxclust")
mixed_mask = dataset == "GSE190971"
print("\n" + "=" * 55)
print(f"Best k by cophenetic: {k_best}")
print("Mixed samples (GSE190971) cluster distribution:",
      pd.Series(clusters[mixed_mask]).value_counts().to_dict())
print("\nInterpretation:")
print("  ARI(onset) ≈ ARI(dataset) -> clusters likely batch-driven, cannot separate")
print("  ARI(onset) > ARI(dataset) or mixed samples distribute across clusters -> biological signal")
print("\nSaved: results/nmf_consensus.csv, nmf_consensus_heatmap_k*.pdf")
