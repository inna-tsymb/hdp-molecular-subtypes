# Molecular subtypes of hypertensive disorders of pregnancy

**Early- and late-onset preeclampsia as molecularly distinct conditions: a stratified transcriptomic analysis of the placenta and the limits of signature transferability to cfRNA**

MSc in Bioinformatics, Kyiv School of Economics · Author: Inna Kucherova · Supervisor: Ruslan Rodrigues · 2026

Preeclampsia (PE) is clinically split into early-onset (<34 weeks) and late-onset (≥34 weeks) forms, but whether these are molecularly distinct conditions or a severity continuum is unresolved. This project tests that question on public placental RNA-seq data, builds a compact predictive panel, and asks whether the panel transfers to cell-free RNA (cfRNA) in maternal plasma.

**Hypotheses**

- **H1** — EOPE and LOPE have distinct transcriptomic signatures. → **supported**
- **H2** — a compact panel separates PE from control in placental tissue across cohorts. → **supported** (AUC 0.861)
- **H3** — the placental panel retains predictive power in plasma cfRNA. → **not supported** (best AUC 0.657)

---

## Headline results

| # | Result |
| --- | --- |
| 1 | EOPE (283 DEGs) and LOPE (117 DEGs) share **4 genes and 0 enriched pathways**. |
| 2 | EOPE hypoxic signature (KEGG HIF-1 signalling, p-adj ≈ 8×10⁻⁴) **replicated independently** in GSE148241: 166/166 shared genes concordant in direction. |
| 3 | LOPE biology is membrane-repolarisation and immune, with **no hypoxic term**, and is detectable **only in male-bearing placentas** (117 vs 15 DEGs). |
| 4 | 25-gene panel: cross-cohort **AUC 0.861** (LASSO) > RF 0.786 > XGBoost 0.780. The LOPE cohort is the hardest fold in all three models. |
| 5 | Transfer to cfRNA is limited (**AUC 0.657**, 13–20 weeks); natively trained cfRNA models do no better (0.613); **0 of 25** panel genes are DE in plasma. |

**Interpretation:** onset-based subtypes are real and biologically separable in tissue, but the tissue→plasma compartment gap means cfRNA predictors should be built directly on plasma rather than transferred from tissue signatures.

---

## Design

Two branches on the same discovery pool, plus an external validation stage.

- **Branch 1 — stratified DEA.** Onset is almost fully confounded with cohort (EOPE ≈ GSE203507, LOPE ≈ GSE306864), so merging with batch correction would erase the onset signal together with the batch effect. Each dataset is analysed separately: DESeq2, `~ category` (GSE306864: `~ sex + category`), apeglm shrinkage, p-adj < 0.05 and |log₂FC| > 1. Enrichment: clusterProfiler (GO BP, KEGG), universe = all tested genes. Replication on GSE148241 with the identical protocol.
- **Branch 2 — predictive panel.** Per-dataset normalisation (log1p → z-scoring within dataset), which makes raw counts, RPKM and pre-normalised values comparable, acts as batch correction, and does not leak across folds. Feature space restricted to genes transferable to cfRNA (6,184). Leave-One-Dataset-Out CV with **feature selection inside each fold** (univariate pre-filter + LASSO). LASSO vs RandomForest vs XGBoost. Interpretation via SHAP (global summary + local waterfall).
- **cfRNA validation.** The panel is trained on the full placental cohort, **frozen**, and applied to cfRNA. Patient-level evaluation with bootstrap CIs, stratified by sampling timepoint. As a controlled contrast, models are also trained natively on cfRNA (LODO across three sub-cohorts) in the same feature space, to separate *"the signal does not transfer"* from *"there is no signal"*. Plasma DE by timepoint (DESeq2, `~ cohort + condition`, one sample per patient).
- **Unsupervised check.** Consensus NMF (100 subsampling runs) on 35 PE samples, with ARI computed separately against onset and against dataset to expose confounding.

---

## Cohorts (datasets actually used)

| GSE | Technology | n | Format | Role | Onset |
| --- | --- | --- | --- | --- | --- |
| GSE203507 | RNA-seq | 20\* | raw counts | DEA + ML | EOPE |
| GSE306864 | RNA-seq | 58 | raw counts | DEA + ML | LOPE |
| GSE190971 | RNA-seq | 13 | raw counts | DEA + ML | mixed |
| GSE186257 | RNA-seq | 44 | DESeq2-norm | ML | severe |
| GSE148241 | RNA-seq | 41 | raw counts | ML + replication | EOSPE |
| GSE192902 | RNA-seq | 199† | raw counts | cfRNA validation | — |

`*` villous-only after tissue filtering (from 57 samples). `†` unique patients after deduplication of longitudinal draws; three sub-cohorts used as folds in the native cfRNA analysis.

**Downloaded but not used in the final analyses:** GSE154377 (second cfRNA cohort) and GSE75010 (microarray) — acquisition and annotation code is retained under `scripts/s1_acquisition/` and `scripts/s2_annotation/`, but they are not integrated into any reported result.
**Excluded:** GSE97320 — peripheral blood microarray, myocardial infarction, n = 6; not PE-relevant.

All data are public from NCBI GEO. No restricted or request-only data are included.

---

## Results in detail

### Branch 1 — stratified differential expression

| Stratum | n (Case/Control) | DEGs |
| --- | --- | --- |
| EOPE (GSE203507) | 15 / 5 | 283 |
| LOPE, both sexes (GSE306864) | 13 / 45 | 15 |
| LOPE, male-bearing only | 5 / 16 | 117 |
| mixed (GSE190971) | 7 / 6 | 452 |
| EOPE ∩ LOPE / shared pathways | — | 4 (4/4 concordant) / 0 |

The contrast between the two LOPE rows shows how pooling fetal sexes masks the signal: three times fewer cases yield eight times more DEGs once female-bearing placentas are removed.

### Branch 2 — Leave-One-Dataset-Out cross-validation (ROC-AUC)

| Held-out dataset | LASSO | RF | XGBoost |
| --- | --- | --- | --- |
| GSE203507 (EOPE) | 1.000 | 0.947 | 1.000 |
| GSE190971 (mixed) | 1.000 | 0.881 | 0.952 |
| GSE186257 (severe) | 0.842 | 0.818 | 0.780 |
| GSE148241 (EOSPE) | 0.783 | 0.609 | 0.526 |
| GSE306864 (LOPE) | 0.679 | 0.675 | 0.643 |
| **Mean AUC** | **0.861** | 0.786 | 0.780 |

The linear model wins as expected under p ≫ n. AUC = 1.000 values come from small held-out sets (n = 13–20) and should be read with caution. Two panels were derived — top-25 and stable-18, sharing 14 genes. Leading SHAP features are **HK2** (canonical HIF-1 target, hypoxic axis) and **CX3CR1** (fractalkine receptor, immune axis), so the panel spans both subtype axes. Local SHAP shows two patients with near-identical model output (0.858 vs 0.854) but different decision structure — one HK2-dominated, the other a distributed non-hypoxic profile.

### cfRNA validation (ROC-AUC [95% CI])

| Timepoint | n (PE) | Transfer | Native (6,184) | Native (7,109) |
| --- | --- | --- | --- | --- |
| ≤12 weeks | 127 (32) | 0.545 | — | — |
| 13–20 weeks | 113 (34) | 0.657 [0.53–0.76] | 0.613 [0.52–0.70] | 0.598 [0.50–0.69] |
| ≥23 weeks | 82 (26) | 0.544 | 0.648 [0.52–0.77] | 0.707 [0.59–0.81] |

Plasma DE confirms the gap at gene level: 93 cfRNA DEGs in total (0 at ≤12 weeks, 76 at 13–20, 26 at ≥23), of which only 2/283 EOPE DEGs, 3/117 LOPE DEGs and **0/25 panel genes** overlap. Pooling timepoints collapses the signal (4 DEGs across 196 samples) — the plasma signature is timepoint-specific. Enrichment of the 13–20-week set returned no significant terms; dominance of XIST and myeloid markers (ITGAX, MPEG1, CD300E) points to cell-composition and fetal-fraction effects.

### Unsupervised check — negative

Consensus NMF gives highly stable clusters (cophenetic > 0.94) aligning perfectly with onset (ARI up to 1.0) — but ARI against onset equals ARI against dataset at every k, and mixed samples form a cohort-specific cluster. Because of this confounding, **clustering cannot serve as independent confirmation of the subtypes**, and is reported as such rather than as supporting evidence.

---

## Limitations

1. **Onset is confounded with cohort.** This is the principal design limitation: it precludes matrix merging and blocks independent verification of the subtypes by clustering.
2. **Small samples.** EOPE results are robust; LOPE (male-only, n = 5 cases) and the mixed cohort are preliminary. AUC = 1.000 on n = 13–20 held-out sets is not a performance claim.
3. **Temporal gap** between tissue (collected at delivery) and cfRNA (early plasma) constrains translational conclusions.
4. **Heterogeneous normalisation** across the ML pool, partially mitigated by per-dataset z-scoring.
5. **No cell-composition deconvolution or fetal-fraction correction** in the cfRNA analysis — a direct specification for future work.
6. **cfRNA validation rests on a single cohort** (three sub-cohorts).

This project does not propose a diagnostic tool. Its contribution is an independent, reproducible confirmation of the subtyping hypothesis on open data, and an empirical argument for building non-invasive predictors directly on plasma.

---

## Repository structure

```
hdp-molecular-subtypes/
├── README.md
├── requirements.txt
├── install_R_deps.R           # R dependencies for Branch 1
├── data/
│   ├── metadata/              # GEO SOFT + phenotype CSV
│   ├── raw_counts/            # per-GSM counts (not tracked)
│   ├── series_counts/         # series matrices
│   ├── moufarrej/             # cfRNA count matrices (GSE192902)
│   ├── master_metadata.csv    # single source of truth for sample metadata
│   ├── cfrna_labels.csv
│   ├── dea_counts.csv         # merged DEA count matrix
│   ├── dea_metadata.csv
│   └── master_gene_list.txt   # common gene universe
├── results/                   # DEG lists, enrichment tables, panels, figures
└── scripts/
    ├── s1_acquisition/        # download and reformatting
    ├── s2_annotation/         # inventory, master metadata, cfRNA labels
    ├── s3_scouting/           # external cfRNA cohort scouting
    ├── s4_branch1_dea/        # merge → DESeq2 → enrichment
    ├── s5_branch2_ml/         # LODO ML pipeline + SHAP
    └── s6_validation/         # cfRNA transfer, native models, plasma DE
```

All scripts run **from the repository root**; paths like `data/...` are relative to it:

```bash
python scripts/s2_annotation/02_build_master_metadata.py
```

`data/raw_counts/` holds raw GEO downloads and is not tracked. Reproduce it with `scripts/s1_acquisition/`.

---

## Run order

```bash
# Setup
pip install -r requirements.txt
Rscript install_R_deps.R

# s1 · Data acquisition
python scripts/s1_acquisition/01_download_geo.py
python scripts/s1_acquisition/02_flatten_per_gsm.py
python scripts/s1_acquisition/03_flatten_final.py
bash   scripts/s1_acquisition/download_series.sh

# s2 · QC and annotation
python scripts/s2_annotation/01_inventory_datasets.py
python scripts/s2_annotation/02_build_master_metadata.py
python scripts/s2_annotation/03_build_cfrna_labels.py

# s3 · Cohort scouting (exploratory)
python scripts/s3_scouting/01_scout_cfrna_cohorts.py

# s4 · Branch 1 — stratified DEA
python  scripts/s4_branch1_dea/01_merge_dea_matrices.py
Rscript scripts/s4_branch1_dea/02_dea_stratified.R
Rscript scripts/s4_branch1_dea/03_enrichment.R

# s5 · Branch 2 — ML panel
# <fill in actual filenames>

# s6 · cfRNA validation
# <fill in actual filenames>
```

## Key output files

| Path | Contents |
| --- | --- |
| `results/...` | DEG lists per stratum (EOPE, LOPE, LOPE-male, mixed, GSE148241) |
| `results/...` | GO BP and KEGG enrichment tables |
| `results/...` | top-25 and stable-18 panel gene lists |
| `results/...` | LODO AUC tables, SHAP plots, cfRNA ROC curves |

<!-- replace with real paths -->

---

## Reproducibility

- **Python:** GEOparse, pandas, scikit-learn, shap, mygene, xgboost — see `requirements.txt`.
- **R:** DESeq2, apeglm, clusterProfiler, org.Hs.eg.db — see `install_R_deps.R`, which writes `R_sessionInfo.txt`.
- All accession IDs and filtering rules are documented in `scripts/s2_annotation/` with explicit column mappings.

## Data availability and licensing

Analysis code in this repository is released under the **MIT License**. The underlying expression data are **not** covered by it: all datasets are public NCBI GEO records belonging to their original authors, and any reuse should cite the primary publications — in particular Leavey et al. (GSE75010), Aisagbonhi et al. 2025 (GSE306864) and Moufarrej et al., *Nature* 2022 (GSE192902).

## Citation

If you use this code, please cite the repository (see `CITATION.cff`) and the original data sources.

## Key references

- Staff AC. The two-stage placental model of preeclampsia: an update. *J Reprod Immunol.* 2019;134–135:1–10.
- Elovitz MA, et al. Molecular subtyping of hypertensive disorders of pregnancy. *Nat Commun.* 2025;16:2948.
- Aisagbonhi O, et al. Sex-specific placental transcriptome alterations in late-onset preeclampsia. *Biol Sex Differ.* 2025.
- Moufarrej MN, et al. Early prediction of preeclampsia in pregnancy with cell-free RNA. *Nature.* 2022;602:689–694.
- Love MI, Huber W, Anders S. DESeq2. *Genome Biol.* 2014;15:550.

Full reference list is in the project report.
