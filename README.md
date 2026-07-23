# Molecular subtypes of hypertensive pregnancy disorders

Course project (MSc Bioinformatics, KSE). Identifying molecular subtypes of preeclampsia (PE) in placenta using differential expression analysis and machine learning, with in silico validation of a predictive panel on non-invasive plasma cfRNA.

---

## Two-branch design

Because of small sample sizes and heterogeneous data formats, the discovery analysis is split into two branches:

* **Branch 1 (DEA)** — rigorous differential expression analysis on raw count data (DESeq2), stratified by onset. Provides biological interpretation with pathway enrichment (KEGG/GO).
* **Branch 2 (ML)** — predictive biomarker panel selection (LASSO/RF + SHAP) using a broader pool of datasets, including normalized data.
* **Validation** — the selected panel is tested on cfRNA plasma (Moufarrej GSE192902) as an external non-invasive validation.

---

## Cohorts (final audited set)

| GSE        | Modality   | n   | Format         | Branch       | Onset |
|------------|------------|-----|----------------|--------------|-------|
| GSE203507  | RNA-seq    | 20* | raw counts     | DEA + ML     | EOPE  |
| GSE306864  | RNA-seq    | 58  | raw counts     | DEA + ML     | LOPE  |
| GSE190971  | RNA-seq    | 13  | raw counts     | DEA + ML     | mixed |
| GSE186257  | RNA-seq    | 44  | DESeq2 norm    | ML           | severe|
| GSE148241  | RNA-seq    | 43  | RPKM           | ML           | EOPE  |
| GSE75010   | microarray | 157 | log2 norm      | ML (limma)   | mixed |
| GSE154377  | RNA-seq    | 17† | raw counts     | cfRNA valid. | —     |
| GSE192902  | RNA-seq    | 199†| raw counts     | cfRNA valid. | —     |

* `*` villous-only placenta after tissue filtering (from 57 samples).
* `†` unique patients after timepoint deduplication (not raw samples).

GSE97320 is excluded because it is only n=6 microarray and contains myocardial infarction samples, not PE.

All data are public from NCBI GEO. No closed or request-only data are included in the repository.

---

## Repository structure

```
course_project/
├── README.md
├── requirements.txt
├── install_R_deps.R          # R dependencies for Branch 1
├── data/
│   ├── metadata/             # GEO SOFT + phenotype CSV files
│   ├── raw_counts/           # per-GSM counts (not tracked in GitHub)
│   ├── series_counts/        # series matrices for 306864/190971/186257
│   ├── moufarrej/            # cfRNA count matrices for GSE192902
│   ├── master_metadata.csv   # single source of truth for sample metadata
│   ├── cfrna_labels.csv      # labels for cfRNA validation
│   ├── dea_counts.csv        # merged DEA count matrix 18325 × 91
│   ├── dea_metadata.csv      # metadata aligned with the DEA matrix
│   ├── master_gene_list.txt  # common gene universe
│   └── ...
├── results/                  # DEG lists, enrichment tables, plots
└── scripts/
    ├── s1_acquisition/       # data download and reformatting
    ├── s2_annotation/        # inventory, master metadata, cfRNA labels
    ├── s3_scouting/          # external cfRNA cohort scouting
    ├── s4_branch1_dea/       # merge -> DESeq2 -> enrichment
    ├── s5_branch2_ml/        # ML pipeline and validation
    └── s6_validation/        # cfRNA validation analyses
```

**Important:** All scripts run from the project root. Paths like `data/...` are relative to the repository root. Example:

```bash
python scripts/s2_annotation/02_build_master_metadata.py
```

**Note:** The `data/raw_counts/` folder contains raw GEO downloads and is not tracked in the repository. To reproduce the data, use the scripts in `scripts/s1_acquisition/` and `bash download_series.sh`.

---

## Run order

### Setup

```bash
pip install -r requirements.txt
Rscript install_R_deps.R
```

### s1 · Data acquisition

```bash
python scripts/s1_acquisition/01_download_geo.py       # GEO SOFT + per-GSM supplements
python scripts/s1_acquisition/02_flatten_per_gsm.py    # normalize file structure
python scripts/s1_acquisition/03_flatten_final.py
bash   scripts/s1_acquisition/download_series.sh       # series and cfRNA matrices
```

### s2 · QC and annotation

```bash
python scripts/s2_annotation/01_inventory_datasets.py   # audit platforms / formats / n
python scripts/s2_annotation/02_build_master_metadata.py # build master_metadata.csv
python scripts/s2_annotation/03_build_cfrna_labels.py   # build cfrna_labels.csv
```

### s3 · Cohort scouting (exploratory)

```bash
python scripts/s3_scouting/01_scout_cfrna_cohorts.py
```

### s4 · Branch 1 (DEA)

```bash
python scripts/s4_branch1_dea/01_merge_dea_matrices.py
Rscript scripts/s4_branch1_dea/02_dea_stratified.R
Rscript scripts/s4_branch1_dea/03_enrichment.R
```

---

## Key result of Branch 1

Stratified DEA indicates that **EOPE and LOPE are molecularly distinct states**: they share only a handful of overlapping DEGs and no shared enriched pathway at the pathway level.

* **EOPE** is driven by a hypoxic signature (HIF-1 signaling, response to hypoxia).
* **LOPE** is driven by membrane and vascular remodeling biology and is present primarily in male placentas.

Caution: the LOPE signature is based on a small male-only subset (n=5 cases), so it should be interpreted as a pathway-level hypothesis rather than a final gene catalog.

---

## Known limitations

1. **Onset is confounded with dataset** (EOPE ≈ GSE203507, LOPE ≈ GSE306864) — therefore the pipeline uses stratified DEA rather than merging with batch correction.
2. **Small discovery sample sizes** — EOPE results are stronger, LOPE results are more preliminary.
3. **Time shift from placenta (at delivery) to cfRNA plasma (early sampling)** — this is expected to reduce validation performance and is explicitly acknowledged.

---

## Reproducibility

* Python: use `pip freeze > requirements.lock.txt` after installing.
* R: `install_R_deps.R` writes `R_sessionInfo.txt`; for full reproducibility consider `renv::init()`.
* All accession IDs and filtering rules are documented in `scripts/s2_annotation/` with explicit column mappings.
