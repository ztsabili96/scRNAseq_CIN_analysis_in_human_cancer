scRNAseq_CIN_analysis_in_human_cancer
# scRNA-seq Chromosomal Instability Analysis in Human Cancer

This repository contains an R-based pipeline for inferring large-scale copy number variations
(CNVs) from single-cell RNA-seq data using copyKAT and quantifying chromosomal instability (CIN)
at the single-cell level.
The workflow is intentionally kept modular and reproducible, following common
practices in single-cell cancer genomics analyses of chromosomal instability.

## Project overview
- CNV inference using copyKAT
- Classification of cells as diploid or aneuploid
- Quantification of chromosomal instability (CIN score)

## Project Structure
```text
scRNAseq_CNV_Human_Cancer/
├── README.md
├── scripts/
│   ├── 01_load_data_and_qc.R
│   ├── 02_run_copyKAT.R
│   └── 03_calculate_CIN_score.R
├── data/
│   └── raw/
├── results/
│   ├── CIN_score_per_cell.csv
│   └── copykat/
│       ├── CID3586_copykat_clustering_results.rds
│       └── CID3586_copykat_prediction.txt
└── figures/
    ├── QC_scatter_plots.png
    ├── QC_violin_plots.png
    └── CIN_score_boxplot.png


## Analysis Workflow

The analysis was carried out in three main steps.

First, raw count data were loaded into Seurat and standard quality control was
performed. Cells were filtered based on the number of detected genes, total UMI
counts, and mitochondrial gene percentage. QC metrics were visualized to assess
filtering thresholds.

Next, copy number variation was inferred using copyKAT. Based on copyKAT
predictions, cells were classified as aneuploid (putative tumor cells) or
diploid (reference cells).

Finally, copyKAT CNV profiles were integrated with per-cell metadata to compute
a Chromosomal Instability (CIN) score for each cell. CIN score distributions were
visualized and compared between aneuploid and diploid populations.

## Results

The main outputs of this analysis include:

- `CIN_score_per_cell.csv`, containing per-cell CIN scores together with copyKAT
  classification labels
- copyKAT CNV prediction files stored in the `results/copykat/` directory
- Quality control (QC) plots summarizing filtering metrics prior to CNV inference
- A boxplot summarizing CIN score differences between aneuploid and diploid cells

### Quality Control of scRNA-seq Data

Cell- and gene-level quality control metrics were evaluated prior to CNV inference
to remove low-quality cells and potential technical artifacts.
![QC scatter plots](figures/QC_scatter_plots.png)

![QC violin plots](figures/QC_violin_plots.png)

### Chromosomal Instability (CIN) Score

Chromosomal instability was quantified at the single-cell level as the mean absolute
deviation of inferred copy number values from the diploid state.

As expected, aneuploid cells show higher CIN scores compared to diploid cells,
consistent with increased chromosomal instability in cancer.

<img src="figures/CIN_score_boxplot.png" width="650">

## Reproducibility

All analysis steps are documented in the `scripts/` directory and can be executed
sequentially. Raw input matrices are not included in the repository due to their
large size, but the required input formats and preprocessing steps are described
in the scripts. The final script reports R session information (`sessionInfo()`)
to facilitate software environment reproducibility.

## Data Availability

Due to GitHub file size limitations, large copyKAT CNA matrices 
(e.g. gene-by-cell CNV profiles) were not uploaded.

The raw scRNA-seq data are publicly available from GEO (GSE176078).

Processed outputs required to reproduce the main analyses 
(copyKAT predictions, CIN scores, and scripts) are provided in this repository.
