# Glioma single-cell analysis code

This repository contains the main analysis scripts used in the manuscript.

The scripts include single-cell preprocessing, quality control, doublet removal, Harmony integration, gene-set scoring, CellChat analysis, BayesPrism deconvolution, Monocle2 pseudotime analysis

## Code files

* `00_Package_installation.R`
* `01_Data Processing_QC_doublet_removal.R`
* `02_Integration_Harmony.R`
* `03_Gene_set_mapping.R`
* `04_Gene_set_score.R`
* `05_ellchat.R`
* `06_Pericyte.R`

## Required input files

The main input files used in the scripts include:

* `Glioma.rds`
* `Pericyte.rds`
* `signature.xlsx`
* `TCGAGBMLGG全数据.xlsx`

`Glioma.rds` is the annotated major-cell-type Seurat object.

`Pericyte.rds` is refined pericyte Seurat object used for pericyte-focused subclustering.

`signature.xlsx` contains the 179-gene, 68-gene, 30-gene, and 83-gene signatures.

`TCGAGBMLGG全数据.xlsx` contains the processed TCGA GBM_LGG bulk RNA-seq matrix with survival information used for BayesPrism and survival analyses.

## Main analyses

The repository contains scripts for the following analyses:

1. Single-cell quality control, doublet removal, normalization, clustering, annotation, and Harmony integration.

2. Gene-set mapping of the 179-gene, 68-gene, 30-gene, and 83-gene signatures.

3. AUCell-based gene-set scoring and UMAP visualization.

4. CellChat-based cell-cell communication analysis.

5. Pericyte subcluster analysis, including cell proportion, Ro/e enrichment, BayesPrism deconvolution, survival analysis, and Monocle2 pseudotime analysis.


Some downstream scripts can be run independently if the required processed Seurat objects are available.

## Notes

The detailed R package versions are provided in `sessionInfo.txt`.
