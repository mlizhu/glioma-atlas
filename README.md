##--------------Glioma single-cell analysis code----------------##
This repository contains the main analysis scripts used in the manuscript

The scripts include single-cell preprocessing, quality control, doublet removal, Harmony integration, gene-set scoring, CellChat analysis, BayesPrism deconvolution, Monocle2 pseudotime analysis, and spatial transcriptomic gene-signature scoring.

##---------------Code files------------------------------------##
00_Package_installation.R
01_Data Processing_QC_doublet_removal.R
02_Integration_Harmony.R
03_Fig3E_Fig4G_Gene_set_mapping.R
04_Fig3F_Fig3G_Fig4H_Fig4I_Gene_set_score.R
05_Fig3C_Spatial_score_WT179_MU68.R
06_Fig4F_Spatial score_MC83_MN30.R
07_Fig5_Cellchat.R
08_Fig6_Pericyte.R
09_Fig7_Myeloid.R
10_Fig8_T_cell.R


##-----------------Required input files----------------------##

The main input files used in the scripts include:

Glioma.rds
Pericyte.rds
Myeloid.rds
T_cell.rds
signature.xlsx
TCGAGBMLGG全数据.xlsx

`Glioma.rds` is the annotated major-cell-type Seurat object.
`Pericyte.rds` is refined pericyte Seurat object used for pericyte-focused subclustering
`Myeloid.rds` is refined myeloid Seurat object used for myeloid subclustering 
`T_cell.rds` is refined T-cell Seurat object used for T-cell subclustering
`signature.xlsx` contains the 179-gene, 68-gene, 30-gene, and 83-gene signatures.
`TCGAGBMLGG全数据.xlsx` contains the processed TCGA GBM_LGG bulk RNA-seq matrix with survival information used for BayesPrism and survival analyses.

##----------------Main analyses------------------------------##
The repository contains scripts for the following analyses:
1. Single-cell quality control, doublet removal, normalization, clustering, annotation, and Harmony integration.
2. Gene-set mapping of the 179-gene, 68-gene, 30-gene, and 83-gene signatures.
3. AUCell-based gene-set scoring and UMAP visualization.
4. CellChat-based cell-cell communication analysis.
5. Pericyte subcluster analysis, including cell proportion, Ro/e enrichment, BayesPrism deconvolution, survival analysis, and Monocle2 pseudotime analysis.
6. Myeloid cell downstream analysis, including subcluster analysis, Ro/e enrichment, BayesPrism deconvolution, survival analysis, AUCell scoring, and Monocle2 pseudotime analysis.
7. T/NK/B lineage cell downstream analysis, including subcluster analysis, Ro/e enrichment, BayesPrism deconvolution, survival analysis, AUCell scoring, and Monocle2 pseudotime analysis.
8. Spatial transcriptomic scoring of the 179-gene, 68-gene, 30-gene, and 83-gene signatures.
Some downstream scripts can be run independently if the required processed Seurat objects are available.

##-----------------Notes--------------------------------------##
The detailed R package versions are provided in `sessionInfo.txt`.




