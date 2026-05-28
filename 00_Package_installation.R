## ================================================================
## Install required R packages
## ================================================================

## CRAN packages
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE)
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

if (!requireNamespace("NMF", quietly = TRUE)) {
  install.packages("NMF")
}

if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (!requireNamespace("survminer", quietly = TRUE)) {
  install.packages("survminer")
}

if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}


## Bioconductor packages
if (!requireNamespace("AUCell", quietly = TRUE)) {
  BiocManager::install("AUCell")
}

if (!requireNamespace("monocle", quietly = TRUE)) {
  BiocManager::install("monocle")
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  BiocManager::install("glmGamPoi")
}


## GitHub packages
if (!requireNamespace("Startrac", quietly = TRUE)) {
  remotes::install_github("Japrin/STARTRAC", quiet = TRUE)
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat", quiet = TRUE)
}

if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  remotes::install_github("Danko-Lab/BayesPrism/BayesPrism", quiet = TRUE)
}

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
}


## Save session information
sink("sessionInfo.txt")
sessionInfo()
sink()
