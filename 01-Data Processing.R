# required R packages
library(Seurat)

seurat <- readRDS("seurat.rds")

seurat <- lapply(seurat, function(seurat) {
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & 
                     nCount_RNA > 500 & nCount_RNA < 35000 & 
                     percent.mt < 20)
  return(seurat)
})






