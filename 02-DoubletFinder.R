# Some required R packages
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(glmGamPoi)

scObject <- readRDS("scObject.rds")

## Homotypic Doublet Proportion Estimate
## ~1000 cells, and a multiplet rate of ~0.8% according to 10X Genomics
Doublet <- sapply(scObject, function(seu)round(ncol(seu)*8*1e-6,digits = 4))

scObject <- lapply(names(scObject), FUN = function(x){
  print(x)
  seurat1 <- scObject[[x]]
  
  print("Clustering...")
  seurat1 <- seurat1 %>% 
    NormalizeData() %>% 
    FindVariableFeatures(nfeatures = 3000) %>% 
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE, seed.use = 111) %>% 
    RunUMAP(dims = 1:20, verbose = FALSE, seed.use = 111) %>% 
    FindNeighbors(dims = 1:20, verbose = FALSE) %>% 
    FindClusters(resolution = 0.5, verbose = FALSE, random.seed = 111)
  
  ## pK Identification
  print("pK Identification...")
  sweep.res.list <- paramSweep(seurat1, PCs = 1:20, sct = TRUE, num.cores = 1)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  saveRDS(list(sweep.res.list, sweep.stats, bcmvn, pK_bcmvn),
          file = paste0("doublet", x, "_sweep.rds"))
  
  ## Automatic doublet rate estimation
  print("Proportion Estimate...")
  DoubletRate <- Doublet[x]
  print(paste("Using DoubletRate =", DoubletRate))
  
  homotypic.prop <- modelHomotypic(seurat1$seurat_clusters)  
  nExp_poi <- round(DoubletRate * nrow(seurat1@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  print("doubletFinder_v3...")
  seurat1 <- doubletFinder(seu = seurat1,
                           PCs = 1:20,
                           pN = 0.25,
                           pK = pK_bcmvn, 
                           nExp = nExp_poi.adj,
                           reuse.pANN = FALSE,
                           sct = F)
  
  ## UMAP plot by doublet classification
  group1 <- names(seurat1@meta.data)[grep("^DF", names(seurat1@meta.data))]
  pdf(file = paste0("doublet", x, "_doublet.pdf"), width = 5, height = 5)
  print(DimPlot(seurat1, reduction = "umap", group.by = group1))
  dev.off()
  
  ## pANN score plot
  group2 <- names(seurat1@meta.data)[grep("^pANN", names(seurat1@meta.data))]
  pdf(file = paste0("doublet", x, "_pANN.pdf"), width = 5, height = 5)
  print(FeaturePlot(seurat1, reduction = "umap", features = group2))
  dev.off()
  
  return(seurat1)
})










