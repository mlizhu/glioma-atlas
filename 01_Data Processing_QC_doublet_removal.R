library(Seurat)
library(tidyverse)
library(patchwork)
library(DoubletFinder)
library(glmGamPoi)


###############################################################################

##                                 Read data

###############################################################################

setwd("D:/project/Glima")
list_files <- list.files()
n <- length(list_files)
dir <- paste("D:/project/Glima", list_files, sep = "/")
sample_name  <- list_files

count  <- NULL
seurat <- NULL
for (i in 1:n){
  count[[i]]  <- Read10X(dir[i]) 
  seurat[[i]] <- CreateSeuratObject(counts = count[[i]],
                                    project = sample_name[i],
                                    min.cells = 3,
                                    min.features = 200)
}

names(seurat) <- sample_name
names(count) <-  sample_name

###############################################################################

##                                 Preprocessing

###############################################################################


seurat <- lapply(seurat, function(seurat) {
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & 
                     nCount_RNA > 500 & nCount_RNA < 35000 & 
                     percent.mt < 20)
  return(seurat)
})

seurat.list <- seurat

Doublet <- sapply(seurat.list, function(seu)round(ncol(seu)*8*1e-6,digits = 4))

seurat.list <- lapply(names(seurat.list), FUN = function(x){
  print(x)
  seurat1 <- seurat.list[[x]]
  
  print("Clustering...")
  seurat1 <- seurat1 %>% 
    NormalizeData() %>% 
    FindVariableFeatures(nfeatures = 2000) %>% 
    ScaleData(verbose = TRUE) %>%
    RunPCA(verbose = FALSE, seed.use = 111) %>% 
    RunUMAP(dims = 1:20, verbose = FALSE, seed.use = 111) %>% 
    FindNeighbors(dims = 1:20, verbose = FALSE) %>% 
    FindClusters(resolution = 0.5, verbose = FALSE, random.seed = 111)
  
  ## pK Identification
  print("pK Identification...")
  sweep.res.list <- paramSweep(seurat1, PCs = 1:20, sct =  FALSE, num.cores = 1)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  saveRDS(list(sweep.res.list, sweep.stats, bcmvn, pK_bcmvn),
          file = paste0("doublet", x, "_sweep.rds"))
  
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


seurat.list <- lapply(seurat.list, FUN = function(x){
  print(table(x@meta.data[, 8]))
  x$"doublet" <- x@meta.data[, 8]
  x <- subset(x, doublet == "Singlet")
  x
})

sample <-names(seurat)
names(seurat.list) <- sample




