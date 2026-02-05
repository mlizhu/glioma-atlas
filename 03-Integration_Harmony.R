# Some required R packages
library(Seurat)

# Step1 : Load seurat object 
seurat.combined <- readRDS("Glioma.Seurat.rds")

# Step2 : Perform standard preprocessing (log-normalization, identify variable features, scale and PCA)
seurat.combined <- seurat.combined %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nFeatures = 3000) %>% 
  ScaleData(verbose = TRUE) %>% 
  RunPCA(npcs = 20, verbose = FALSE, seed.use = 111) %>% 
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111)

## Step3 : Dimensional reduction and clustering
seurat.combined <- IntegrateLayers(
  object = seurat.combined,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

seurat.combined <- JoinLayers(seurat.combined)

seurat.combined <- seurat.combined %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(resolution = 0.5) %>% 
  RunUMAP(reduction = "harmony", reduction.name = "umap.harmony", dims = 1:20, seed.use = 111)

# Step4 : Finding differentially expressed genes
DefaultAssay(object = seurat.combined ) <- "RNA"

Idents(seurat.combined) <- seurat.combined$"RNA_snn_res.0.5"

markers  <- FindAllMarkers(seurat.combined, only.pos = TRUE, min.pct = 0.25, 
                           logfc.threshold = 0.25)





