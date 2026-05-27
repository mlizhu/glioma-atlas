Figure4F:
library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library('RColorBrewer')
######
col1<-c("grey90", "#D4B4DE", "#800080", "#4B0082")
pal1<-colorRampPalette(col1)
col2<-c("grey90", "#D2B48C", "#8B4513", "#5C4033")
pal2<-colorRampPalette(col2)
IDHmt_cod_83 <- read.table("/home/Cancer_Cell_10xSpatial/genelist/IDHmt_cod_83.txt",header=T)
IDHmt_cod_83 <- as.list(IDHmt_cod_83)
IDHmt_noncod_30 <- read.table("/home/Cancer_Cell_10xSpatial/genelist/IDHmt_noncod_30.txt",header=T)
IDHmt_noncod_30 <- as.list(IDHmt_noncod_30)

######UKF270T-IM-N######
UKF270T-IM-N <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF270_IDHMutant_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF270T-IM-N) <- "Spatial"
UKF270T-IM-N <- NormalizeData(UKF270T-IM-N)
expr_matrix <- GetAssayData(
  UKF270T-IM-N, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHmt_cod_83
cells_AUC <- AUCell_calcAUC(
  IDHmt_cod_83, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF270T-IM-N@meta.data$IDHmt_cod_83 <- auc_matrix[, "cod83"]
P1<-SpatialFeaturePlot(UKF270T-IM-N, features = "IDHmt_cod_83",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_noncod_30
cells_AUC <- AUCell_calcAUC(
  IDHmt_noncod_30, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF270T-IM-N@meta.data$IDHmt_noncod_30 <- auc_matrix[, "noncod"]
P2<-SpatialFeaturePlot(UKF270T-IM-N, features = "IDHmt_noncod_30",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF270T-IM-N_IDHmt_cod_83.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF270T-IM-N_IDHmt_noncod_30.png", plot = P2, width = 6, height = 8, dpi = 300)





######UKF268T-IM-C######
UKF268T-IM-C <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF268_IDHMutant_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF268T-IM-C) <- "Spatial"
UKF268T-IM-C <- NormalizeData(UKF268T-IM-C)
expr_matrix <- GetAssayData(
  UKF268T-IM-C, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHmt_cod_83
cells_AUC <- AUCell_calcAUC(
  IDHmt_cod_83, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF268T-IM-C@meta.data$IDHmt_cod_83 <- auc_matrix[, "cod83"]
P1<-SpatialFeaturePlot(UKF268T-IM-C, features = "IDHmt_cod_83",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_noncod_30
cells_AUC <- AUCell_calcAUC(
  IDHmt_noncod_30, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF268T-IM-C@meta.data$IDHmt_noncod_30 <- auc_matrix[, "noncod"]
P2<-SpatialFeaturePlot(UKF268T-IM-C, features = "IDHmt_noncod_30",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF268T-IM-C_IDHmt_cod_83.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF268T-IM-C_IDHmt_noncod_30.png", plot = P2, width = 6, height = 8, dpi = 300)




######UKF243T-WT######
UKF243T-WT <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF243_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF243T-WT) <- "Spatial"
UKF243T-WT <- NormalizeData(UKF243T-WT)
expr_matrix <- GetAssayData(
  UKF243T-WT, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHmt_cod_83
cells_AUC <- AUCell_calcAUC(
  IDHmt_cod_83, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF243T-WT@meta.data$IDHmt_cod_83 <- auc_matrix[, "cod83"]
P1<-SpatialFeaturePlot(UKF243T-WT, features = "IDHmt_cod_83",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_noncod_30
cells_AUC <- AUCell_calcAUC(
  IDHmt_noncod_30, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF243T-WT@meta.data$IDHmt_noncod_30 <- auc_matrix[, "noncod"]
P2<-SpatialFeaturePlot(UKF243T-WT, features = "IDHmt_noncod_30",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF243T-WT_IDHmt_cod_83.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF243T-WT_IDHmt_noncod_30.png", plot = P2, width = 6, height = 8, dpi = 300)



######UKF334T-WT######
UKF334T-WT <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF334_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF334T-WT) <- "Spatial"
UKF334T-WT <- NormalizeData(UKF334T-WT)
expr_matrix <- GetAssayData(
  UKF334T-WT, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHmt_cod_83
cells_AUC <- AUCell_calcAUC(
  IDHmt_cod_83, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF334T-WT@meta.data$IDHmt_cod_83 <- auc_matrix[, "cod83"]
P1<-SpatialFeaturePlot(UKF334T-WT, features = "IDHmt_cod_83",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_noncod_30
cells_AUC <- AUCell_calcAUC(
  IDHmt_noncod_30, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF334T-WT@meta.data$IDHmt_noncod_30 <- auc_matrix[, "noncod"]
P2<-SpatialFeaturePlot(UKF334T-WT, features = "IDHmt_noncod_30",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF334T-WT_IDHmt_cod_83.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF334T-WT_IDHmt_noncod_30.png", plot = P2, width = 6, height = 8, dpi = 300)
