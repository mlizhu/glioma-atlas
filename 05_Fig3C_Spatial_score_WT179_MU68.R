Figure3C:
library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library('RColorBrewer')
######
col1<-c("grey90", "lightpink", "red", "darkred")
pal1<-colorRampPalette(col1)
col2<-c("grey90", "lightblue", "dodgerblue", "darkblue")
pal2<-colorRampPalette(col2)
IDHwt_179 <- read.table("/home/Cancer_Cell_10xSpatial/genelist/IDHwt_179.txt",header=T)
IDHwt_179 <- as.list(IDHwt_179)
IDHmt_68 <- read.table("/home/Cancer_Cell_10xSpatial/genelist/IDHmt_68.txt",header=T)
IDHmt_68 <- as.list(IDHmt_68)

######UKF270T-IDH-M-N######
UKF270T-IDH-M-N <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF270_IDHMutant_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF270T-IDH-M-N) <- "Spatial"
UKF270T-IDH-M-N <- NormalizeData(UKF270T-IDH-M-N)
expr_matrix <- GetAssayData(
  UKF270T-IDH-M-N, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHwt_179
cells_AUC <- AUCell_calcAUC(
  IDHwt_179, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF270T-IDH-M-N@meta.data$IDHwt_179 <- auc_matrix[, "IDHwt179"]
P1<-SpatialFeaturePlot(UKF270T-IDH-M-N, features = "IDHwt_179",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_68
cells_AUC <- AUCell_calcAUC(
  IDHmt_68, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF270T-IDH-M-N@meta.data$IDHmt_68 <- auc_matrix[, "IDHmt68"]
P2<-SpatialFeaturePlot(UKF270T-IDH-M-N, features = "IDHmt_68",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF270T-IDH-M-N_IDHwt_179.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF270T-IDH-M-N_IDHmt_68.png", plot = P2, width = 6, height = 8, dpi = 300)




######UKF268T-IDH-M-C######
UKF268T-IDH-M-C <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF268_IDHMutant_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF268T-IDH-M-C) <- "Spatial"
UKF268T-IDH-M-C <- NormalizeData(UKF268T-IDH-M-C)
expr_matrix <- GetAssayData(
  UKF268T-IDH-M-C, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHwt_179
cells_AUC <- AUCell_calcAUC(
  IDHwt_179, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF268T-IDH-M-C@meta.data$IDHwt_179 <- auc_matrix[, "IDHwt179"]
P1<-SpatialFeaturePlot(UKF268T-IDH-M-C, features = "IDHwt_179",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_68
cells_AUC <- AUCell_calcAUC(
  IDHmt_68, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF268T-IDH-M-C@meta.data$IDHmt_68 <- auc_matrix[, "IDHmt68"]
P2<-SpatialFeaturePlot(UKF268T-IDH-M-C, features = "IDHmt_68",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF268T-IDH-M-C_IDHwt_179.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF268T-IDH-M-C_IDHmt_68.png", plot = P2, width = 6, height = 8, dpi = 300)




######UKF334T-IDH-WT######
UKF334T-IDH-WT <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF334_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF334T-IDH-WT) <- "Spatial"
UKF334T-IDH-WT <- NormalizeData(UKF334T-IDH-WT)
expr_matrix <- GetAssayData(
  UKF334T-IDH-WT, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHwt_179
cells_AUC <- AUCell_calcAUC(
  IDHwt_179, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF334T-IDH-WT@meta.data$IDHwt_179 <- auc_matrix[, "IDHwt179"]
P1<-SpatialFeaturePlot(UKF334T-IDH-WT, features = "IDHwt_179",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_68
cells_AUC <- AUCell_calcAUC(
  IDHmt_68, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF334T-IDH-WT@meta.data$IDHmt_68 <- auc_matrix[, "IDHmt68"]
P2<-SpatialFeaturePlot(UKF334T-IDH-WT, features = "IDHmt_68",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF334T-IDH-WT_IDHwt_179.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF334T-IDH-WT_IDHmt_68.png", plot = P2, width = 6, height = 8, dpi = 300)




######UKF243T-IDH-WT######
UKF243T-IDH-WT <- Load10X_Spatial(
  '/home/Cancer_Cell_10xSpatial/#UKF243_T_ST/outs/',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE
)
DefaultAssay(UKF243T-IDH-WT) <- "Spatial"
UKF243T-IDH-WT <- NormalizeData(UKF243T-IDH-WT)
expr_matrix <- GetAssayData(
  UKF243T-IDH-WT, 
  assay = "Spatial", 
  layer = "data"
)
set.seed(42)
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,    
  splitByBlocks = TRUE
)
#######IDHwt_179
cells_AUC <- AUCell_calcAUC(
  IDHwt_179, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF243T-IDH-WT@meta.data$IDHwt_179 <- auc_matrix[, "IDHwt179"]
P1<-SpatialFeaturePlot(UKF243T-IDH-WT, features = "IDHwt_179",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal1(n=100))
#######IDHmt_68
cells_AUC <- AUCell_calcAUC(
  IDHmt_68, 
  cells_rankings,
  aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
  nCores = 4
)
auc_matrix <- t(getAUC(cells_AUC))
UKF243T-IDH-WT@meta.data$IDHmt_68 <- auc_matrix[, "IDHmt68"]
P2<-SpatialFeaturePlot(UKF243T-IDH-WT, features = "IDHmt_68",pt.size.factor= 2.4, stroke=0, image.alpha = 0) & scale_fill_gradientn(limits=c(0, 0.5), oob = scales::squish, colours=pal2(n=100))
######save_plot
ggsave("UKF243T-IDH-WT_IDHwt_179.png", plot = P1, width = 6, height = 8, dpi = 300)
ggsave("UKF243T-IDH-WT_IDHmt_68.png", plot = P2, width = 6, height = 8, dpi = 300)