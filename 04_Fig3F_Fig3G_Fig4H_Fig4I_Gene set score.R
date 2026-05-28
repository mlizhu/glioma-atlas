############################################################
## Script: Fig3F_Fig3G_Fig4H_Fig4I_Gene_set_score.R
## The scores corresponding to the glioma subtypes of gliomas and umap
## Purpose:
##   Calculate AUCell scores for glioma subtype-related gene sets
##   and reproduce the UMAP and gene-set score plots shown in
##   Fig. 1D, Fig. 3F,Fig. 3G, Fig. 4H, Fig. 4I.
##
## Required input files:
##   1. Glioma.rds
##      Annotated major-cell-type Seurat object containing:
##      - RNA assay
##      - celltype annotation
##      - group annotation
##      - umap.harmony reduction
##
##   2. signature.xlsx
##      Gene signatures including:
##      - IDHWT
##      - IDHMU
##      - Noncodel
##      - Codel
##
## Main outputs:
##   Fig. 1D: UMAP of all glioma cells
##   Fig. 3F: UMAP of GBM cells
##   Fig. 3F: UMAP of IDH-mutant glioma cells
##   Fig. 4H: UMAP of astrocytoma cells
##   Fig. 4H: UMAP of oligodendroglioma cells
##   Fig. 4I: 30-gene score in astrocytoma
##   Fig. 4I: 83-gene score in oligodendroglioma
############################################################

library(AUCell)
library(Matrix)
library(ggplot2)
library(scales)
library(Seurat)
library(openxlsx)

cell_colors <- c(
  "Tum"         = "#A6CEE3",
  "Pro.Tum"     = "#67a8cd",
  "OD"          = "#AF98B5",
  "T"           = "#EEBDD1",
  "Mφ"          = "#E4CEBB",
  "MD"          = "#cf9f88",
  "Peri"        = "#B3E2D5",
  "MO"          = "#FDBF6F",
  "DC"          = "#62B2A1",
  "Endo"        = "#17BECF",
  "Neu"         = "#f36569"
)



##Read the 179, 68, 30, and 83 gene sets
signature  <- read.xlsx("signature.xlsx",1)
gene_179<- na.omit(signature$IDHWT)
gene_68 <- na.omit(signature$IDHMU)
gene_30 <- na.omit(signature$Noncodel)
gene_83 <- na.omit(signature$Codel)


##Read the annotated major subgroups
data <- readRDS("Glioma.rds")
expr <- GetAssayData(data, assay = "RNA", slot = "counts")
gene_sets <- list(
  gene_179 = gene_179,   
  gene_68  = gene_68,
  gene_30  = gene_30,
  gene_83  = gene_83
)


cells_AUC <- AUCell_run(expr, gene_sets)
auc_mat <- t(getAUC(cells_AUC))[colnames(data), ]
data@meta.data <- cbind(data@meta.data, auc_mat)

GBM <- subset(data, group %in% c("GBM"))
OD  <- subset(data, group %in% c("Oligodendroglioma"))
AC  <- subset(data, group %in% c("Astrocytoma"))
IDHMU  <- subset(data, group %in% c("Astrocytoma","Oligodendroglioma"))


##Fig.1D
p_ALL <- DimPlot(data,reduction = "umap.harmony",group.by = "celltype",label = T,cols = cell_colors,raster = TRUE)

##Fig.3F
P_GBM <- DimPlot(GBM,reduction = "umap.harmony",group.by = "celltype",label = T, cols = cell_colors,raster = TRUE)

##Fig.3F
p_IDHMU  <- DimPlot(IDHMU, reduction = "umap.harmony",group.by = "celltype",label = T,cols = cell_colors,raster = TRUE)

##Fig.4H
P_AC  <- DimPlot(AC, reduction = "umap.harmony",group.by = "celltype",label = T,cols = cell_colors,raster = TRUE)

##Fig.4H
P_OD  <- DimPlot(OD, reduction = "umap.harmony",group.by = "celltype",label = T,cols = cell_colors,raster = TRUE)

##----------------------------------------------------------------------##

##Fig.3G --------The 179 gene score of GBM————————------————————————————##

##-----------------------------------------------------------------------##

umap_df <- as.data.frame(GBM@reductions$umap.harmony@cell.embeddings)
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
plot_df <- data.frame(
  GBM@meta.data,
  umap_df
)
p_GBM_score <- ggplot(
  plot_df,
  aes(x = umap_1, y = umap_2, color = gene_179)
) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(
    option = "A",
    limits = c(0, 0.3),
    oob = scales::squish,
    name = "AUC"
  ) +
  labs(
    title = "179-gene score",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_light(base_size = 7, base_line_size = 0.5) +
  theme(
    text = element_text(colour = "black", size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_text(colour = "black", size = 12),
    plot.title = element_text(
      hjust = 0.5,
      colour = "black",
      size = 12
    )
  )

##--------------------------------------------------------------##

##Fig.3G --  The 68 gene score of IDHMU————————————————————————##

##-------------------------------------------------------------##

umap_df <- as.data.frame(IDHMU@reductions$umap.harmony@cell.embeddings)
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
plot_df <- data.frame(
  IDHMU@meta.data,
  umap_df
)
p_IDHMU_score <- ggplot(
  plot_df,
  aes(x = umap_1, y = umap_2, color = gene_68)
) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(
    option = "A",
    limits = c(0, 0.3),
    oob = scales::squish,
    name = "AUC"
  ) +
  labs(
    title = "68-gene score",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_light(base_size = 7, base_line_size = 0.5) +
  theme(
    text = element_text(colour = "black", size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_text(colour = "black", size = 12),
    plot.title = element_text(
      hjust = 0.5,
      colour = "black",
      size = 12
    )
  )

##---------------------------------------------------------##

##----Fig.4I --  The 30 gene score of AC-------------------##

##----------------------------------------------------------##


umap_df <- as.data.frame(AC@reductions$umap.harmony@cell.embeddings)
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
plot_df <- data.frame(
  AC@meta.data,
  umap_df
)
p_AC_score <- ggplot(
  plot_df,
  aes(x = umap_1, y = umap_2, color = gene_30)
) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(
    option = "A",
    limits = c(0, 0.3),
    oob = scales::squish,
    name = "AUC"
  ) +
  labs(
    title = "30-gene score",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_light(base_size = 7, base_line_size = 0.5) +
  theme(
    text = element_text(colour = "black", size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_text(colour = "black", size = 12),
    plot.title = element_text(
      hjust = 0.5,
      colour = "black",
      size = 12
    )
  )

##----------------------------------------------------##

## Fig.4I --  The 83 gene score of OD—----------------##

##-----------------------------------------------------##

umap_df <- as.data.frame(OD@reductions$umap.harmony@cell.embeddings)
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
plot_df <- data.frame(
  OD@meta.data,
  umap_df
)
p_OD_score <- ggplot(
  plot_df,
  aes(x = umap_1, y = umap_2, color = gene_83)
) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(
    option = "A",
    limits = c(0, 0.3),
    oob = scales::squish,
    name = "AUC"
  ) +
  labs(
    title = "83-gene score",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_light(base_size = 7, base_line_size = 0.5) +
  theme(
    text = element_text(colour = "black", size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "lines"),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_text(colour = "black", size = 12),
    plot.title = element_text(
      hjust = 0.5,
      colour = "black",
      size = 12
    )
  )
