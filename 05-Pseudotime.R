## Pseudotime analysis - IDHWT Pericytes

library(Seurat)      
library(monocle)     
library(ggplot2)
library(AUCell)     
library(openxlsx)


signature <- read.xlsx("signature.xlsx",1)
IDHWT     <-  na.omit(signature$IDHWT)

# Subset Pericyte from GBM samples
Pericyte  <- readRDS("Pericyte_celltype.rds")
Pericyte <- subset(Pericyte, subset = sample %in% c("GBM"))

##Monocle pseudotime analysis results for Pericyte 
Cds      <- readRDS("Pericyte_monocle_result.rds")
Pericyte$Pseudotime <- 
  Biobase::pData(Cds)[colnames(Pericyte), "Pseudotime"]

p <- plot_cell_trajectory(Cds,color_by="Pseudotime") 


##AUC Score
exprMatrix <-  GetAssayData(Pericyte, assay = "RNA", slot = "counts")
cells_AUC2 <-  AUCell_run(exprMatrix, IDHWT)
AUCell_auc <-  cells_AUC2@assays@data@listData[["AUC"]]
Pericyte$AUC <- AUCell_auc
common_cells <- intersect(
  rownames(pData(Cds)),
  colnames(Pericyte)
)
Biobase::pData(Cds)[common_cells, "AUC"] <- Pericyte$AUC[common_cells]


## Monocle trajectory colored by 179-gene AUC score
p <-  plot_cell_trajectory(Cds, color_by = "AUC") +
  scale_color_viridis_c() +
  ggtitle("179gene score")


##179-gene signature score along pseudotime
library(ggplot2)
library(dplyr)


df <- data.frame(
  cell = rownames(Biobase::pData(Cds)),
  pseudotime = Biobase::pData(Cds)$Pseudotime,
  AUC = Biobase::pData(Cds)$AUC,
  celltype = Biobase::pData(Cds)$celltype
)

df <- na.omit(df)  

cell_colors <- c(
  "Pericyte_Transport"      =    "#A6CEE3",
  "Pericyte_PTPRZ1"         =    "#67a8cd",
  "Pericyte_ECM"            =    "#B3E2D5",
  "Pericyte_PGF"            =    "#EDCFD8",
  "Pericyte_Proliferation"  =    "#E4CEBB",
)

ggplot(df, aes(x = pseudotime, y = AUC)) +
  geom_point(aes(color = celltype), alpha = 0.9, size = 1)+ 
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  scale_color_manual(values = cell_colors) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加方框
    panel.background = element_blank(),  # 确保背景干净
    plot.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),  # 刻度线
    axis.ticks.length = unit(0.2, "cm"),  # 刻度线长度
  ) +labs(
    title = "179 Gene Score",
    x = "Pseudotime",
    y = "179 Gene Score",
    color = "celltype"
  )




