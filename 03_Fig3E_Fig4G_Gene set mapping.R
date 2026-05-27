##############################################################################
## Script: Fig3E_Fig4G_Gene_set_mapping.R
## 179 gene, 68gene, 30gene, 83gene mapping
## Purpose:
##   Reproduce the gene-set mapping heatmaps shown in Fig. 3E and Fig. 4G.
##
## Required input files:
##   1. Glioma.rds
##      Annotated major-cell-type Seurat object.
##   2. signature.xlsx
##      Gene sets including IDHWT, IDHMU, Noncodel, and Codel signatures.
##
## Main outputs:
##   Fig. 3E: 179-gene mapping heatmap
##   Fig. 3E: 68-gene mapping heatmap
##   Fig. 4G: 30-gene mapping heatmap
##   Fig. 4G: 83-gene mapping heatmap
##
## Notes:
##   Neurons were retained in the annotated Seurat object, but excluded
##   only from the downstream gene-mapping analysis because of their
##   limited cell number.
###############################################################################


if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE)
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(Seurat)
library(openxlsx)
library(tidyverse)

##Read the annotated major subgroups
data <- readRDS("Glioma.rds")

## Neurons were retained in the annotated Seurat object.
## Due to their limited cell number, neurons were excluded only from
## the downstream gene-mapping analysis to avoid unstable estimates.

data <- subset(x = data, celltype != "Neu")

##Read the 179, 68, 30, and 83 gene sets
signature  <- read.xlsx("signature.xlsx",1)
signature1 <- na.omit(signature$IDHWT)
signature2 <- na.omit(signature$IDHMU)
signature3 <- na.omit(signature$Noncodel)
signature4 <- na.omit(signature$Codel)

###--------------------------------------------------------------------##
#----------Fig.3E----- 179 gene mapping---------------------------------#
##---------------------------------------------------------------------##
signature <- as.data.frame(signature1)
names(signature ) <- "SYMBOL"

cluster.ave <- AverageExpression(data,group.by = "celltype",layer = "data")
cluster.ave1  <- cluster.ave 
cluster.ave <- cluster.ave1[["RNA"]]

cluster.ave <- as.data.frame(cluster.ave)
cluster.ave <- cluster.ave %>% rownames_to_column(var = "SYMBOL") 
cluster.ave <- dplyr::left_join(signature,cluster.ave,by="SYMBOL")
cluster.ave <- na.omit(cluster.ave)
rownames(cluster.ave) <- cluster.ave[,1]
cluster.ave <- cluster.ave[,-1]

highest_celltype <- apply(cluster.ave, 1, function(gene_expr) {
  cell_type <- names(which.max(gene_expr))
  return(cell_type)
})
sorted_celltype <- table(highest_celltype) %>% sort(decreasing = TRUE) %>% names()
num_celltype <- length(sorted_celltype)
expr <- cluster.ave[,sorted_celltype]
highest_celltype <- as.data.frame(highest_celltype)
highest_celltype <- highest_celltype %>% rownames_to_column (var = "signature")

df <- expr %>% rownames_to_column (var = "signature")
df <- left_join(highest_celltype,df,by="signature")

df1 <- NULL
df2 <- NULL
for (v in sorted_celltype) {
  df1[[v]] <- df[df$highest_celltype == v, ] 
}
for (i in 1:num_celltype){
  df2[[i]] <- df1[[i]][order(-df1[[i]][,i+2]),]
}
names(df2) <- sorted_celltype
df3 <- as.data.frame(data.table::rbindlist(df2))
rownames(df3) <- df3$signature
df3 <- df3[,c(-1,-2)]
colours <- colorRampPalette(c("#64a5f8","#FFFFFF","#e7110f"))(100)
ph  <- pheatmap::pheatmap(df3,
                          border_color = NA,
                          border=F,
                          color = colours,
                          scale = "row",
                          fontsize_row = 10,
                          fontsize_col = 10,
                          drop_levels=F,
                          show_rownames = F,
                          show_colnames = T,
                          cluster_cols = F,
                          cluster_rows = F
)

##-----------------------------------------------------------------------------##

###Fig 3E---------------------------68 gene mapping-----------------------------##

##------------------------------------------------------------------------------##


signature <- as.data.frame(signature2)
names(signature ) <- "SYMBOL"

cluster.ave <- AverageExpression(data,group.by = "celltype",layer = "data")
cluster.ave1  <- cluster.ave 
cluster.ave <- cluster.ave1[["RNA"]]
cluster.ave <- as.data.frame(cluster.ave)
cluster.ave <- cluster.ave %>% rownames_to_column(var = "SYMBOL") 


cluster.ave <- dplyr::left_join(signature,cluster.ave,by="SYMBOL")
cluster.ave <- na.omit(cluster.ave)
rownames(cluster.ave) <- cluster.ave[,1]
cluster.ave <- cluster.ave[,-1]

highest_celltype <- apply(cluster.ave, 1, function(gene_expr) {
  cell_type <- names(which.max(gene_expr))
  return(cell_type)
})
sorted_celltype <- table(highest_celltype) %>% sort(decreasing = TRUE) %>% names()
num_celltype <- length(sorted_celltype)

expr <- cluster.ave[,sorted_celltype]
highest_celltype <- as.data.frame(highest_celltype)
highest_celltype <- highest_celltype %>% rownames_to_column (var = "signature")

df <- expr %>% rownames_to_column (var = "signature")
df <- left_join(highest_celltype,df,by="signature")

df1 <- NULL
df2 <- NULL
for (v in sorted_celltype) {
  df1[[v]] <- df[df$highest_celltype == v, ] 
}
for (i in 1:num_celltype){
  df2[[i]] <- df1[[i]][order(-df1[[i]][,i+2]),]
}
names(df2) <- sorted_celltype
df3 <- as.data.frame(data.table::rbindlist(df2))

## IDHMU supplements cells without the highest gene expression: MG and Mφ
cluster.ave_1 <- cluster.ave %>% rownames_to_column(var = "signature") 
Def<- cluster.ave_1[,c(1,4,6)]
df4 <- left_join(df3,Def,by="signature")
rownames(df4) <- df4$signature
df4 <- df4[,c(-1,-2)]

ph  <- pheatmap::pheatmap(df4,
                          border_color = NA,
                          border=F,
                          color = colours,
                          scale = "row",
                          fontsize_row = 10,
                          fontsize_col = 10,
                          drop_levels=F,
                          show_rownames = F,
                          show_colnames = T,
                          cluster_cols = F,
                          cluster_rows = F
)


##------------------------------------------------------------------------------##

###Fig 4G---------------------------30 gene mapping-----------------------------##

##------------------------------------------------------------------------------##

signature <- as.data.frame(signature3)
names(signature ) <- "SYMBOL"

cluster.ave <- AverageExpression(data,group.by = "celltype",layer = "data")
cluster.ave1  <- cluster.ave 
cluster.ave <- cluster.ave1[["RNA"]]
cluster.ave <- as.data.frame(cluster.ave)
cluster.ave <- cluster.ave %>% rownames_to_column(var = "SYMBOL") 


cluster.ave <- dplyr::left_join(signature,cluster.ave,by="SYMBOL")
cluster.ave <- na.omit(cluster.ave)
rownames(cluster.ave) <- cluster.ave[,1]
cluster.ave <- cluster.ave[,-1]

highest_celltype <- apply(cluster.ave, 1, function(gene_expr) {
  cell_type <- names(which.max(gene_expr))
  return(cell_type)
})

sorted_celltype <- table(highest_celltype) %>% sort(decreasing = TRUE) %>% names()
num_celltype <- length(sorted_celltype)
expr <- cluster.ave[,sorted_celltype]
highest_celltype <- as.data.frame(highest_celltype)
highest_celltype <- highest_celltype %>% rownames_to_column (var = "signature")

df <- expr %>% rownames_to_column (var = "signature")
df <- left_join(highest_celltype,df,by="signature")

df1 <- NULL
df2 <- NULL
for (v in sorted_celltype) {
  df1[[v]] <- df[df$highest_celltype == v, ] 
}
for (i in 1:num_celltype){
  df2[[i]] <- df1[[i]][order(-df1[[i]][,i+2]),]
}
names(df2) <- sorted_celltype
df3 <- as.data.frame(data.table::rbindlist(df2))

##AC supplements cells without the highest gene expression: Endo and Mφ
cluster.ave_1 <- cluster.ave %>% rownames_to_column(var = "signature") 
Def<- cluster.ave_1[,c(1,3,6)]
df4 <- left_join(df3,Def,by="signature")
rownames(df4) <- df4$signature
df4 <- df4[,c(-1,-2)]

ph  <- pheatmap::pheatmap(df4,
                          border_color = NA,
                          border=F,
                          color = colours,
                          scale = "row",
                          fontsize_row = 10,
                          fontsize_col = 10,
                          drop_levels=F,
                          show_rownames = F,
                          show_colnames = T,
                          cluster_cols = F,
                          cluster_rows = F
)

##----------------------------------------------------------------------------##

###Fig 4G---------------------------83 gene mapping---------------------------##

##----------------------------------------------------------------------------##

signature <- as.data.frame(signature4)
names(signature ) <- "SYMBOL"

cluster.ave <- AverageExpression(data,group.by = "celltype",layer = "data")
cluster.ave1  <- cluster.ave 
cluster.ave <- cluster.ave1[["RNA"]]
cluster.ave <- as.data.frame(cluster.ave)
cluster.ave <- cluster.ave %>% rownames_to_column(var = "SYMBOL") 


cluster.ave <- dplyr::left_join(signature,cluster.ave,by="SYMBOL")
cluster.ave <- na.omit(cluster.ave)
rownames(cluster.ave) <- cluster.ave[,1]
cluster.ave <- cluster.ave[,-1]

highest_celltype <- apply(cluster.ave, 1, function(gene_expr) {
  cell_type <- names(which.max(gene_expr))
  return(cell_type)
})
sorted_celltype <- table(highest_celltype) %>% sort(decreasing = TRUE) %>% names()
num_celltype <- length(sorted_celltype)
expr <- cluster.ave[,sorted_celltype]
highest_celltype <- as.data.frame(highest_celltype)
highest_celltype <- highest_celltype %>% rownames_to_column (var = "signature")

df <- expr %>% rownames_to_column (var = "signature")
df <- left_join(highest_celltype,df,by="signature")

df1 <- NULL
df2 <- NULL
for (v in sorted_celltype) {
  df1[[v]] <- df[df$highest_celltype == v, ] 
}
for (i in 1:num_celltype){
  df2[[i]] <- df1[[i]][order(-df1[[i]][,i+2]),]
}
names(df2) <- sorted_celltype
df3 <- as.data.frame(data.table::rbindlist(df2))

##OD supplements cells without the highest gene expression: Endo and Mφ
cluster.ave_1 <- cluster.ave %>% rownames_to_column(var = "signature") 
Def<- cluster.ave_1[,c(1,6)]
df4 <- left_join(df3,Def,by="signature")
rownames(df4) <- df4$signature
df4 <- df4[,c(-1,-2)]
ph  <- pheatmap::pheatmap(df4,
                          border_color = NA,
                          border=F,
                          color = colours,
                          scale = "row",
                          fontsize_row = 10,
                          fontsize_col = 10,
                          drop_levels=F,
                          show_rownames = F,
                          show_colnames = T,
                          cluster_cols = F,
                          cluster_rows = F
)











