###############################################################################
##
## Purpose:
##   Reproduce the CellChat-based cell-cell communication analyses shown in
##   Fig. 5A-G.
##
## Analyses:
##   1. Infer global intercellular communication networks in GBM, IDHMU,
##      oligodendroglioma, and astrocytoma using CellChat.
##   2. Identify ligand-receptor interactions involving subtype-associated
##      signature genes.
##   3. Visualize overall communication strength using circle plots.
##   4. Visualize selected ligand-receptor interactions using bubble plots.
##
## Required input files:
##   1. Glioma.rds
##      Annotated major-cell-type Seurat object containing glioma cells and
##      metadata columns including:
##        - group
##        - celltype
##
##      The object should include the following diagnostic groups:
##        - GBM
##        - Astrocytoma
##        - Oligodendroglioma
##
##   2. signature.xlsx
##      Gene signature file containing subtype-associated gene sets, including:
##        - IDHWT
##        - IDHMU
##        - Noncodel
##        - Codel
##
## Main outputs:
##   CellChat communication tables:
##     1. GBM_CellChat_LR_interactions.csv
##     2. GBM_CellChat_LR_interactions_net_pathway.csv
##     3. IDHMU_CellChat_LR_interactions.csv
##     4. IDHMU_CellChat_LR_interactions_net_pathway.csv
##     5. Oligodendroglioma_CellChat_LR_interactions.csv
##     6. Oligodendroglioma_CellChat_LR_interactions_net_pathway.csv
##     7. Astrocytoma_CellChat_LR_interactions.csv
##     8. Astrocytoma_CellChat_LR_interactions_net_pathway.csv
##
##   Signature-related ligand-receptor interaction tables:
##     1. Fig5E_GBM_CellChat_LR_interactions_use.csv
##     2. Fig5F_IDHMU_CellChat_LR_interactions_use.csv
##     3. Fig5G_Oligodendroglioma_CellChat_LR_interactions_use.csv
##     4. Astrocytoma_CellChat_LR_interactions_use.csv
##
##   Main figures:
##     Fig. 5A: Global CellChat communication network in GBM.
##     Fig. 5B: Global CellChat communication network in IDHMU gliomas.
##     Fig. 5C: Global CellChat communication network in astrocytoma.
##     Fig. 5D: Global CellChat communication network in oligodendroglioma.
##     Fig. 5E: Signature-related ligand-receptor interactions in GBM.
##     Fig. 5F: Signature-related ligand-receptor interactions in IDHMU gliomas.
##     Fig. 5G: Signature-related ligand-receptor interactions in oligodendroglioma.
##
## Notes:
##   IDHMU gliomas were defined as astrocytoma and oligodendroglioma in this
##   analysis. Ligand-receptor pairs were retained if either the ligand or the
##   receptor was included in the corresponding subtype-associated gene
##   signature.
##
##   The alphabetical reordering of ligand-receptor pairs in bubble plots was
##   performed only for visualization and did not modify the CellChat object,
##   communication probabilities, or statistical results.
###############################################################################

##----------------------Fig5 cell-cell communication-----------------------##
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("NMF", quietly = TRUE)) {
  install.packages("NMF")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  devtools::install_github("jokergoo/circlize")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  devtools::install_github("jokergoo/ComplexHeatmap")
}
if (!requireNamespace("CellChat", quietly = TRUE)) {
  devtools::install_github("sqjin/CellChat")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE)
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}


library(openxlsx)
library(dplyr)
library(CellChat) 

Astrocytoma <-  subset(x = data, group == "Astrocytoma")
Oligodendroglioma <- subset(x = data,  group== "Oligodendroglioma")
GBM <- subset(x = data, group == "GBM")
IDHMU <- subset( x = data,subset = group %in% c("Astrocytoma", "Oligodendroglioma"))

##---------------------------------------------------------------------------##
##------------ Fig. 5E | CellChat LR interactions in GBM -------------------##
##---------------------------------------------------------------------------##
data.input = GetAssayData(GBM, assay = "RNA", slot = "data")
meta = GBM@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 1) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")  
cellchat <- filterCommunication(cellchat, min.cells = 10) 
df.net   <- subsetCommunication(cellchat)
write.csv(df.net, "GBM_CellChat_LR_interactions.csv")

cellchat <- computeCommunProbPathway(cellchat)
df1.net1  <-  subsetCommunication(cellchat, slot.name = "netP")
write.csv(df1.net1, "GBM_CellChat_LR_interactions_net_pathway.csv")

cellchat = aggregateNet(cellchat)

signature <- read.xlsx("signature.xlsx",1)
signature <- na.omit(signature$IDHWT)

##ligand–receptor pairs containing signature genes
LR.filtered <- df.net %>%
  dplyr::filter(ligand %in% signature | receptor %in% signature)

pair.use    <- unique(LR.filtered$interaction_name)
pair.use.df <- data.frame(
  interaction_name = pair.use,
  stringsAsFactors = FALSE
)
write.csv(
  LR.filtered,
  "Fig5E_GBM_CellChat_LR_interactions_use.csv",
  row.names = FALSE
)


##Fig5E------ Pericytes as source-------------##

p <- netVisual_bubble(
  cellchat,
  sources.use = c(8),
  targets.use = c(1,2,3,4,5,6,7,8,9,10,11),
  pairLR.use = pair.use.df,
  remove.isolate = T
)

y_order <- sort(unique(as.character(p$data$interaction_name_2)))
p$data$interaction_name_2 <- factor(
  as.character(p$data$interaction_name_2),
  levels = y_order
)
p <- p + scale_y_discrete(limits = y_order)

##Fig5E---Endothelia  as source-----##
p <- netVisual_bubble(
  cellchat,
  sources.use = c(2),
  targets.use = c(1,2,3,4,5,6,7,8,9,10,11),
  pairLR.use = pair.use.df,
  remove.isolate = T
)
y_order <- sort(unique(as.character(p$data$interaction_name_2)))
p$data$interaction_name_2 <- factor(
  as.character(p$data$interaction_name_2),
  levels = y_order
)
p <- p + scale_y_discrete(limits = y_order)

##------Fig5A-----------##
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



##---------------------------------------------------------------------------##
##------------ Fig. 5F | CellChat LR interactions in IDHMU -------------------##
##---------------------------------------------------------------------------##

data.input = GetAssayData(IDHMU, assay = "RNA", slot = "data")
meta = IDHMU@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 1) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")  
cellchat <- filterCommunication(cellchat, min.cells = 10) 
df.net   <- subsetCommunication(cellchat)
write.csv(df.net, "IDHMU_CellChat_LR_interactions.csv")

cellchat <- computeCommunProbPathway(cellchat)
df1.net1  <-  subsetCommunication(cellchat, slot.name = "netP")
write.csv(df1.net1, "IDHMU_CellChat_LR_interactions_net_pathway.csv")

cellchat = aggregateNet(cellchat)

signature <- read.xlsx("signature.xlsx",1)
signature <- na.omit(signature$IDHMU)

##ligand–receptor pairs containing signature genes
LR.filtered <- df.net %>%
  dplyr::filter(ligand %in% signature | receptor %in% signature)

pair.use    <- unique(LR.filtered$interaction_name)
pair.use.df <- data.frame(
  interaction_name = pair.use,
  stringsAsFactors = FALSE
)

write.csv(
  LR.filtered,
  "Fig5F_IDHMU_CellChat_LR_interactions_us.csv",
  row.names = FALSE
)

##Fig5F------ Oligodendrocytes as source-------------##

p <- netVisual_bubble(
  cellchat,
  sources.use = c(7),
  targets.use = c(1,2,3,4,5,6,7,8,9,10,11),
  pairLR.use = pair.use.df,
  remove.isolate = T
)

y_order <- sort(unique(as.character(p$data$interaction_name_2)))
p$data$interaction_name_2 <- factor(
  as.character(p$data$interaction_name_2),
  levels = y_order
)
p <- p + scale_y_discrete(limits = y_order)


##Fig5F---Neurons  as source-----##
p <- netVisual_bubble(
  cellchat,
  sources.use = c(6),
  targets.use = c(1,2,3,4,5,6,7,8,9,10,11),
  pairLR.use = pair.use.df,
  remove.isolate = T
)
y_order <- sort(unique(as.character(p$data$interaction_name_2)))
p$data$interaction_name_2 <- factor(
  as.character(p$data$interaction_name_2),
  levels = y_order
)
p <- p + scale_y_discrete(limits = y_order)


##-----------Fig5B-------------##
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



##---------------------------------------------------------------------------##
##---------- Fig. 5G | CellChat LR interactions in  Oligodendroglioma--------##
##---------------------------------------------------------------------------##

library(CellChat)
Oligodendroglioma <- subset(x = data,  group== "Oligodendroglioma")

data.input = GetAssayData(Oligodendroglioma, assay = "RNA", slot = "data")
meta = Oligodendroglioma@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 1) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")  
cellchat <- filterCommunication(cellchat, min.cells = 10) 
df.net   <- subsetCommunication(cellchat)

write.csv(df.net, "Oligodendroglioma_CellChat_LR_interactions.csv")

cellchat <- computeCommunProbPathway(cellchat)

df1.net1  <-  subsetCommunication(cellchat, slot.name = "netP")
write.csv(df1.net1, "Oligodendroglioma_CellChat_LR_interactions_net_pathway.csv")

cellchat = aggregateNet(cellchat)
groupSize = as.numeric(table(cellchat@idents))

signature <- read.xlsx("signature.xlsx",1)
signature <- na.omit(signature$Codel)

##ligand–receptor pairs containing signature genes
LR.filtered <- df.net %>%
  dplyr::filter(ligand %in% signature | receptor %in% signature)

pair.use    <- unique(LR.filtered$interaction_name)
pair.use.df <- data.frame(interaction_name = pair.use)

write.csv(
  LR.filtered,
  "Fig5G_Oligodendroglioma_CellChat_LR_interactions_use.csv",
  row.names = FALSE
)

##Fig5G--- Neurons as source-----##
p <- netVisual_bubble(
  cellchat,
  sources.use = c(6),
  targets.use = c(1,2,3,4,5,6,7,8,9,10,11),
  pairLR.use = pair.use.df,
  remove.isolate = T
)
y_order <- sort(unique(as.character(p$data$interaction_name_2)))
p$data$interaction_name_2 <- factor(
  as.character(p$data$interaction_name_2),
  levels = y_order
)
p <- p + scale_y_discrete(limits = y_order)

##Fig5G--- OD as source----##
p <- netVisual_bubble(
  cellchat,
  sources.use = c(7),
  targets.use = c(1,2,3,4,5,6,7,8,9,10,11),
  pairLR.use = pair.use.df,
  remove.isolate = T
)
y_order <- sort(unique(as.character(p$data$interaction_name_2)))
p$data$interaction_name_2 <- factor(
  as.character(p$data$interaction_name_2),
  levels = y_order
)
p <- p + scale_y_discrete(limits = y_order)

##Fig5D----##
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")




##---------------------------------------------------------------------------##
##---------- Fig. 5C | CellChat LR interactions in  Astrocytoma--------##
##---------------------------------------------------------------------------##

library(CellChat)
Astrocytoma<- subset(x = data,  group== "Astrocytoma")

data.input = GetAssayData(Astrocytoma, assay = "RNA", slot = "data")
meta = Astrocytoma@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 1) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")  
cellchat <- filterCommunication(cellchat, min.cells = 10) 
df.net   <- subsetCommunication(cellchat)

write.csv(df.net, "Astrocytoma_CellChat_LR_interactions.csv")

cellchat <- computeCommunProbPathway(cellchat)

df1.net1  <-  subsetCommunication(cellchat, slot.name = "netP")
write.csv(df1.net1, "Astrocytoma_CellChat_LR_interactions_net_pathway.csv")

cellchat = aggregateNet(cellchat)
groupSize = as.numeric(table(cellchat@idents))

signature <- read.xlsx("signature.xlsx",1)
signature <- na.omit(signature$Noncodel)

##ligand–receptor pairs containing signature genes
LR.filtered <- df.net %>%
  dplyr::filter(ligand %in% signature | receptor %in% signature)

pair.use    <- unique(LR.filtered$interaction_name)
pair.use.df <- data.frame(interaction_name = pair.use)

write.csv(
  LR.filtered,
  "Astrocytoma_CellChat_LR_interactions_use.csv",
  row.names = FALSE
)

##Fig5C---##
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

