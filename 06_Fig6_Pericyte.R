###############################################################################
## Script: Fig6_Pericyte_subcluster_analysis.R
##
## Purpose:
##   Reproduce the pericyte subcluster analyses shown in Fig. 6.
##
## Analyses:
##   1. Visualize pericyte subclusters using tSNE.
##   2. Compare the proportions of pericyte subclusters across glioma subtypes.
##   3. Calculate Ro/e scores to evaluate the relative enrichment of pericyte
##      subclusters across tumor types.
##   4. Perform BayesPrism-based bulk RNA-seq deconvolution to infer pericyte
##      subcluster fractions in external glioma bulk RNA-seq samples.
##   5. Evaluate the survival association of BayesPrism-inferred pericyte
##      subcluster fractions.
##   6. Visualize precomputed AUCell-based subtype-associated gene signature
##      scores in pericytes.
##   7. Perform Monocle2 trajectory analysis for pericytes from GBM,
##      astrocytoma, and oligodendroglioma.
##   8. Visualize subtype-associated gene signature scores along pericyte
##      pseudotime trajectories.
##
## Required input files:
##   1. Pericyte.rds
##      Annotated pericyte Seurat object containing pericyte subclusters and
##      metadata columns including:
##        - group
##        - celltype
##        - tumor_type
##
##      The object should include pericytes from the following diagnostic groups:
##        - GBM
##        - Astrocytoma
##        - Oligodendroglioma
##
##      The tumor_type column should contain:
##        - WT
##        - M-N
##        - M-C
##
##   2. signature.xlsx
##      Gene signature file containing subtype-associated gene sets, including:
##        - IDHWT
##        - IDHMU
##        - Noncodel
##        - Codel
##
##   3. TCGA_GBM_LGG_bulk_count_matrix_with_survival.xlsx
##      Bulk RNA-seq expression matrix with survival information used for
##      BayesPrism-based deconvolution and survival analysis.
##
## Main outputs:
##   Main figures:
##     Fig. 6A: tSNE visualization of pericyte subclusters.
##     Fig. 6B: Cell proportion of pericyte subclusters across tumor types.
##     Fig. 6C: Ro/e enrichment heatmap of pericyte subclusters across tumor types.
##     Fig. 6D: BayesPrism-inferred pericyte subcluster fractions and survival
##              associations.
##     Fig. 6E: AUCell-based 179-, 68-, 30-, and 83-gene signature scores in
##              pericytes.
##     Fig. 6F: Pseudotime density distribution of pericyte subclusters in GBM,
##              astrocytoma, and oligodendroglioma.
##     Fig. 6G: Subtype-associated gene signature scores along pericyte
##              pseudotime.
##     Fig. 6H: Monocle trajectories colored by subtype-associated gene
##              signature scores.
##
##   Output files:
##     1. BayesPrism-inferred pericyte fraction.csv
##        BayesPrism-inferred fractions of pericyte subclusters in bulk RNA-seq
##        samples, combined with survival metadata.
##
##     2. survival/*.pdf
##        Kaplan-Meier survival plots stratified by the median value of each
##        BayesPrism-inferred pericyte subcluster fraction.
##
##     3. Pericyte_GBM_monocle_cds.rds
##        Monocle2 CDS object for GBM-derived pericytes.
##
##     4. Pericyte_AC_monocle_cds.rds
##        Monocle2 CDS object for astrocytoma-derived pericytes.
##
##     5. Pericytes_OD_monocle_cds.rds
##        Monocle2 CDS object for oligodendroglioma-derived pericytes.
##
## Notes:
##   Pericyte subclusters were analyzed using the annotated pericyte Seurat
##   object. Ro/e scores were calculated using tumor_type as the tissue/group
##   variable and celltype as the pericyte subcluster variable.
##
##   For Fig. 6D, BayesPrism was used to deconvolute external bulk RNA-seq
##   profiles and infer the fractions of pericyte subclusters. These inferred
##   fractions were further used for Kaplan-Meier survival analysis.
##
##   For trajectory analysis, pericytes were analyzed separately within GBM,
##   astrocytoma, and oligodendroglioma. Ordering genes were selected based on
##   differential expression across pericyte subclusters using Monocle2.
##
##   Subtype-associated gene signatures were scored along pseudotime using:
##     - IDHWT signature for GBM-derived pericytes
##     - Noncodel signature for astrocytoma-derived pericytes
##     - Codel signature for oligodendroglioma-derived pericytes
###############################################################################


if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("Startrac", quietly = TRUE)) {
  remotes::install_github("Japrin/STARTRAC", quiet = TRUE)
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(dplyr)
library(ggplot2)
library(openxlsx)
library(Seurat)
library(Startrac)
library(RColorBrewer)
library(pheatmap)


Pericyte <- readRDS("Pericyte.rds")

cell_colors <- c(
  "Transport"      =    "#A6CEE3",
  "PTPRZ1"         =    "#67a8cd",
  "ECM"            =    "#B3E2D5",
  "PGF"            =    "#EDCFD8",
  "Proliferation"  =    "#E4CEBB"
)
Pericyte$tumor_type <- factor(
  Pericyte$tumor_type,
  levels = c("WT", "M-N", "M-C")
)


##--Fig6A--------tsne-----------------##

p <- DimPlot(Pericyte, reduction = "tsne",label=T,cols = cell_colors)

##--Fig6B----- cell proportion-----------##

meta.data <- Pericyte@meta.data
meta.data$celltype <- factor(
  meta.data$celltype,
  levels = c(
    "ECM",
    "PGF",
    "Proliferation",
    "PTPRZ1",
    "Transport"
  )
)

sub_data <- table(meta.data$tumor_type, meta.data$celltype) %>% 
  as.data.frame() %>% 
  group_by(Var1) %>% 
  mutate(
    cell_proportion = Freq/sum(Freq)
  )
colnames(sub_data) <- c("sample", "cluster", "cell_num", "cell_proportion")

p <- ggplot(sub_data)+
  geom_col(aes(x = sample, y = cell_proportion, fill = cluster), position = 'stack', width = 0.6)+
  geom_hline(aes(yintercept = 0.5), linetype = 2)+
  theme_bw()+
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(
      angle = 90
    )
  )+scale_fill_manual(values = cell_colors)


##-Fig6C-------Roe Score-----------##

Roe <- calTissueDist(meta.data,
                     colname.cluster = "celltype", 
                     colname.tissue = "tumor_type", 
                     method = "chisq", 
                     min.rowSum = 0) 

colours  <- colorRampPalette(brewer.pal(9, "OrRd"))(50)

ph  <- pheatmap(Roe,
                border_color = "black",
                border=F,
                color = colours,
                scale = 'none',
                fontsize_row = 12,
                fontsize_col = 6,
                #angle_col = 45,
                drop_levels=F,
                show_rownames = T,
                show_colnames = T,
                cluster_cols = F,
                cluster_rows = F,
                display_numbers = TRUE,
                cellwidth = 20, 
                cellheight = 20,
                
)


##--------------------------------------------------------------------------##

##Fig6D-----------Bulk RNA-seq deconvolution using BayesPrism---------------##

##---------------------------------------------------------------------------##
if (!requireNamespace("survminer", quietly = TRUE)) {
  install.packages("survminer")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")
}

library(devtools)
library(BayesPrism)
library(survival)
library(survminer)

species = "hs"  
sc.dat <-  as.matrix(t(GetAssayData(Pericyte, assay = "RNA", slot = "counts")))
cell.type.labels<- Pericyte@meta.data$celltype  
cell.state.labels<-NULL   

#1,QC of cell type and state labels
plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              cexRow=0.5, cexCol=0.5,
)

#2,Filter outlier genes
sc.stat <- plot.scRNA.outlier(
  input=sc.dat,
  cell.type.labels=cell.type.labels,
  species=species, 
  return.raw=TRUE 
)

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species=species, 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

exp <- read.xlsx("整理过的TCGAGBMLGG全数据_需要log2.xlsx")
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- as.data.frame(t(exp))

metadata <- exp[,c(4:5)]
names(metadata) <- c("OS.time","OS")
metadata$OS     <- as.numeric (metadata$OS)
metadata$OS.time <- as.numeric (metadata$OS.time)


n <- ncol(exp)
exp<- exp[,c(13:n)]
exp <- as.data.frame(t(exp))
Symbol <- rownames(exp)
Sample  <- names(exp)
exp <- as.data.frame(lapply(exp,as.numeric))
rownames(exp) <- Symbol

bk.dat<-as.matrix(t(exp))

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,
  sc.input=sc.dat, 
  cell.type.labels=cell.type.labels,
  species=species, 
  return.raw=TRUE
)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
)

myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = NULL,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1)

bp.res <- run.prism(prism = myPrism, n.cores=30)

score <- get.fraction(bp=bp.res,
                      which.theta="final",
                      state.or.type="type")
score <- as.data.frame(score)
score <- rownames_to_column(score, "id")

metadata <- rownames_to_column(metadata, "id")
metadata_score <-left_join(metadata,score,by="id")
write.csv(metadata_score,"BayesPrism-inferred pericyte fraction.csv",row.names = F)


dir.create("survival", showWarnings = FALSE)
for (i in colnames(metadata_score)[4:ncol(metadata_score)]) {
  rt <- metadata_score
  med <- median(rt[[i]], na.rm = TRUE)
  rt$group <- ifelse(rt[[i]] > med, "High", "Low")
  rt$group <- factor(rt$group, levels = c("High", "Low"))
  diff <- survdiff(Surv(OS.time, OS) ~ group, data = rt)
  p <- 1 - pchisq(diff$chisq, df = 1)
  
  if (p < 0.0001) {
    pValue <- "p < 0.0001"
  } else {
    pValue <- paste0("p = ", sprintf("%.4f", p))
  }
  
  fit <- survfit(Surv(OS.time, OS) ~ group, data = rt)
  surPlot <- ggsurvplot(
    fit,
    data = rt,
    pval = pValue,
    risk.table = FALSE,
    ggtheme = theme_classic() +
      theme(
        axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black"),
        panel.grid = element_blank()
      ),
    censor.shape = 124,
    censor.size = 2,
    conf.int = FALSE,
    legend.labs = c("High", "Low"),
    legend.title = i,
    xlab = "Overall survival (months)",
    ylab = "Survival probability",
    palette = c("#E31A1C", "#1F78B4")
  )
  pdf(
    paste0("survival/", i, ".pdf"),
    onefile = FALSE,
    width = 5.5,
    height = 5
  )
  print(surPlot)
  dev.off()
}

##--------------------------------------------------------------------------##

##Fig6E---------------------179,68,30,83 gene score----------------------##

##---------------------------------------------------------------------------##

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("AUCell", quietly = TRUE)) {
  BiocManager::install("AUCell")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

library(AUCell)
library(Matrix)
library(ggplot2)
library(scales)

signature  <- read.xlsx("signature.xlsx",1)

gene_179<- na.omit(signature$IDHWT)
gene_68 <- na.omit(signature$IDHMU)
gene_30 <- na.omit(signature$Noncodel)
gene_83 <- na.omit(signature$Codel)


expr <- GetAssayData(Pericyte, assay = "RNA", slot = "counts")
gene_sets <- list(
  gene_179 = gene_179,   
  gene_68  = gene_68,
  gene_30  = gene_30,
  gene_83  = gene_83
)
cells_AUC <- AUCell_run(expr, gene_sets)
auc_mat <- t(getAUC(cells_AUC))[colnames(Pericyte), ]
Pericyte@meta.data <- cbind(Pericyte@meta.data, auc_mat)

##The drawing code is already available in Gene set score


##--------------------------------------------------------------------------##

##Fig6F,G,H------ GBM,Astrocytoma,Oligodendroglioma--monocle-----------------##

##---------------------------------------------------------------------------##

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("monocle", quietly = TRUE)) {
  BiocManager::install("monocle", ask = FALSE, update = FALSE)
}
if (!requireNamespace("AUCell", quietly = TRUE)) {
  BiocManager::install("AUCell", ask = FALSE, update = FALSE)
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = TRUE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}


library(monocle)
library(Seurat)
library(ggplot2)
library(AUCell)     
library(openxlsx)
library(dplyr)
library(ggridges)

##-------------------01GBM-------------------------##
Pericytes_GBM     <-  subset(Pericyte, subset = group %in% c("GBM"))
data <- GetAssayData(Pericytes_GBM , assay = "RNA", slot = "counts")
data <- as(as.matrix(data),'sparseMatrix')


pd <- new('AnnotatedDataFrame', data = Pericytes_GBM@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)


HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM,cores=30)
HSMM <- detectGenes(HSMM, min_expr = 0.1 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
head(pData(HSMM))


diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],      
                                      fullModelFormulaStr = "~ celltype",
                                      cores = 30)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) 


HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)
table(pData(HSMM)$State,pData(HSMM)$celltype)     

cds <- orderCells(HSMM, root_state = 7)

Pericytes_GBM$Pseudotime <- 
  Biobase::pData(cds)[colnames(Pericytes_GBM), "Pseudotime"]

signature <- read.xlsx("signature.xlsx",1)
IDHWT     <-  na.omit(signature$IDHWT)

##AUC Score
exprMatrix <-  GetAssayData(Pericytes_GBM, assay = "RNA", slot = "counts")
cells_AUC2 <-  AUCell_run(exprMatrix, IDHWT)
AUCell_auc <-  cells_AUC2@assays@data@listData[["AUC"]]
Pericytes_GBM$AUC <- AUCell_auc
common_cells <- intersect(
  rownames(pData(cds)),
  colnames(Pericytes_GBM)
)
Biobase::pData(cds)[common_cells, "AUC"] <- Pericytes_GBM$AUC[common_cells]


##179-gene signature score along pseudotime
df <- data.frame(
  cell = rownames(Biobase::pData(cds)),
  pseudotime = Biobase::pData(cds)$Pseudotime,
  AUC = Biobase::pData(cds)$AUC,
  celltype = Biobase::pData(cds)$celltype
)


##Fig.6G(G-1)
p <- ggplot(df, aes(x = pseudotime, y = AUC)) +
  geom_point(aes(color = celltype), alpha = 0.9, size = 1)+ 
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  scale_color_manual(values = cell_colors) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    panel.background = element_blank(),  
    plot.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),
  ) +labs(
    title = "179 Gene Score",
    x = "Pseudotime",
    y = "179 Gene Score",
    color = "celltype"
  )

##Fig.6F(F-1)
df <- data.frame(
  pseudotime = pData(cds)$Pseudotime,
  celltype = pData(cds)$celltype
)

p <- ggplot(df, aes(x = pseudotime, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9, color = "white") +
  scale_fill_manual(values = cell_colors) +
  # scale_y_discrete(limits = c( "Tumor","Inflamed","Normal_control"))+
  labs(x = "Pseudotime", y = "celltype", title = "") +
  theme_bw(base_size = 14) +  
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5) 
  )


##Fig.6H(H-1)
## Monocle trajectory colored by 179-gene AUC score
p <-  plot_cell_trajectory(cds, color_by = "AUC") +
  scale_color_viridis_c() +
  ggtitle("179gene score")


saveRDS(cds,"Pericyte_GBM_monocle_cds.rds") 


##-------------------02Astrocytoma-------------------------##
Pericytes_Astro    <-  subset(Pericyte, subset = group %in% c("Astrocytoma"))
data <- GetAssayData(Pericytes_Astro, assay = "RNA", slot = "counts")
data <- as(as.matrix(data),'sparseMatrix')


pd <- new('AnnotatedDataFrame', data = Pericytes_Astro@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM,cores=30)
HSMM <- detectGenes(HSMM, min_expr = 0.1 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],      
                                      fullModelFormulaStr = "~ celltype",
                                      cores = 30)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) 


HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)

table(pData(HSMM)$State,pData(HSMM)$celltype)     

cds <- orderCells(HSMM, root_state = 1) 

Pericytes_Astro$Pseudotime <- 
  Biobase::pData(cds)[colnames(Pericytes_Astro), "Pseudotime"]

p <- plot_cell_trajectory(cds,color_by="Pseudotime") 

signature <- read.xlsx("signature.xlsx",1)
Noncodel    <-  na.omit(signature$Noncodel)

##AUC Score
exprMatrix <-  GetAssayData(Pericytes_Astro, assay = "RNA", slot = "counts")
cells_AUC2 <-  AUCell_run(exprMatrix, Noncodel)
AUCell_auc <-  cells_AUC2@assays@data@listData[["AUC"]]
Pericytes_Astro$AUC <- AUCell_auc
common_cells <- intersect(
  rownames(pData(cds)),
  colnames(Pericytes_Astro)
)
Biobase::pData(cds)[common_cells, "AUC"] <- Pericytes_Astro$AUC[common_cells]


##30-gene signature score along pseudotime
df <- data.frame(
  cell = rownames(Biobase::pData(cds)),
  pseudotime = Biobase::pData(cds)$Pseudotime,
  AUC = Biobase::pData(cds)$AUC,
  celltype = Biobase::pData(cds)$celltype
)


##Fig.6G(G-2)
p <- ggplot(df, aes(x = pseudotime, y = AUC)) +
  geom_point(aes(color = celltype), alpha = 0.9, size = 1)+ 
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  scale_color_manual(values = cell_colors) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    panel.background = element_blank(),  
    plot.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),
  ) +labs(
    title = "30 Gene Score",
    x = "Pseudotime",
    y = "30 Gene Score",
    color = "celltype"
  )

##Fig.6F(F-2)
df <- data.frame(
  pseudotime = pData(cds)$Pseudotime,
  celltype = pData(cds)$celltype
)

p <- ggplot(df, aes(x = pseudotime, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9, color = "white") +
  scale_fill_manual(values = cell_colors) +
  labs(x = "Pseudotime", y = "celltype", title = "") +
  theme_bw(base_size = 14) +  
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.title = element_text(hjust = 0.5)  
  )


##Fig.6H(H-2)
## Monocle trajectory colored by 30-gene AUC score
p <-  plot_cell_trajectory(cds, color_by = "AUC") +
  scale_color_viridis_c() +
  ggtitle("30gene score")

saveRDS(cds,"Pericyte_AC_monocle_cds.rds") 



##-------------------03 Oligodendroglioma ---------------------##

Pericytes_OD     <-  subset(Pericyte, subset = group %in% c("Oligodendroglioma"))
data <- GetAssayData(Pericytes_OD , assay = "RNA", slot = "counts")
data <- as(as.matrix(data),'sparseMatrix')


pd <- new('AnnotatedDataFrame', data = Pericytes_OD @meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       # lowerDetectionLimit = 1, #默认0.1, zx0.5
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM,cores=30)
HSMM <- detectGenes(HSMM, min_expr = 0.1 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],      
                                      fullModelFormulaStr = "~ celltype",
                                      cores = 30)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) 

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)

table(pData(HSMM)$State,pData(HSMM)$celltype)
cds <- orderCells(HSMM, root_state = 1) 


Pericytes_OD$Pseudotime <- 
  Biobase::pData(cds)[colnames(Pericytes_OD), "Pseudotime"]

p <- plot_cell_trajectory(cds,color_by="Pseudotime") 

signature <- read.xlsx("signature.xlsx",1)
Codel     <-  na.omit(signature$Codel)

##AUC Score
exprMatrix <-  GetAssayData(Pericytes_OD, assay = "RNA", slot = "counts")
cells_AUC2 <-  AUCell_run(exprMatrix, Codel)
AUCell_auc <-  cells_AUC2@assays@data@listData[["AUC"]]
Pericytes_OD$AUC <- AUCell_auc
common_cells <- intersect(
  rownames(pData(cds)),
  colnames(Pericytes_OD)
)
Biobase::pData(cds)[common_cells, "AUC"] <- Pericytes_OD$AUC[common_cells]


##83-gene signature score along pseudotime
df <- data.frame(
  cell = rownames(Biobase::pData(cds)),
  pseudotime = Biobase::pData(cds)$Pseudotime,
  AUC = Biobase::pData(cds)$AUC,
  celltype = Biobase::pData(cds)$celltype
)


##Fig.6G(G-3)
p <- ggplot(df, aes(x = pseudotime, y = AUC)) +
  geom_point(aes(color = celltype), alpha = 0.9, size = 1)+ 
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "blue", size = 1) +
  scale_color_manual(values = cell_colors) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    panel.background = element_blank(),  
    plot.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),
  ) +labs(
    title = "83 Gene Score",
    x = "Pseudotime",
    y = "83 Gene Score",
    color = "celltype"
  )

##Fig.6F(F-3)
df <- data.frame(
  pseudotime = pData(cds)$Pseudotime,
  celltype = pData(cds)$celltype
)

p <- ggplot(df, aes(x = pseudotime, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9, color = "white") +
  scale_fill_manual(values = cell_colors) +
  labs(x = "Pseudotime", y = "celltype", title = "") +
  theme_bw(base_size = 14) + 
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5)  
  )


##Fig.6H(H-1)
## Monocle trajectory colored by 83-gene AUC score
p <-  plot_cell_trajectory(cds, color_by = "AUC") +
  scale_color_viridis_c() +
  ggtitle("83gene score")


saveRDS(cds,"Pericytes_OD_monocle_cds.rds") 


