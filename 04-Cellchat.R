## Using IDHWT as an example
cellchat <- readRDS("IDHWT_cellchat.rds")

## CellChat output: ligand–receptor interaction table (IDHWT)
LR <- subsetCommunication(cellchat)

##IDHWT gene signature (179 gene sets)
signature <- read.xlsx("signature.xlsx",1)
signature <- na.omit(signature$IDHWT)

##ligand–receptor pairs containing signature genes
LR.filtered <- LR %>%
  dplyr::filter(ligand %in% signature | receptor %in% signature)

pair.use    <- unique(LR.filtered$interaction_name)
pair.use.df <- data.frame(interaction_name = pair.use)

p <- netVisual_bubble(cellchat,
                      pairLR.use = pair.use.df,
                      remove.isolate = TRUE)


#Pericytes: Top 10 outgoing signaling pathways

df.peri <- subsetCommunication(cellchat, sources.use = "Pericytes")
sort(unique(df.peri$pathway_name))

top.pathways <- names(sort(tapply(df.peri$prob, df.peri$pathway_name, sum), decreasing=TRUE))[1:10]

p <- netVisual_bubble(
  cellchat,
  sources.use = "Pericytes",
  signaling = top.pathways,
  remove.isolate = TRUE
)



