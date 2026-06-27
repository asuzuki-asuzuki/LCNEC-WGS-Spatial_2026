library(Seurat) # v5.4.0
library(dplyr)
library(patchwork)
library(ggplot2)
library(CellChat)
library(tidyverse)
library(reshape2)
library(gplots)

set.seed(1234)

packageVersion("Seurat")
packageVersion("CellChat")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100 * 1024^3)

# Please read the sample list (Sample ID, Subtype)
DATA <- read.table("<LIST>", header=T)
case <- DATA$Case

data <- list()

# Load the CellChat objects
for(i in case$Case) {
      data[[i]] <- readRDS(paste0("cellchat/", i, "_cellchat.rds"))
      data[[i]] <- aggregateNet(data[[i]])
      data[[i]] <- netAnalysis_computeCentrality(data[[i]], slot.name = "netP")
}

# Merging the CellChat objects
cellchat <- mergeCellChat(data, add.names = names(data))
cellchat <- liftCellChat(cellchat, unique(unlist(cellchat@idents)))

# Save a merged object
saveRDS(cellchat, "cellchat/merge_cellchat.rds")

# Set colors for merge_annotation
COL <- DiscretePalette(14, palette = "glasbey", shuffle = FALSE)
COL2 <- c(COL[1:13], "brown", COL[14])
names(COL2) <- c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC/Pericyte", "EC", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Mastocyte", "Others")


###############
#Circle (TGFb)#
###############
pathways.show <- c("TGFb")
weight <- c()
for (i in 1:length(data)) {
    if (pathways.show %in% data[[i]]@netP$pathways) {
       weight <- c(weight, max(data[[i]]@netP$prob[,,pathways.show]))
    }
}
weight.max <- max(weight)
pdf("cellchat/TGFb.pdf")
par(mfrow = c(3,2), xpd=TRUE)
for (i in 1:length(data)) {
    #color order
    current_cell_order <- levels(data[[i]]@idents)
    COL <- COL2[current_cell_order]
    if (pathways.show %in% data[[i]]@netP$pathways) {
    netVisual_aggregate(data[[i]], signaling = pathways.show, layout = "circle",
    				   label.edge= F, edge.weight.max = weight.max, edge.width.max = 12,
    				   color.use = COL)
    title(main = names(data)[i], line = 1)
    }
}
dev.off()

rm(list = ls())
gc();gc()


################
#Circle (Notch)#
################
pathways.show <- c("NOTCH")
weight <- c()
for (i in 1:length(data)) {
    if (pathways.show %in% data[[i]]@netP$pathways) {
       weight <- c(weight, max(data[[i]]@netP$prob[,,pathways.show]))
    }
}
weight.max <- max(weight)
pdf("cellchat/NOTCH.pdf")
par(mfrow = c(3,2), xpd=TRUE)
for (i in 1:length(data)) {
    #color order
    current_cell_order <- levels(data[[i]]@idents)
    COL <- COL2[current_cell_order]
    if (pathways.show %in% data[[i]]@netP$pathways) {
    netVisual_aggregate(data[[i]], signaling = pathways.show, layout = "circle",
                                   label.edge= F, edge.weight.max = weight.max, edge.width.max = 12,
                                   color.use = COL)
    title(main = names(data)[i], line = 1)
    }
}
dev.off()


#############################
#Chord (TGFb to tumor cells)#
#############################
pdf("cellchat/TGFb_to_tumor_chord.pdf", width = 12, height = 16, pointsize=14)
par(mfrow = c(3,2), xpd=TRUE)
for (i in 1:length(data)) {
    #color order
    current_cell_order <- levels(data[[i]]@idents)
    COL <- COL2[current_cell_order]
    if (pathways.show %in% data[[i]]@netP$pathways) {
    netVisual_aggregate(data[[i]], signaling = pathways.show, layout = "chord",
                                   color.use = COL, targets.use = c("Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor"))
    title(main = names(data)[i], line = 1)
    }
}
dev.off()


##############################
#Chord (Notch to tumor cells)#
##############################
pdf("cellchat/NOTCH_to_tumor_chord.pdf", width = 12, height = 16, pointsize=14)
par(mfrow = c(3,2), xpd=TRUE)
for (i in 1:length(data)) {
    #color order
    current_cell_order <- levels(data[[i]]@idents)
    COL <- COL2[current_cell_order]
    if (pathways.show %in% data[[i]]@netP$pathways) {
    netVisual_aggregate(data[[i]], signaling = pathways.show, layout = "chord",
                                   color.use = COL, targets.use = c("Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor"))
    title(main = names(data)[i], line = 1)
    }
}
dev.off()


###############################
#Bubble plots (to tumor cells)#
###############################
g <- list()
OKcell <- c(<Sample IDs with the interaction>)
for (i in 1:length(data)) {
    if (names(data)[i] %in% OKcell){
    g[[i]] <- netVisual_bubble(data[[i]], signaling = c("TGFb"), remove.isolate = TRUE, targets.use = c("Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor"),  title.name = names(data)[i])
}
}
pdf("cellchat/TGFb_bubble.pdf", width = 5, height = 5)
print (g)
dev.off()


g <- list()
OKcell <- c(<Sample IDs with the interaction>)
for (i in 1:length(data)) {
    if (names(data)[i] %in% OKcell){
    g[[i]] <- netVisual_bubble(data[[i]], signaling = c("NOTCH"), remove.isolate = TRUE, targets.use = c("Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor"),  title.name = names(data)[i])
}
}
pdf("cellchat/NOTCH_bubble.pdf", width = 5, height = 5.5)
print (g)
dev.off()


rm(list = ls())
gc();gc()
