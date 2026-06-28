library(Seurat) # v5.0.3
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

# Set colors
COL <- c("blue", "red", "pink", "#FF4FB5", "black", "green", "purple", "brown", "orange", "yellow", "grey")
names(COL) <- c("NE tumor", "Non-NE tumor", "Tumor", "Alveolar", "Bronchiolar", "Basal", "Macrophage", "Plasma cell", "Fibroblast", "Endothelial cell", "Others")


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
    # color order
    current_cell_order <- levels(data[[i]]@idents)
    COL2 <- COL[current_cell_order]
    if (pathways.show %in% data[[i]]@netP$pathways) {
    netVisual_aggregate(data[[i]], signaling = pathways.show, layout = "circle",
    				   label.edge= F, edge.weight.max = weight.max, edge.width.max = 12,
    				   color.use = COL2)
    title(main = names(data)[i], line = 1)
    }
}
dev.off()

rm(list = ls())
gc();gc()


#############################
#Chord (TGFb to tumor cells)#
#############################
pdf("cellchat/TGFb_to_tumor_chord.pdf", width = 12, height = 16, pointsize=14)
par(mfrow = c(3,2), xpd=TRUE)
for (i in 1:length(data)) {
    #color order
    current_cell_order <- levels(data[[i]]@idents)
    COL2 <- COL[current_cell_order]
    if (pathways.show %in% data[[i]]@netP$pathways) {
    netVisual_aggregate(data[[i]], signaling = pathways.show, layout = "chord",
                                   color.use = COL2, targets.use = c("Tumor", "Non-NE tumor", "NE tumor"))
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
    g[[i]] <- netVisual_bubble(data[[i]], signaling = c("TGFb"), remove.isolate = TRUE, targets.use = c("Tumor", "Non-NE tumor", "NE tumor"),  title.name = names(data)[i])
}
}
pdf("cellchat/TGFb_bubble.pdf", width = 5, height = 5)
print (g)
dev.off()


rm(list = ls())
gc();gc()
