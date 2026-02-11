library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3)
set.seed(1234)

packageVersion("Seurat")


# Load a Seurat object of a sample
lung <- readRDS("<SAMPLE DIR>/lung_xenium.rds")

# Load cell type annotation and subannotation for the sample
anno <- read.table("<SAMPLE>_annotation.txt"), header = T, sep = "\t", row.names=1)

lung <- subset(lung, subset = nCount_Xenium > 10)
lung

# check the number of cells
nrow(anno)
length(lung$orig.ident)

# Set annotation and subannotation to the object
lung[[c("annotation", "subannotation")]] <- anno[names(lung$orig.ident),]

# Set colors
SUB <- c("NE-1", "NE-2", "NE-3", "NE-4", "NE-5", "NE-6", "NE-7", "NE-8", "NE-9", "Non-NE-p-1", "Non-NE-p-2", "Non-NE-p-3", "Non-NE-p-4", "Non-NE-p-5", "Non-NE-y-1", "Non-NE-y-2", "Non-NE-y-3", "NE/Non-NE-p-1", "NE/Non-NE-y-1", "NE/Non-NE-y-2", "NE/Non-NE-y-3", "Alveolar-1", "Alveolar-2", "Bronchiolar-1", "Basal-1", "Fibrotic-1", "T cell-1", "T cell-2", "T cell-3", "T cell-4", "T cell-5", "T cell-6", "T cell-7", "B cell-1", "Plasma cell-1", "Plasma cell-2", "Plasma cell-3", "Plasma cell-4", "Macrophage-1", "Macrophage-2", "Macrophage-3", "Macrophage-4", "Macrophage-5", "Granulocyte-1", "Fibroblast-1", "Fibroblast-2", "Fibroblast-3", "Fibroblast-4", "Fibroblast-5", "Fibroblast-6", "Fibroblast-7", "Endothelial cell-1", "Endothelial cell-2", "Endothelial cell-3", "Endothelial cell-4", "LEC-1", "SMC-1", "Others")
lung$subannotation <- factor(x = lung$subannotation, levels = SUB)
COL2 <- c("#0000FF", "#33CCFF", "#000099", "#0099FF", "#003399", "#3366FF", "#3399FF", "#6699FF", "#0066FF", "#FF0000", "#660000", "#990000", "#CC3300", "#FF0033", "#FF6600", "#FF9900", "#FFCC00", "#663399", "#9900FF", "#660099", "#6600CC", "#FF00CC", "#FF0099", "#330033", "#00FF00", "#99FF33", "#33FFFF", "#99FFFF", "#66FFFF", "#00FFFF", "#00FCCF", "#99FCCF", "#33FF99", "#7742c0", "#109093", "#006666", "#336666", "#669999", "#FFCCFF", "#FF99FF", "#FF99CC", "#FF999F", "#FFCCCC", "#663300", "#006600", "#003300", "#336633", "#336600", "#006633", "#009933", "#339900",  "#00CCFF", "#0099FF", "#66CCFF", "#3399FF", "#a16057", "#FFD300", "#B1CC71")
names(COL2) <- SUB
okCOL2 <- COL2[sort(unique(lung$subannotation))]

ANNO <- c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC", "Endothelial cell", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others")
lung$annotation <- factor(x = lung$annotation, levels = ANNO)
COL <- DiscretePalette(14, palette = "glasbey", shuffle = FALSE)
names(COL) <- ANNO
okCOL <- COL[sort(unique(lung$annotation))]

# Spatial plots
pdf("spatialplot_annotation.pdf", width = 18, height = 8)
p1 <- ImageDimPlot(lung, size = 0.3, border.size = NA, axes = TRUE, group.by = c("annotation"), cols = okCOL) + NoGrid()
p2 <- ImageDimPlot(lung, size = 0.3, border.size = NA, axes = TRUE, group.by = c("subannotation"), cols = okCOL2) + NoGrid()
p1 + p2
dev.off()

# Save the object
saveRDS(lung, file = "lung_annotation.rds")

rm(list = ls())
gc();gc()
