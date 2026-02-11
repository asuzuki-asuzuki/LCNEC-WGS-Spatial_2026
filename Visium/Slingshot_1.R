library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(monocle3)
library(magrittr)
library(slingshot) # v2.10.0

set.seed(1234)

packageVersion("Seurat")
packageVersion("slingshot")

# Load a Seurat object
lung <- readRDS("lung_sub_merge_harmony_tumor.rds")

# Spatial assay
DefaultAssay(lung) <- "Spatial"
lung

# Join layers
lung[["Spatial"]] <- JoinLayers(lung[["Spatial"]])
lung

# Preparation of SCE
sce <- as.SingleCellExperiment(lung, assay = "Spatial")

# Slingshot
sce <- slingshot(sce, clusterLabels = 'harmony_clusters', reducedDim = 'UMAP.HARMONY')

# Set colors
DATA <- read.table("id_type_color.txt", header = T, comment.char = "") # please prepare a list of sample names, subtypes and colors
# Set colors for 19 samples
COL <- DiscretePalette(19, palette = "polychrome", shuffle = FALSE)
names(COL) <- sort(DATA$ID)
# Set colors for 4 ANPY subtypes
COL2 <- DATA$Type4Col
names(COL2) <- DATA$ID
# Set colors for new subtypes (NE, non-NE (ASCL2), non-NE (YAP1), Mixed, DP)
COL3 <- DATA$TypeCol
names(COL3) <- DATA$ID
# Set colors for annotation
COL4 <- c("blue", "red", "pink", "#FF4FB5", "black", "green", "purple", "brown", "orange", "yellow", "grey")
names(COL4) <- c("NE tumor", "Non-NE tumor", "Tumor", "Alveolar", "Bronchiolar", "Basal", "Macrophage", "Plasma cell", "Fibroblast", "Endothelial cell", "Others")

# Plot with annotation
clus <- lung$annotation
pdf("slingshot_harmony_tumor.pdf", width = 7, height = 8)
plot(reducedDims(sce)$UMAP.HARMONY, col = COL4[clus], cex = 0.5, pch = 16)
for (i in levels(sce@colData$annotation)) {
    text(mean(reducedDims(sce)$UMAP.HARMONY[sce@colData$annotation == i, 1]), mean(reducedDims(sce)$UMAP.HARMONY[sce@colData$annotation == i, 2]), labels = i, font = 2)
}
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# Plot with samples
sub <- lung$original
pdf("slingshot_original_tumor.pdf", width = 7, height = 8)
plot(reducedDims(sce)$UMAP.HARMONY, col = COL[sub], cex = 0.5, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# Plot with ANPT subtypes
sub <- lung$original
pdf("slingshot_ANPY_tumor.pdf", width = 7, height = 8)
plot(reducedDims(sce)$UMAP.HARMONY, col = COL2[sub], cex = 0.5, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# Plot with new subtypes
sub <- lung$original
pdf("slingshot_Type_tumor.pdf", width = 7, height = 8)
plot(reducedDims(sce)$UMAP.HARMONY, col = COL3[sub], cex = 0.5, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# Save the object
saveRDS(sce, file = "lung_slingshot_1_harmony_tumor.rds")

rm(list = ls())
gc();gc()
