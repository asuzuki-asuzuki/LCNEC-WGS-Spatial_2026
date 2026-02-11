library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 100 * 1024^3)

set.seed(1234)
packageVersion("Seurat")

# Load a Seurat object
lung.merge <- readRDS("lung_visium_merge_harmony.rds")

# Spatial assay
DefaultAssay(lung.merge) <- "Spatial"
lung.merge

# Assign annotations into clusters
new.cluster.ids <- c("NE tumor", "Tumor", "Plasma cell", "Macrophage", "Non-NE tumor", "Fibroblast", "Alveolar", "NE tumor", "Basal", "Non-NE tumor", "NE tumor", "NE tumor", "Endothelial cell", "Fibroblast", "NE tumor", "Others", "NE tumor", "Bronchiolar", "Non-NE tumor", "NE tumor", "NE tumor", "Non-NE tumor", "NE tumor")
names(new.cluster.ids) <- levels(lung.merge)
lung.merge <- RenameIdents(lung.merge, new.cluster.ids)

# Set metadata "annotation"
lung.merge[["annotation"]] <- Idents(lung.merge)

# Set colors for annotation
lung.merge$annotation <- factor(x = lung.merge$annotation, levels = c("NE tumor", "Non-NE tumor", "Tumor", "Alveolar", "Bronchiolar", "Basal", "Macrophage", "Plasma cell", "Fibroblast", "Endothelial cell", "Others"))
COL <- c("blue", "red", "pink", "#FF4FB5", "black", "green", "purple", "brown", "orange", "yellow", "grey")
names(COL) <- c("NE tumor", "Non-NE tumor", "Tumor", "Alveolar", "Bronchiolar", "Basal", "Macrophage", "Plasma cell", "Fibroblast", "Endothelial cell", "Others")

# UMAP plot 
pdf("UMAP_merge_harmony_annotation.pdf", width = 8, height = 7)
p <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("annotation"), cols = COL)
p
dev.off()

# Spatial plot
pdf("SpatialPlot_merge_harmony_annotation.pdf", width = 40, height = 28)
p <- SpatialDimPlot(lung.merge, label = FALSE, ncol = 5, group.by = c("annotation"), pt.size.factor = 1.3, image.alpha = 0, cols = COL, stroke=NA, crop=F) & NoLegend()
p
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_visium_merge_harmony_annotation.rds")

rm(list = ls())
gc();gc()
