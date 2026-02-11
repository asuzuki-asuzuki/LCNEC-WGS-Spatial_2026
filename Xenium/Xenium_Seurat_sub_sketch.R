library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3)
set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung <- readRDS("lung_xenium_merge_sketch_ProjectData_annotation.rds")

# Subset epithelial cells and tumor cells
lung <- subset(lung, subset = annotation %in% c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, alveolar", "Epithelial cell, bronchiolar", "Epithelial cell, basal"))
lung

# Sketch assay
DefaultAssay(lung) <- "sketch"
lung

# Dimentional reduction, clustering and UMAP
lung <- FindVariableFeatures(lung)
lung <- ScaleData(lung)
lung <- RunPCA(lung, verbose = F, reduction.name="sub_pca")
lung <- FindNeighbors(lung, dims = 1:20, reduction = "sub_pca")
lung <- FindClusters(lung, resolution = 0.8)
lung <- RunUMAP(lung, dims = 1:20, return.model = T, verbose = F, reduction = "sub_pca", reduction.name = "sub_umap")

# Set colors for 37 sections
COL <- c(DiscretePalette(36, palette = "polychrome", shuffle = FALSE), "#0000FF")

# UMAP plots
pdf("sketch_UMAP_sub_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung, reduction = "sub_umap", label = TRUE)
p2 <- DimPlot(lung, reduction = "sub_umap", group.by = c("original"), cols = COL)
p1 + p2
dev.off()

# Save the object
saveRDS(lung, file = "lung_sub_sketch.rds")

rm(list = ls())
gc();gc()
