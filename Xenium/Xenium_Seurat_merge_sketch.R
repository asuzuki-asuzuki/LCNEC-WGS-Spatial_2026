library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

options(future.globals.maxSize = 100 * 1024^3)

# Load a merged Seurat object
lung.merge <- readRDS("<DIR>/lung_xenium_merge.rds")
lung.merge

lung.merge <- subset(lung.merge, subset = nCount_Xenium > 10)
lung.merge

# Normalization
DefaultAssay(lung.merge) <- "Xenium"
lung.merge <- NormalizeData(lung.merge)
lung.merge <- FindVariableFeatures(lung.merge, verbose = FALSE)

# Sketch data (5000 cells per sample)
lung.merge <- SketchData(object = lung.merge, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")

# Sketch assay
DefaultAssay(lung.merge) <- "sketch"
lung.merge

# Normalization, dimentional reduction, clustering and UMAP
lung.merge <- FindVariableFeatures(lung.merge)
lung.merge <- ScaleData(lung.merge)
lung.merge <- RunPCA(lung.merge, verbose = F, features = rownames(lung.merge))
lung.merge <- FindNeighbors(lung.merge, dims = 1:30)
lung.merge <- FindClusters(lung.merge, resolution = 1)
lung.merge <- RunUMAP(lung.merge, dims = 1:30, min.dist = 0.1, n.neighbors = 10L, return.model = T)

# Set colors for 37 sections
COL <- c(DiscretePalette(36, palette = "polychrome", shuffle = FALSE), "#0000FF")

# UMAP plots
pdf("sketch_UMAP_merge_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung.merge, reduction = "umap", label = TRUE)
p2 <- DimPlot(lung.merge, reduction = "umap", group.by = c("original"), cols = COL)
p1 + p2
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge_sketch.rds")

rm(list = ls())
gc();gc()
