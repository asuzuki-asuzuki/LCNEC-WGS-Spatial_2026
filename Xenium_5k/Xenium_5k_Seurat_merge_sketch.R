library(Seurat) # v5.4.0
library(ggplot2)
library(patchwork)
library(dplyr)
library(BPCells)
library(SeuratObject)

set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung.merge <- readRDS("lung_xenium_merge.rds")
lung.merge

# Normalization
DefaultAssay(lung.merge) <- "Xenium"
lung.merge <- NormalizeData(lung.merge)
lung.merge <- FindVariableFeatures(lung.merge, verbose = FALSE)
lung.merge

# Sketch data (5000 cells per sample)
lung.merge <- SketchData(object = lung.merge, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")

# Sketch assay
DefaultAssay(lung.merge) <- "sketch"
lung.merge

# Normalization, dimentional reduction, clustering and UMAP
lung.merge <- FindVariableFeatures(lung.merge, verbose = FALSE)
lung.merge <- ScaleData(lung.merge)
lung.merge[["sketch"]] <- JoinLayers(lung.merge[["sketch"]])
lung.merge <- RunPCA(lung.merge, verbose = F)
lung.merge <- FindNeighbors(lung.merge, dims = 1:30)
lung.merge <- FindClusters(lung.merge, resolution = 1.5)
lung.merge <- RunUMAP(lung.merge, dims = 1:30, return.model = T)


# Set colors for cell-type annotation in each section
lung.merge$annotation <- factor(x = lung.merge$annotation, levels = c("Tumor cell, NE", "Tumor cell, non-NE", "Tumor cell (Fetal-Ad)", "Alveolar epithelium", "Bronchiolar epithelium", "Fibroblast", "SMC", "Pericyte", "EC", "LEC", "T cell", "B cell", "Plasma cell", "NK cell", "Macrophage", "DC", "Mastocyte", "Others"))
COL <- c(DiscretePalette(14, palette = "glasbey", shuffle = FALSE), "cyan3", "darkorange", "hotpink", "cyan4")
names(COL) <- c("Tumor cell, NE", "Tumor cell, non-NE", "Basal epithelium", "Bronchiolar epithelium", "Alveolar epithelium", "Fibroblast", "SMC", "EC", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others", "DC", "Mastocyte", "Pericyte", "NK cell")

# UMAP plots
pdf("sketch_UMAP_merge_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung.merge, reduction = "umap", label = TRUE)
p2 <- DimPlot(lung.merge, reduction = "umap", group.by = c("annotation"), cols = COL)
p1 + p2
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge_sketch.rds")

rm(list = ls())
gc();gc()
