library(Seurat) # v5.0.3
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3)

set.seed(1234)
packageVersion("Seurat")

# Load a Seurat object
lung.merge <- readRDS("lung_visium_merge.rds")

# Spatial assay
DefaultAssay(lung.merge) <- "Spatial"
lung.merge

# Normalization, dimentional reduction, clustering and UMAP
lung.merge <- NormalizeData(lung.merge)
lung.merge <- FindVariableFeatures(lung.merge, verbose = FALSE)
lung.merge <- ScaleData(lung.merge)
lung.merge <- RunPCA(lung.merge)
lung.merge <- FindNeighbors(lung.merge, dims = 1:30, reduction = "pca")
lung.merge <- FindClusters(lung.merge)
lung.merge <- RunUMAP(lung.merge, dims = 1:30, reduction = "pca")

# Harmony integration
lung.merge <- IntegrateLayers(object = lung.merge, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
lung.merge <- FindNeighbors(lung.merge, dims = 1:30, reduction = "harmony")
lung.merge <- FindClusters(lung.merge, cluster.name = "harmony_clusters")
lung.merge <- RunUMAP(lung.merge, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

# Set colors for 19 samples
COL <- DiscretePalette(19, palette = "polychrome", shuffle = FALSE)

# UMAP plots
pdf("UMAP_cluster_harmony.pdf", width = 16, height = 7)
p1 <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("original"), cols = COL)
p2 <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("harmony_clusters"))
p1 + p2
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_visium_merge_harmony.rds")

rm(list = ls())
gc();gc()
