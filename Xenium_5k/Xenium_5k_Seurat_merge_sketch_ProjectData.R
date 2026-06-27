library(Seurat) # v5.4.0
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung.merge <- readRDS("lung_xenium_merge_sketch.rds")

# Projection of full data to the sketch assay
lung.merge <- ProjectData(object = lung.merge, assay = "Xenium", full.reduction = "pca.full", sketched.assay = "sketch", sketched.reduction = "pca", umap.model = "umap", dims = 1:30, refdata = list(cluster_full = "seurat_clusters"))

# Clusters
lung.merge$cluster_full <- factor(x = lung.merge$cluster_full, levels = unique(sort(lung.merge$seurat_clusters)))
table(lung.merge$cluster_full)

# Set colors for annotation in each section
COL <- c(DiscretePalette(14, palette = "glasbey", shuffle = FALSE), "cyan3", "darkorange", "hotpink", "cyan4")
names(COL) <- c("Tumor cell, NE", "Tumor cell, non-NE", "Basal epithelium", "Bronchiolar epithelium", "Alveolar epithelium", "Fibroblast", "SMC", "EC", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others", "DC", "Mastocyte", "Pericyte", "NK cell")

# UMAP plots
pdf("full_UMAP_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("cluster_full"), alpha = 0.2)
p2 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("annotation"), alpha = 0.2, cols = COL)
p1 + p2
dev.off()

# Set colors for 13 sections
big_palette <- c(
  DiscretePalette(78, palette = "parade"),
  DiscretePalette(26, palette = "alphabet"),
  DiscretePalette(26, palette = "alphabet2"),
  DiscretePalette(36, palette = "polychrome"),
  DiscretePalette(32, palette = "glasbey")
)

COL2 <- big_palette[1:length(unique(lung.merge$original))]

pdf("full_UMAP_cluster_ori.pdf", width = 9, height = 8)
p1 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("original"), alpha = 0.2, cols = COL2)
p1
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge_sketch_ProjectData.rds")

rm(list = ls())
gc();gc()
