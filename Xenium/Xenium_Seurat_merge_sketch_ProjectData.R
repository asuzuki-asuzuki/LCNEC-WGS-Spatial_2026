library(Seurat) #v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

options(future.globals.maxSize = 100 * 1024^3)

# Load a merged & sketched Seurat object
lung.merge <- readRDS("lung_xenium_merge_sketch.rds")
lung.merge

# Sketch assay
DefaultAssay(lung.merge) <- "sketch"

# Projection of full data to the sketch assay
lung.merge <- ProjectData(object = lung.merge, assay = "Xenium", full.reduction = "pca.full", sketched.assay = "sketch", sketched.reduction = "pca", umap.model = "umap", dims = 1:30, refdata = list(cluster_full = "seurat_clusters"))

# Xenium assay (full data)
DefaultAssay(lung.merge) <- "Xenium"
lung.merge

# 34 clusters
lung.merge$cluster_full <- factor(x = lung.merge$cluster_full, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33"))
table(lung.merge$cluster_full)

# Set colors for 37 sections
COL <- c(DiscretePalette(36, palette = "polychrome", shuffle = FALSE), "#0000FF")

# UMAP plots
pdf("full_UMAP_merge_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("cluster_full"), alpha = 0.1)
p2 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("original"), alpha = 0.1, cols = COL)
p1 + p2
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge_sketch_ProjectData.rds")

rm(list = ls())
gc();gc()
