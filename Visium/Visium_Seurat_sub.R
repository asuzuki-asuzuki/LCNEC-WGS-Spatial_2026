library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object
lung.merge <- readRDS("lung_visium_merge_harmony_annotation.rds")

# Subset epithelial cell and tumor cell spots
lung.merge <- subset(lung.merge, subset = annotation %in% c("NE tumor", "Non-NE tumor", "Tumor", "Alveolar", "Bronchiolar", "Basal"))

# Number of spots
ncol(lung.merge@assays$Spatial$counts)

# Normalization and dimentional reduction
lung.merge <- NormalizeData(lung.merge)
lung.merge <- FindVariableFeatures(lung.merge, verbose = FALSE)
lung.merge <- ScaleData(lung.merge)
lung.merge <- RunPCA(lung.merge)

# Harmony integration
lung.merge <- IntegrateLayers(object = lung.merge, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
lung.merge <- FindNeighbors(lung.merge, dims = 1:30, reduction = "harmony")
lung.merge <- FindClusters(lung.merge, cluster.name = "harmony_clusters")
lung.merge <- RunUMAP(lung.merge, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

# Set colors for 19 samples
COL <- DiscretePalette(19, palette = "polychrome", shuffle = FALSE)

pdf("UMAP_sub_cluster_harmony.pdf", width = 16, height = 7)
p1 <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("original"), cols = COL)
p2 <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("harmony_clusters"))
p1 + p2
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_sub_merge_harmony.rds")

# Join layers
lung.merge[["Spatial"]] <- JoinLayers(lung.merge[["Spatial"]])

# DEG analysis
lung.markers <- FindAllMarkers(lung.merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_sub_merge_harmony.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)
write.table(lung.markers, "all_sub_merge_harmony.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

rm(list = ls())
gc();gc()
