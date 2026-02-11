library(Seurat) # v5.0.3
library(patchwork)
library(dplyr)
library(ggplot2)

set.seed(1234)
packageVersion("Seurat")

# Please input the output directory of Xenium Analyzer
DATA <- "<INPUT>"

lung <- LoadXenium(DATA, fov = "fov")
lung <- subset(lung, subset = nCount_Xenium > 0)
lung

# QC plot
pdf("Xenium_nCount.pdf")
VlnPlot(lung, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()

# Normalization, dimentional reduction, clustering and UMAP
lung <- SCTransform(lung, assay = "Xenium")
lung <- RunPCA(lung, npcs = 30, features = rownames(lung))
lung <- RunUMAP(lung, dims = 1:30)
lung <- FindNeighbors(lung, reduction = "pca", dims = 1:30)
lung <- FindClusters(lung, resolution = 0.3)

# UMAP and spatial plots
pdf("UMAP_cluster.pdf", width = 18, height = 8)
p1 <- DimPlot(lung, reduction = "umap", label = TRUE)
p2 <- ImageDimPlot(lung, size = 0.3, border.size = NA, axes = TRUE) + NoGrid()
p1 + p2
dev.off()

# Save the object
saveRDS(lung, file = "lung_xenium.rds")

## DEG analysis
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_xenium.csv",append = FALSE, quote = TRUE, sep = "," ,row.names = TRUE, col.names = TRUE)

# Save the object
saveRDS(lung.markers, file = "lung_marker.rds")

rm(list = ls())
gc();gc()
