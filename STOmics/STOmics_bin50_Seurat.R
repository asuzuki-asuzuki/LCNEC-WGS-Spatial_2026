library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object
lung <- readRDS("bin50.rds")

# Number of spots
ncol(lung@assays$Spatial$counts)

# Normalization, dimensional reduction, clustering and UMAP
lung <- SCTransform(lung, assay = "Spatial", verbose = FALSE)
lung <- RunPCA(lung, assay = "SCT", verbose = FALSE)
lung <- FindNeighbors(lung, reduction = "pca", dims = 1:30)
lung <- FindClusters(lung, verbose = FALSE)
lung <- RunUMAP(lung, reduction = "pca", dims = 1:30)

# UMAP and spatial plots
pdf("UMAP_cluster_bin50.pdf", width = 18, height = 8)
p1 <- DimPlot(lung, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(lung, label = FALSE, pt.size.factor = 70, stroke=NA)
p1 + p2
dev.off()

# DEG analysis
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_stomics_gene50.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)
write.table(lung.markers, "all_stomics_gene50.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

# Save the object
saveRDS(lung, file = "lung_stomics_bin50.rds")

rm(list = ls())
gc();gc()
