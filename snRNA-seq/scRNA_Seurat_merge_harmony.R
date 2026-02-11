library(Seurat) # v5.0.3
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3) # 100 Gb

set.seed(1234)
packageVersion("Seurat")

# Load s Seurat object
lung.merge <- readRDS("lung_merge_annotation.rds")

# RNA assay
DefaultAssay(lung.merge) <- "RNA"
lung.merge

# Normalization and dimensional reduction
lung.merge <- NormalizeData(lung.merge)
lung.merge <- FindVariableFeatures(lung.merge, verbose = FALSE)
lung.merge <- ScaleData(lung.merge)
lung.merge <- RunPCA(lung.merge)

# Harmony integration
lung.merge <- IntegrateLayers(object = lung.merge, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE, npcs = 30)
lung.merge <- FindNeighbors(lung.merge, dims = 1:30, reduction = "harmony")
lung.merge <- FindClusters(lung.merge, cluster.name = "harmony_clusters", resolution = 1.2)
lung.merge <- RunUMAP(lung.merge, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

# Set colors for 6 samples
COL <- DiscretePalette(6, palette = "parade", shuffle = FALSE)
names(COL) <- c("<SAMPLE4>", "<SAMPLE6>", "<SAMPLE9>", "<SAMPLE11>", "<SAMPLE23>", "<SAMPLE28>")

# Set colors for annotation
COL2 <- c("pink", "#FF4FB5", "orange", "yellow", "green", "blue", "brown", "purple", "grey")
names(COL2) <- c("Tumor cell", "Alveolar epithelial cell", "Fibroblast", "Endothelial cell", "T cell", "B cell", "Plasma cell", "Myeloid cell", "Others")

# UMAP plot
pdf("UMAP_merge_cluster_harmony_annotation.pdf", width = 24, height = 7)
p1 <- DimPlot(lung.merge, reduction = "umap.harmony", group.by = c("original"), cols = COL)
p2 <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("harmony_clusters"))
p3 <- DimPlot(lung.merge, reduction = "umap.harmony", label = TRUE, group.by = c("annotation"), cols = COL2)
p1 + p2 + p3
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_merge_annotation_harmony.rds")

rm(list = ls())
gc();gc()
