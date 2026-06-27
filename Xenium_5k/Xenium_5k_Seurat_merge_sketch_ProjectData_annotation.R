library(Seurat) # v5.4.0
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung.merge <- readRDS("lung_xenium_merge_sketch_ProjectData.rds")

DefaultAssay(lung.merge) <- "Xenium"
lung.merge

# Assign "merge_annotation" into clusters
Idents(object = lung.merge) <- "cluster_full"
new.clutser.id <- c("T cell", "Epithelial cell, NE tumor", "Fibroblast", "Epithelial cell, non-NE tumor", "EC", "Others", "Macrophage", "Plasma cell", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, alveolar", "SMC/Pericyte", "Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor", "Macrophage", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Macrophage", "Epithelial cell, NE tumor", "Fibroblast", "Fibroblast", "T cell", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "B cell", "Epithelial cell, NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, bronchiolar", "T cell", "Epithelial cell, non-NE tumor", "Mastocyte", "T cell")
names(new.clutser.id) <- levels(lung.merge)
lung.merge <- RenameIdents(lung.merge, new.clutser.id)

# Set metadata "merge_annotation"
lung.merge[["merge_annotation"]] <- Idents(lung.merge)

lung.merge$merge_annotation <- factor(x = lung.merge$merge_annotation, levels = c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC/Pericyte", "EC", "T cell", "B cell", "Plasma cell", "Macrophage", "Mastocyte", "Others"))

# Set colors for merge annotations
COL <- DiscretePalette(14, palette = "glasbey", shuffle = FALSE)
COL2 <- c(COL[1:13], "brown", COL[14])
names(COL2) <- c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC/Pericyte", "EC", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Mastocyte", "Others")

# UMAP plots
pdf("full_UMAP_cluster_merge_annotation.pdf", width = 10, height = 8)
p1 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("merge_annotation"), alpha = 0.2, cols = COL2)
p1
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge_sketch_ProjectData_annotation.rds")

rm(list = ls())
gc();gc()
