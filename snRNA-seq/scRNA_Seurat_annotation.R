library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object
lung <- readRDS("scRNA.rds")

# Assign annotations into clusters
new.clutser.id <- c("Tumor cell", "Alveolar epithelial cell", "Tumor cell", "Fibroblast", "Tumor cell", "Tumor cell", "Endothelial cell", "Others", "Myeloid cell", "Others")
names(new.clutser.id) <- levels(lung)
lung <- RenameIdents(lung, new.clutser.id)

# Set metadata "annotation"
lung[["annotation"]] <- Idents(lung)

# Set colors for annotation
COL <- c("pink", "#FF4FB5", "orange", "yellow", "green", "blue", "brown", "purple", "grey")
names(COL) <- c("Tumor cell", "Alveolar epithelial cell", "Fibroblast", "Endothelial cell", "T cell", "B cell", "Plasma cell", "Myeloid cell", "Others")

# UMAP plot 
pdf("UMAP_cluster_annotation.pdf", width = 10, height = 7)
DimPlot(lung, reduction = "umap", group.by = "annotation", label = TRUE, label.size = 2.5, repel = TRUE, cols = COL) + guides(color = guide_legend(override.aes = list(size=4), ncol=1))
dev.off()

# Save the object
saveRDS(lung, file = "scRNA_annotation.rds")

rm(list = ls())
gc();gc()
